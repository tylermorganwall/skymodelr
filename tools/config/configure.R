# Prepare your package for installation here.
# Use 'define()' to define configuration variables.
# Use 'configure_file()' to substitute configuration values.

# Common: Find C/C++ compilers, deal with ccache, find the architecture
# and find CMake.
is_windows = identical(.Platform$OS.type, "windows")
is_macos = identical(Sys.info()[['sysname']], "Darwin")
is_linux = identical(Sys.info()[["sysname"]], "Linux")

# This script configures Makevars explicitly; disable the helper's
# post-script auto-configuration pass to avoid regenerating both files.
options(
  configure.common = FALSE,
  configure.platform = FALSE
)

check_atomic_libs = function() {
  # Keep this workaround restricted to Linux.
  if (!identical(Sys.info()[["sysname"]], "Linux")) {
    return("")
  }

  build_dir = tempfile("skymodelr-atomic-check-")
  if (!dir.create(build_dir)) {
    stop(
      "Cannot create temporary directory for libatomic check.",
      call. = FALSE
    )
  }

  old_wd = getwd()
  on.exit(
    {
      setwd(old_wd)
      unlink(build_dir, recursive = TRUE)
    },
    add = TRUE
  )
  setwd(build_dir)

  # Test library availability, not whether a particular atomic operation
  # happens to require an external runtime call.
  writeLines(
    'extern "C" int atomic_link_probe(void) { return 0; }',
    "atomic_link_probe.cpp"
  )

  dll = paste0("atomic_link_probe", .Platform$dynlib.ext)
  log_file = "atomic-check.log"

  status = system2(
    command = file.path(R.home("bin"), "R"),
    args = c(
      "CMD",
      "SHLIB",
      "--preclean",
      "-o",
      shQuote(dll),
      "atomic_link_probe.cpp",
      "-latomic"
    ),
    stdout = log_file,
    stderr = log_file
  )

  available = identical(status, 0L) && file.exists(dll)

  message(
    "checking whether -latomic can be linked... ",
    if (available) "yes" else "no"
  )

  if (!available) {
    message(paste(readLines(log_file, warn = FALSE), collapse = "\n"))
    return("")
  }

  "-latomic"
}

define(ATOMIC_LIBS = check_atomic_libs())

TARGET_ARCH = Sys.info()[["machine"]]
PACKAGE_BASE_DIR = normalizePath(getwd(), winslash = "/")

IMATH_INCLUDE_DIR = system.file(
  "include",
  "Imath",
  package = "libimath",
  mustWork = TRUE
)
IMATH_LIB_ARCH = normalizePath(
  sprintf(
    "%s/%s",
    system.file(
      "lib",
      package = "libimath",
      mustWork = TRUE
    ),
    Sys.info()[["machine"]]
  ),
  winslash = "/"
)

OPENEXR_INCLUDE_DIR = system.file(
  "include",
  "OpenEXR",
  package = "libopenexr",
  mustWork = TRUE
)
OPENEXR_LIB_ARCH = normalizePath(
  sprintf(
    "%s/%s",
    system.file(
      "lib",
      package = "libopenexr",
      mustWork = TRUE
    ),
    Sys.info()[["machine"]]
  ),
  winslash = "/"
)

define(
  PACKAGE_BASE_DIR = PACKAGE_BASE_DIR,
  TARGET_ARCH = TARGET_ARCH,
  IMATH_INCLUDE_DIR = IMATH_INCLUDE_DIR,
  IMATH_LIB_ARCH = IMATH_LIB_ARCH,
  OPENEXR_INCLUDE_DIR = OPENEXR_INCLUDE_DIR,
  OPENEXR_LIB_ARCH = OPENEXR_LIB_ARCH
)

if (!is_windows) {
  configure_file("src/Makevars.in")
} else {
  configure_file("src/Makevars.win.in")
}
