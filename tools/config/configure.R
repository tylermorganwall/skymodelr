# Prepare your package for installation here.
# Use 'define()' to define configuration variables.
# Use 'configure_file()' to substitute configuration values.

# Common: Find C/C++ compilers, deal with ccache, find the architecture
# and find CMake.
is_windows = identical(.Platform$OS.type, "windows")
is_macos = identical(Sys.info()[['sysname']], "Darwin")

# This script configures Makevars explicitly; disable the helper's
# post-script auto-configuration pass to avoid regenerating both files.
options(
	configure.common = FALSE,
	configure.platform = FALSE
)

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
