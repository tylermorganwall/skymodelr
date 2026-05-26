#' Download Prague Sky Model Coefficient Data
#'
#' This model which allows for sun angles below the horizon,
#' wide spectral ranges including infrared, polarization, as well
#' rendering at altitude) need to download a coefficient file. There are three versions
#' of this dataset, listed here in order of increasing size:
#'
#' | Argument combination               | File                                 | Size |
#' |------------------------------------|--------------------------------------|----------------|
#' | `sea_level = TRUE,  wide_spectrum = FALSE` | `SkyModelDatasetGround.dat`            |   107MB |
#' | `sea_level = TRUE,  wide_spectrum = TRUE`  | `PragueSkyModelDatasetGroundInfra.dat` |   574MB |
#' | `sea_level = FALSE`                          | `SkyModelDataset.dat`                  | 2.4GB |
#'
#'
#' @param sea_level Default `TRUE`. Download the sea-level-only data.
#' Set to `FALSE` to download the full-altitude dataset.
#' @param wide_spectrum Default `FALSE`. If `TRUE`, downloads the
#' wide-spectrum (55-channel, polarised) version. Valid only when `sea_level = TRUE`.
#'
#' @returns Invisibly, the full path to the data file.
#' @export
#'
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' # Standard (11-channel, sea-level) coefficients
#' download_sky_data()
#'
#' # Wide-spectrum sea-level coefficients
#' download_sky_data(wide_spectrum = TRUE)
#'
#' # Full altitude-range coefficients
#' download_sky_data(sea_level = FALSE)
download_sky_data = function(sea_level = TRUE, wide_spectrum = FALSE) {
  if (wide_spectrum && !sea_level) {
    stop("`wide_spectrum = TRUE` is only valid when `sea_level = TRUE`.")
  }

  file_info = prague_coef_info(
    sea_level = sea_level,
    wide_spectrum = wide_spectrum
  )

  cache_dir = prague_coef_cache_dir()
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }

  tmp_path = tempfile(pattern = paste0(file_info$name, "."), tmpdir = cache_dir)
  on.exit(if (file.exists(tmp_path)) unlink(tmp_path), add = TRUE)
  final_dest_path = file.path(cache_dir, file_info$name)

  file_present = FALSE
  if (file.exists(final_dest_path)) {
    message("File exists, comparing checksum.")
    if (identical(unname(tools::md5sum(final_dest_path)), file_info$md5)) {
      file_present = TRUE
    } else {
      message(
        "File exists in cache folder but checksum not correct: removing ",
        final_dest_path
      )
      file.remove(final_dest_path)
    }
  }

  if (file_present) {
    message(
      "File already present at '",
      final_dest_path,
      "' with correct checksum."
    )
  } else {
    message("Downloading ", file_info$name)
    old_timeout = getOption("timeout")
    on.exit(options(timeout = old_timeout), add = TRUE)
    options(timeout = 3600)
    download_error = new.env(parent = emptyenv())
    status = tryCatch(
      utils::download.file(
        file_info$url,
        destfile = tmp_path,
        mode = "wb",
        quiet = FALSE
      ),
      error = function(e) {
        download_error$message = conditionMessage(e)
        NA_integer_
      }
    )

    if (!identical(status, 0L) && !identical(status, 0)) {
      if (file.exists(tmp_path)) {
        unlink(tmp_path)
      }
      detail = if (
        exists("message", envir = download_error, inherits = FALSE)
      ) {
        paste0(": ", download_error$message)
      } else {
        paste0(" with status ", status)
      }
      stop(
        "Failed to download Prague coefficient file '",
        file_info$name,
        "'",
        detail,
        ".",
        call. = FALSE
      )
    }

    if (!file.exists(tmp_path)) {
      stop(
        "Download reported success but no temporary file was created for '",
        file_info$name,
        "'.",
        call. = FALSE
      )
    }

    if (!identical(unname(tools::md5sum(tmp_path)), file_info$md5)) {
      unlink(tmp_path)
      stop(
        "Checksum mismatch after downloading Prague coefficient file '",
        file_info$name,
        "'. The partial file was deleted.",
        call. = FALSE
      )
    }

    if (file.exists(final_dest_path)) {
      unlink(final_dest_path)
    }
    moved = file.rename(tmp_path, final_dest_path)
    if (!moved) {
      moved = file.copy(tmp_path, final_dest_path, overwrite = TRUE)
      unlink(tmp_path)
    }
    if (!moved) {
      stop("Failed to move Prague coefficient file into place.", call. = FALSE)
    }

    if (!identical(unname(tools::md5sum(final_dest_path)), file_info$md5)) {
      unlink(final_dest_path)
      stop(
        "Checksum mismatch after installing Prague coefficient file '",
        file_info$name,
        "'. The bad file was deleted.",
        call. = FALSE
      )
    }
    message("Saved to ", final_dest_path)
  }
  invisible(final_dest_path)
}
