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
#' @param sea_level Default `TRUE`. Download the sea‑level–only data.
#' Set to `FALSE` to download the full‑altitude dataset.
#' @param wide_spectrum Default `FALSE`. If `TRUE`, downloads the
#' wide‑spectrum (55‑channel, polarised) version. Valid only when `sea_level = TRUE`.
#'
#' @returns Invisibly, the full path to the data file.
#' @export
#'
#' @examples
#' \dontrun{
#'   # Standard (11‑channel, sea‑level) coefficients
#'   download_sky_data()
#'
#'   # Wide‑spectrum sea‑level coefficients
#'   download_sky_data(wide_spectrum = TRUE)
#'
#'   # Full altitude‑range coefficients
#'   download_sky_data(sea_level = FALSE)
#' }
download_sky_data = function(sea_level = TRUE, wide_spectrum = FALSE) {
  if (wide_spectrum && !sea_level) {
    stop("`wide_spectrum = TRUE` is only valid when `sea_level = TRUE`.")
  }

  file_info = if (sea_level) {
    if (wide_spectrum) {
      list(
        name = "PragueSkyModelDatasetGroundInfra.dat",
        url = "https://skydata.tylermw.com/PragueSkyModelDatasetGroundInfra.dat",
        md5 = "33b97729ce8cd7fbdfc60317ff2805e7"
      )
    } else {
      list(
        name = "SkyModelDatasetGround.dat",
        url = "https://skydata.tylermw.com/SkyModelDatasetGround.dat",
        md5 = "953864f62c9434269b7f909cdff6c7cd"
      )
    }
  } else {
    list(
      name = "SkyModelDataset.dat",
      url = "https://skydata.tylermw.com/SkyModelDataset.dat",
      md5 = "dcfa2bd165cf803c1791c469e05d134c"
    )
  }

  cache_dir = tools::R_user_dir("skymodelr", "data")
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }

  tmp_path = tempfile()
  on.exit(if (file.exists(tmp_path)) unlink(tmp_path), add = TRUE)
  final_dest_path = file.path(cache_dir, file_info$name)

  file_present = FALSE
  if (file.exists(final_dest_path)) {
    message("File exists, comparing checksum.")
    if (tools::md5sum(final_dest_path) == file_info$md5) {
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
    success = utils::download.file(
      file_info$url,
      destfile = tmp_path,
      mode = "wb",
      quiet = FALSE
    )

    if (success == 0) {
      if (!file.rename(tmp_path, final_dest_path)) {
        stop("Failed to move file into place.")
      }
      message("Saved to ", final_dest_path)
    } else {
      message("Download failed.")
    }
  }
  invisible(final_dest_path)
}
