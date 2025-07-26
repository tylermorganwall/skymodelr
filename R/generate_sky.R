#' Write a Hosek–Wilkie sky dome to EXR
#'
#' @param outfile            Default `""`. Path to `.exr` file.
#' @param albedo             Default `0.5`. 0–1 ground albedo.
#' @param turbidity          Default `3`. 1.7–10 atmospheric turbidity. Only valid for Hosek model.
#' @param elevation          Default `10`. Solar elevation above the horizon (°).
#' @param azimuth            Default `90`. Solar azimuth (°). Defaults South.
#' @param altitude 	         Default `0`. Altitude of the viewer in meters. Valid [0,15000]. Only valid for the
#' Prague model.
#' @param resolution         Default `2048`. Height of the image. Width is 2x this number.
#' @param numbercores        Default `1`. Number of threads to use in computation.
#' @param square_projection  Default `FALSE`. If \code{TRUE} use equal‑area square mapping,
#' else latitude–longitude square mapping.
#' @param hosek              Default `TRUE`. Use `"prague"` to enable the Prague 2021-22 spectral sky model.
#' @param prg_dataset        Default `""`. Full path to the Prague binary dataset (`*.dat`) when `model = "prague"`.
#' @param visibility         Default `50`. Meteorological range in kilometres for Prague model.
#'
#' @return Invisible `NULL`.  The EXR is written to `outfile`.
#' @export
#' @examples
#' #Generate a basic atmosphere with the Hosek model
#' temp_exr = tempfile(fileext=".exr")
#' generate_sky(temp_exr, turbidity = 3, elevation = 2, numbercores = 1)
generate_sky = function(
  outfile,
  albedo = 0.5,
  turbidity = 3,
  elevation = 10,
  azimuth = 90,
  altitude = 0,
  resolution = 2048,
  numbercores = 1,
  hosek = TRUE,
  wide_spectrum = FALSE,
  visibility = 50,
  square_projection = FALSE
) {
  sea_level = altitude == 0
  filesize = ""
  if (sea_level & !wide_spectrum) {
    filesize = "107MB"
  } else if (sea_level & wide_spectrum) {
    filesize = "574MB"
  } else if (!sea_level & !wide_spectrum) {
    filesize = "2.4GB"
  }
  check_coef_file = function(filename) {
    coef_file = file.path(tools::R_user_dir("skymodelr", "data"), filename)
    if (!file.exists(coef_file)) {
      response = readline(
        prompt = sprintf(
          " Coefficient file for this setting not yet present: this is a large file (%s), download? [y/n] ",
          filesize
        )
      )
      if (response == "y") {
        download_sky_data(sea_level, wide_spectrum)
      } else if (response == "n") {
        return("")
      } else {
        stop("Input not recognized.")
      }
    }
    return(coef_file)
  }
  coef_file = ""
  if (!hosek) {
    stopifnot(altitude >= 0 && altitude <= 15000)
    model = "prague"
    if (!wide_spectrum) {
      if (altitude == 0) {
        coef_file = check_coef_file("SkyModelDatasetGround.dat")
      } else {
        coef_file = check_coef_file("SkyModelDataset.dat")
      }
    } else {
      if (altitude == 0) {
        coef_file = check_coef_file("PragueSkyModelDatasetGroundInfra.dat")
      } else {
        stop("`wide_spectrum = TRUE` is only valid when `altitude == 0`.")
      }
    }
    if (coef_file == "") {
      stop("No coefficient file downloaded for this set of inputs.")
    }
  } else {
    model = "hosek"
  }
  makesky_rcpp(
    outfile,
    albedo,
    turbidity,
    elevation,
    azimuth,
    resolution,
    numbercores,
    square_projection,
    visibility = visibility,
    prg_dataset = coef_file,
    altitude,
    model = model
  )
  invisible(NULL)
}

#' Write a Hosek–Wilkie sky dome to EXR
#'
#' @param outfile            Default `""`. Path to `.exr` file.
#' @param albedo             Default `0.5`. 0–1 ground albedo.
#' @param turbidity          Default `3`. 1.7–10 atmospheric turbidity. Only valid for Hosek model.
#' @param elevation          Default `10`. Solar elevation above the horizon (°).
#' @param azimuth            Default `90`. Solar azimuth (°). Defaults South.
#' @param altitude 	         Default `0`. Altitude of the viewer in meters. Valid [0,15000]. Only valid for the
#' Prague model.
#' @param resolution         Default `2048`. Height of the image. Width is 2x this number.
#' @param numbercores        Default `1`. Number of threads to use in computation.
#' @param square_projection  Default `FALSE`. If \code{TRUE} use equal‑area square mapping,
#' else latitude–longitude square mapping.
#' @param hosek              Default `TRUE`. Use `"prague"` to enable the Prague 2021-22 spectral sky model.
#' @param prg_dataset        Default `""`. Full path to the Prague binary dataset (`*.dat`) when `model = "prague"`.
#' @param visibility         Default `50`. Meteorological range in kilometres for Prague model.
#'
#' @return Invisible `NULL`.  The EXR is written to `outfile`.
#' @export
#' @examples
#' #Generate a basic atmosphere with the Hosek model
#' temp_exr = tempfile(fileext=".exr")
#' generate_sky(temp_exr, turbidity = 3, elevation = 2, numbercores = 1)
calculate_sky_values = function(
  phi,
  theta,
  altitude = 0,
  elevation = 10,
  visibility = 50,
  albedo = 0.5,
  azimuth = 90,
  resolution = 2048,
  numbercores = 1,
  wide_spectrum = FALSE
) {
  sea_level = all(altitude == 0)
  filesize = ""
  if (sea_level & !wide_spectrum) {
    filesize = "107MB"
  } else if (sea_level & wide_spectrum) {
    filesize = "574MB"
  } else if (!sea_level & !wide_spectrum) {
    filesize = "2.4GB"
  }
  check_coef_file = function(filename) {
    coef_file = file.path(tools::R_user_dir("skymodelr", "data"), filename)
    if (!file.exists(coef_file)) {
      response = readline(
        prompt = sprintf(
          " Coefficient file for this setting not yet present: this is a large file (%s), download? [y/n] ",
          filesize
        )
      )
      if (response == "y") {
        download_sky_data(sea_level, wide_spectrum)
      } else if (response == "n") {
        return("")
      } else {
        stop("Input not recognized.")
      }
    }
    return(coef_file)
  }
  coef_file = ""
  stopifnot(all(altitude >= 0 & altitude <= 15000))
  model = "prague"
  if (!wide_spectrum) {
    if (sea_level) {
      coef_file = check_coef_file("SkyModelDatasetGround.dat")
    } else {
      coef_file = check_coef_file("SkyModelDataset.dat")
    }
  } else {
    if (sea_level) {
      coef_file = check_coef_file("PragueSkyModelDatasetGroundInfra.dat")
    } else {
      stop("`wide_spectrum = TRUE` is only valid when `altitude == 0`.")
    }
  }
  if (coef_file == "") {
    stop("No coefficient file downloaded for this set of inputs.")
  }
  return(calculate_raw_prague(
    phi,
    theta,
    elevation,
    albedo,
    altitude,
    visibility,
    azimuth,
    numbercores,
    coef_file
  ))
}
