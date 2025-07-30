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

#' Generate a location‑ and time‑specific sky dome
#'
#' Convenience wrapper around [generate_sky()] that first determines the Sun’s
#' apparent position for the supplied date‑time and geographic coordinates
#' (using **suncalc**), then renders the corresponding Hosek–Wilkie or Prague
#' spectral sky model to an OpenEXR file.  Call this when you would rather
#' think in *where* and *when* than in raw solar angles.
#'
#' @param outfile            Default `""`. Path to the `.exr` file to write.
#' @param datetime           Default `"2025-07-29 18:00:00"`. POSIX‑compatible
#'   date‑time string (e.g. YYYY-MM-DD HH:mm::ss) or [POSIXct] object.  If the string lacks a time‑zone
#'   suffix it is interpreted in the system time‑zone.
#' @param lat                Observer latitude in decimal
#'   degrees; positive north of the equator.
#' @param lon                Observer longitude in decimal
#'   degrees; positive east of Greenwich, negative west.
#' @param albedo             Default `0.5`. Ground albedo fraction, range [0, 1].
#' @param turbidity          Default `3`. Atmospheric turbidity, range
#'   [1.7, 10].  Hosek model only.
#' @param altitude           Default `0`. Observer altitude in metres, range
#'   [0, 15000].  Prague model only.
#' @param resolution         Default `2048`. Height of the output image in
#'   pixels; width is always twice this value.
#' @param numbercores        Default `1`. Number of CPU threads to use.
#' @param hosek              Default `TRUE`. Set to `FALSE` to use the Prague
#'   2021‑22 spectral sky model.
#' @param wide_spectrum      Default `FALSE`. Use the 55‑channel, polarised
#'   Prague coefficients.  Ignored for Hosek.
#' @param visibility         Default `50`. Meteorological range in kilometres;
#'   relevant only for the Prague model.
#' @param square_projection  Default `FALSE`. If `TRUE`, write an equal‑area
#'   square projection; otherwise write a latitude–longitude (equirectangular)
#'   projection.
#'
#' @details The Sun’s altitude (degrees above the horizon) and azimuth (degrees clockwise
#'   from east, such that 90° = south) are calculated internally from
#'   `datetime`, `lat`, and `lon`.  This keeps the high‑level interface clean
#'   while delegating heavy lifting to [generate_sky()].
#'
#' @return Invisible `NULL`; the EXR file is created at `outfile`.
#' @export
#'
#' @examples
#' # Render the evening sky over Washington, DC on 30 July 2025
#' temp_exr = tempfile(fileext = ".exr")
#' generate_sky_latlong(
#'   outfile      = temp_exr,
#'   datetime     = "2025-07-30 20:20:00",
#'   lat          = 38.9072,
#'   lon          = -77.0369,
#'   turbidity    = 3,
#'   numbercores  = 2
#' )
generate_sky_latlong = function(
  outfile,
  datetime = "2025-07-29 18:00:00",
  lat = 38.9072,
  lon = -77.0369,
  albedo = 0.5,
  turbidity = 3,
  altitude = 0,
  resolution = 2048,
  numbercores = 1,
  hosek = TRUE,
  wide_spectrum = FALSE,
  visibility = 50,
  square_projection = FALSE
) {
  sun_altitude_azimuth = suncalc::getSunlightPosition(datetime, lat, lon)
  elevation = sun_altitude_azimuth$altitude * 180 / pi
  azimuth = 90 + sun_altitude_azimuth$azimuth * 180 / pi
  generate_sky(
    outfile = outfile,
    albedo = albedo,
    turbidity = turbidity,
    altitude = altitude,
    elevation = elevation,
    azimuth = azimuth,
    resolution = resolution,
    numbercores = numbercores,
    hosek = hosek,
    wide_spectrum = wide_spectrum,
    visibility = visibility,
    square_projection = square_projection
  )
}

#' Sample a direction from the Prague model.
#'
#' @param phi                Horizontal angle of the sample, degrees. Vectorized. Range [0, 360].
#' @param theta              Vertical angle of the sample, degrees. Vectorized. Range [-90, 90].
#' @param altitude 	         Default `0`, vectorized. Altitude of the viewer in meters. Range [0, 15000].
#' @param elevation          Default `10`, vectorized. Solar elevation angle above/below the horizon (°). Range [-4.2, 90].
#' @param visibility         Default `50`, vectorized. Range [20, 131.8]. Meteorological range in kilometers for Prague model.
#' @param albedo             Default `0.5`, vectorized. Range [0, 1]. Ground albedo.
#' @param azimuth            Default `90`, single value. Solar azimuth (°). Defaults South.
#' @param numbercores        Default `1`. Number of threads to use in computation.
#' @param wide_spectrum      Default `FALSE`. Whether to use the wide‑spectrum (55‑channel, polarised) coefficients.
#'
#' @return 3 column RGB matrix.
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
  numbercores = 1,
  wide_spectrum = FALSE
) {
  sea_level = all(altitude == 0)
  stopifnot(all(phi <= 360 & phi >= 0))
  stopifnot(all(theta <= 90 & theta >= -90))
  stopifnot(all(altitude <= 15000 & altitude >= 0))
  stopifnot(all(visibility <= 15000 & altitude >= 0))
  filesize = ""
  df_values = as.list(data.frame(
    phi = phi,
    theta = theta,
    altitude = altitude,
    elevation = elevation,
    visibility = visibility,
    albedo = albedo
  ))
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
    df_values$phi,
    df_values$theta,
    df_values$elevation,
    df_values$albedo,
    df_values$altitude,
    df_values$visibility,
    azimuth[1],
    numbercores,
    coef_file
  ))
}
