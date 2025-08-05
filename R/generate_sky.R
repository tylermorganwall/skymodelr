#' Write a Hosek–Wilkie sky dome to EXR
#'
#' @param outfile            Default `""`. Path to `.exr` file.
#' @param albedo             Default `0.5`. 0–1 ground albedo.
#' @param turbidity          Default `3`. 1.7–10 atmospheric turbidity. Only valid for Hosek model.
#' @param elevation          Default `10`. Solar elevation above the horizon (degrees).
#' @param azimuth            Default `90`. Solar azimuth (degrees). Defaults South.
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
  #If sun below valid region, write black image
  if (model == "hosek") {
    if (elevation < 0.0) {
      message(
        "Drawing black image as Hosek model does not produce valid output for elevation < 0."
      )
      black_sky = array(0, dim = c(resolution, resolution * 2, 4))
      black_sky[,, 4] = 1
      rayimage::ray_write_image(black_sky, outfile)
      return(invisible(NULL))
    }
  } else {
    if (elevation < -4.2) {
      message(
        "Drawing black image as Prague model does not produce valid output for elevation < -4.2."
      )

      black_sky = array(0, dim = c(resolution, resolution * 2, 4))
      black_sky[,, 4] = 1
      rayimage::ray_write_image(black_sky, outfile)
      return(invisible(NULL))
    }
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
  return(invisible(NULL))
}

#' Generate a location- and time-specific sky dome (optionally with stars)
#'
#' Convenience wrapper around [generate_sky()] that:
#' 1. Computes the Sun’s apparent position for `datetime`, `lat`, and `lon`
#'    (via **suncalc**).
#' 2. Renders the corresponding sky model.
#' 3. Optionally overlays a star field using [generate_stars()].
#'
#' A “black-sky” optimisation is applied when *all* of the following hold:
#' * Sun elevation < 0 degrees (−4.2 degrees for Prague model), and
#' * `stars = TRUE`.
#'
#' In that case the sky contribution is effectively zero, so only the star
#' field is written.
#'
#' @param outfile            Default `""`. Path to the `.exr` file to write.
#' @param datetime           Default `"2025-07-29 18:00:00"`. POSIX-compatible
#'   date-time; if missing a time-zone it is assumed to be local.
#' @param lat                Default `38.9072`. Observer latitude (degrees N).
#' @param lon                Default `-77.0369`. Observer longitude (degrees E; west < 0).
#' @param albedo             Default `0.5`. Ground albedo, range [0, 1].
#' @param turbidity          Default `3`. Atmospheric turbidity, range
#'   [1.7, 10] (*Hosek only*).
#' @param altitude           Default `0`. Observer altitude (m), range
#'   [0, 15000] (*Prague only*).
#' @param resolution         Default `2048`. Image height in pixels (width = 2 × height).
#' @param numbercores        Default `1`. CPU threads to use.
#' @param hosek              Default `TRUE`. `FALSE` selects the Prague model.
#' @param wide_spectrum      Default `FALSE`. 55-channel Prague coefficients (altitude = 0 m only).
#' @param visibility         Default `50`. Meteorological range (km); *Prague only*.
#' @param square_projection  Default `FALSE`. `TRUE` results in an equal-area square mapping.
#' @param stars              Default `FALSE`. If `TRUE`, composite a star field
#'   using [generate_stars()].  Automatically falls back to a pure star render
#'   when the black-sky condition is met (see *Details*).
#' @param ...                Additional arguments passed to [generate_stars()]
#'
#' @details
#' *Solar angles* — altitude (degrees above the horizon) and azimuth (degrees clockwise from
#' east, so 90 degrees = south) — are derived internally; you never have to supply
#' them directly.
#'
#' *Black-sky rule* — With the Prague model the sky radiance is defined only
#' down to −4.2 degrees, and with the Hosek model it is defined only about 0 degrees.
#' Below that the function skips the sky render and writes **only stars** when `stars = TRUE`.
#'
#' @return Invisible `NULL`; the EXR file is created at `outfile`.
#' @export
#'
#' @examples
#' # Evening twilight with stars over Washington, DC
#' temp_exr = tempfile(fileext = ".exr")
#' generate_sky_latlong(
#'   outfile     = temp_exr,
#'   datetime    = "2025-07-30 20:20:00",
#'   lat         = 38.9072,
#'   lon         = -77.0369,
#'   turbidity   = 3,
#'   stars       = TRUE,
#'   numbercores = 2
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
  square_projection = FALSE,
  stars = FALSE,
  moon = FALSE,
  ...
) {
  tmp_sky = tempfile(fileext = ".exr")
  sun_altitude_azimuth = suncalc::getSunlightPosition(datetime, lat, lon)
  elevation = sun_altitude_azimuth$altitude * 180 / pi
  azimuth = 90 + sun_altitude_azimuth$azimuth * 180 / pi

  if (moon) {
    moon_altitude_azimuth = suncalc::getMoonPosition(datetime, lat, lon)
    moon_elevation = moon_altitude_azimuth$altitude * 180 / pi
    moon_azimuth = 90 + moon_altitude_azimuth$azimuth * 180 / pi
    moon_phase_info = suncalc::getMoonIllumination(datetime) # 0.5 = full
    moon_phase = ((moon_phase_info$phase - 0.5) * 2 * pi) * 180 / pi
    #sun is 10,000 footcandles/100,760 lux at noon
    sun_brightness_ftcndl = 10760
    m = -12.73 + 0.026 * abs(moon_phase) + 4 * 10^(-9) * moon_phase^4 # m = -12.7 at angle == 0, which is correct
    opposition_surge_multiplier = max(1.35 - 0.05 * moon_phase, 1) #1.35 at max, 1.0 for continuity, just lerping
    moon_brightness_ftcndl = opposition_surge_multiplier *
      10^(-0.4 * (m + 16.57))

    scale_sun_to_moon = moon_brightness_ftcndl / sun_brightness_ftcndl

    tmp_moon = tempfile(fileext = ".exr")
    generate_sky(
      outfile = tmp_moon,
      albedo = albedo,
      turbidity = turbidity,
      altitude = altitude,
      elevation = moon_elevation,
      azimuth = moon_azimuth,
      resolution = resolution,
      numbercores = numbercores,
      hosek = hosek,
      wide_spectrum = wide_spectrum,
      visibility = visibility,
      square_projection = square_projection
    )
    moon_exr = rayimage::ray_read_image(tmp_moon)
    moon_exr[,, 1:3] = moon_exr[,, 1:3] * scale_sun_to_moon
    rayimage::ray_write_image(moon_exr, tmp_moon)
  }
  # Just add up all three
  generate_sky(
    outfile = tmp_sky,
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
  if (moon) {
    tmp_sky_exr = rayimage::ray_read_image(tmp_sky)
    tmp_moon_exr = rayimage::ray_read_image(tmp_moon)
    tmp_combined_sky_moon = tmp_sky_exr + tmp_moon_exr
    tmp_combined_sky_moon[,, 4] = 1
    rayimage::ray_write_image(tmp_combined_sky_moon, tmp_sky)
  }
  # }
  if (stars) {
    tmp_stars = tempfile(fileext = ".exr")
    generate_stars(
      outfile = tmp_stars,
      resolution = resolution,
      lon_deg = lat,
      lat_deg = lon,
      time_utc = as.POSIXct(datetime),
      turbidity = turbidity,
      altitude = altitude,
      numbercores = numbercores,
      ...
    )
    sky_exr = rayimage::ray_read_image(tmp_sky)
    stars_exr = rayimage::ray_read_image(tmp_stars)
    tmp_stars_sky = sky_exr + stars_exr
    tmp_stars_sky[,, 4] = 1
    rayimage::ray_write_image(tmp_stars_sky, filename = tmp_sky)
  }
  file.copy(tmp_sky, outfile, overwrite = TRUE)
  return(invisible())
}

#' Sample a direction from the Prague model.
#'
#' @param phi                Horizontal angle of the sample, degrees. Vectorized. Range [0, 360].
#' @param theta              Vertical angle of the sample, degrees. Vectorized. Range [-90, 90].
#' @param altitude 	         Default `0`, vectorized. Altitude of the viewer in meters. Range [0, 15000].
#' @param elevation          Default `10`, vectorized. Solar elevation angle above/below the horizon (degrees). Range [-4.2, 90].
#' @param visibility         Default `50`, vectorized. Range [20, 131.8]. Meteorological range in kilometers for Prague model.
#' @param albedo             Default `0.5`, vectorized. Range [0, 1]. Ground albedo.
#' @param azimuth            Default `90`, single value. Solar azimuth (degrees). Defaults South.
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
