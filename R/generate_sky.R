#' Write a Hosek–Wilkie sky dome to EXR
#'
#' @param outfile            Default `NA`. Path to `.exr` file. If not given, the data will be returned instead.
#' @param albedo             Default `0.5`. 0–1 ground albedo.
#' @param turbidity          Default `3`. 1.7–10 atmospheric turbidity. Only valid for Hosek model.
#' @param elevation          Default `10`. Solar elevation above the horizon (degrees).
#' @param azimuth            Default `90`, directly east. Solar azimuth (degrees). Left side of the image represents North.
#' Halfway across the image represents South.
#' @param altitude 	         Default `0`. Altitude of the viewer in meters. Valid [0,15000]. Only valid for the
#' Prague model.
#' @param resolution         Default `2048`. Height of the image. Width is 2x this number.
#' @param numbercores        Default `1`. Number of threads to use in computation.
#' @param square_projection  Default `FALSE`. If \code{TRUE} use equal‑area square mapping,
#' else latitude–longitude square mapping.
#' @param hosek              Default `TRUE`. Use `"prague"` to enable the Prague 2021-22 spectral sky model.
#' @param visibility         Default `50`. Meteorological range in kilometres for Prague model.
#' @param verbose            Default `FALSE`. Whether to print progress bars/diagnostic info.
#' @param render_solar_disk  Default `TRUE`. Whether to render the solar disk in addition to the atmosphere.
#'
#' @return Either the raw data, or the data is invisibly returned if outfile is given. The EXR is written to `outfile`.
#' @export
#' @examples
#' #Generate a basic atmosphere with the Hosek model
#' temp_exr = tempfile(fileext=".exr")
#' generate_sky(temp_exr, turbidity = 3, elevation = 2, numbercores = 1)
generate_sky = function(
  outfile = NA,
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
  square_projection = FALSE,
  verbose = FALSE,
  render_solar_disk = TRUE
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
  tmp_outfile = tempfile(fileext = ".exr")
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
      if (verbose) {
        message(
          "Drawing black image as Hosek model does not produce valid output for elevation < 0."
        )
      }
      black_sky = array(0, dim = c(resolution, resolution * 2, 4))
      black_sky[,, 4] = 1
      if (!is.na(outfile)) {
        rayimage::ray_write_image(black_sky, outfile)
        return(invisible(black_sky))
      } else {
        return(black_sky)
      }
    }
  } else {
    if (elevation < -4.2) {
      if (verbose) {
        message(
          "Drawing black image as Prague model does not produce valid output for elevation < -4.2."
        )
      }
      black_sky = array(0, dim = c(resolution, resolution * 2, 4))
      black_sky[,, 4] = 1
      if (!is.na(outfile)) {
        rayimage::ray_write_image(black_sky, outfile)
        return(invisible(black_sky))
      } else {
        return(black_sky)
      }
    }
  }

  makesky_rcpp(
    outfile = tmp_outfile,
    albedo = albedo,
    turbidity = turbidity,
    elevation = elevation,
    azimuth_deg = azimuth,
    resolution = resolution,
    numbercores = numbercores,
    square_projection = square_projection,
    model = model,
    prg_dataset = coef_file,
    altitude = altitude,
    visibility = visibility,
    render_solar_disk = render_solar_disk
  )
  generated_sky = rayimage::ray_read_image(tmp_outfile)

  if (!is.na(outfile)) {
    file.copy(tmp_outfile, outfile, overwrite = TRUE)
    return(invisible(generated_sky))
  } else {
    return(generated_sky)
  }
}

#' Generate a location and time-specific sky dome (optionally with stars)
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
#' @param outfile            Default `NA`. Path to the `.exr` file to write.
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
#' @param verbose            Default `FALSE`. Whether to print progress bars/diagnostic info.
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
#' @return Either the raw data, or the data is invisibly returned if outfile is given. The EXR is written to `outfile`.
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
  outfile = NA,
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
  verbose = FALSE,
  ...
) {
  sun_altitude_azimuth = suncalc::getSunlightPosition(datetime, lat, lon)
  elevation = sun_altitude_azimuth$altitude * 180 / pi
  azimuth = 90 + sun_altitude_azimuth$azimuth * 180 / pi
  if (verbose) {
    message(sprintf(
      "Sun: %0.1f° elevation, %0.1f° azimuth",
      elevation,
      azimuth
    ))
  }
  # Just add up all three
  sky_exr = generate_sky(
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
    square_projection = square_projection,
    verbose = verbose
  )
  if (moon) {
    moon_exr = generate_moon_latlong(
      outfile = outfile,
      datetime = datetime,
      lat = lat,
      lon = lon,
      albedo = albedo,
      turbidity = turbidity,
      altitude = altitude,
      resolution = resolution,
      numbercores = numbercores,
      hosek = hosek,
      wide_spectrum = wide_spectrum,
      visibility = visibility,
      square_projection = square_projection,
      verbose = FALSE,
      ...
    )
    sky_exr = sky_exr + moon_exr
  }
  if (stars) {
    stars_exr = generate_stars(
      resolution = resolution,
      lon_deg = lon,
      lat_deg = lat,
      time_utc = as.POSIXct(datetime),
      turbidity = turbidity,
      altitude = altitude,
      numbercores = numbercores,
      ...
    )
    sky_exr = sky_exr + stars_exr
  }
  sky_exr[,, 4] = 1
  if (!is.na(outfile)) {
    rayimage::ray_write_image(sky_exr, outfile, clamp = FALSE)
    return(invisible(sky_exr))
  } else {
    return(sky_exr)
  }
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
