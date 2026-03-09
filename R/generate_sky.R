normalize_render_mode = function(render_mode) {
  if (is.logical(render_mode)) {
    stopifnot(length(render_mode) == 1)
    return(if (render_mode) "all" else "atmosphere")
  }
  render_mode = tolower(as.character(render_mode))
  match.arg(render_mode, c("all", "atmosphere", "sun"))
}

#' Generate a Hosek-Wilkie sky dome array
#'
#' @description Evaluate either the Hosek-Wilkie or Prague analytic sky models
#' and return a high-dynamic-range image array for the given solar
#' configuration. An image file is written only when `filename` is supplied.
#'
#' @param filename           Default `NA`. Path to an image file to write. If not given, the array is returned instead.
#' @param albedo             Default `0.1`. 0.0-1.0 ground albedo. Grass has an albedo of about 0.09, while a landscape covered in
#' snow will have an albedo of 1.0.
#' @param turbidity          Default `3`. 1.7-10 atmospheric turbidity. Only valid for Hosek model.
#' @param elevation          Default `10`. Solar elevation above the horizon (degrees).
#' @param azimuth            Default `90`, sun directly east. Solar azimuth (degrees). The left edge of the image faces north and the middle faces south.
#' @param altitude           Default `0`. Altitude of the viewer in meters. Valid range: 0 to 15000. Only valid for the Prague model.
#' @param resolution         Default `2048`. Height of the image. Width is twice this number.
#' @param numbercores        Default `1`. Number of threads to use in computation.
#' @param hosek              Default `TRUE`. Set to `FALSE` to enable the Prague 2021-22 spectral sky model.
#' @param wide_spectrum      Default `FALSE`. Use the 55-channel Prague coefficients (sea level only).
#' @param visibility         Default `50`. Meteorological range in kilometres for Prague model.
#' @param verbose            Default `FALSE`. Whether to print progress bars/diagnostic info.
#' @param render_mode        Default `"all"`. One of `"all"`, `"atmosphere"`, or `"sun"`.
#'   Use `"all"` for atmosphere + solar disk, `"atmosphere"` for atmospheric radiance only, or `"sun"` for the solar disk only.
#' @param below_horizon      Default `TRUE`. Whether to sample atmospheric scattering below the horizon, which is non-zero when altitude > 0.
#'
#' @return Either the image array, or the array is invisibly returned if a file
#'   is written. The array has dimensions `(resolution, 2 * resolution, 4)`.
#' @note Writing to non-EXR formats will introduce precision loss because HDR
#'   data are quantised to the destination format, and low dynamic range outputs like PNG
#'   and JPEG files will not represent the true luminosity values encoded in the array.
#' @export
#' @examples
#' if(run_documentation()) {
#' # Hosek model (default): clear morning, Sun SE, with solar disk
#' generate_sky(
#'   resolution = 400,
#'   elevation  = 15,
#'   azimuth    = 135,
#'   turbidity  = 3,
#'   render_mode = "all"
#' ) |>
#'   rayimage::plot_image()
#' }
#' if(run_documentation()) {
#' # Same view but hazier and without the solar disk
#' generate_sky(
#'   resolution = 400,
#'   elevation  = 15,
#'   azimuth    = 135,
#'   turbidity  = 6,
#'   render_mode = "atmosphere"
#' ) |>
#'   rayimage::plot_image()
#' }
#' # Prague model (may prompt to download coefficients on first use)
#' \donttest{
#' if(run_documentation()) {
#' generate_sky(
#'   resolution = 400,
#'   hosek      = FALSE,
#'   altitude   = 0,
#'   visibility = 80,
#'   albedo     = 0.2,
#'   elevation  = 5,
#'   azimuth    = 220,
#'   numbercores = 2
#' ) |>
#'   rayimage::plot_image()
#' }
#' }
generate_sky = function(
  filename = NA,
  albedo = 0.1,
  turbidity = 3,
  elevation = 10,
  azimuth = 90,
  altitude = 0,
  resolution = 2048,
  numbercores = 1,
  hosek = TRUE,
  wide_spectrum = FALSE,
  visibility = 50,
  verbose = FALSE,
  render_mode = "all",
  below_horizon = TRUE
) {
  render_mode = normalize_render_mode(render_mode)
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
      if (verbose) {
        message(
          "Drawing black image as Hosek model does not produce valid output for elevation < 0."
        )
      }
      black_sky = array(0, dim = c(resolution, resolution * 2, 4))
      black_sky[,, 4] = 1
      if (!is.na(filename)) {
        warn_precision_loss(filename)
        rayimage::ray_write_image(black_sky, filename)
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
      if (!is.na(filename)) {
        warn_precision_loss(filename)
        rayimage::ray_write_image(black_sky, filename)
        return(invisible(black_sky))
      } else {
        return(black_sky)
      }
    }
  }

  lambda_values = seq(360, 720, by = 40)

  generated_rgb = makesky_rcpp(
    albedo = albedo,
    turbidity = turbidity,
    elevation = elevation,
    azimuth_deg = azimuth,
    resolution = resolution,
    numbercores = numbercores,
    model = model,
    prg_dataset = coef_file,
    altitude = altitude,
    visibility = visibility,
    render_mode = render_mode,
    lambda_nm = lambda_values,
    below_horizon = below_horizon
  )
  band = attr(generated_rgb, "L_band")

  generated_sky = array(0, dim = c(resolution, resolution * 2, 4))
  generated_sky[,, 1:3] = generated_rgb
  generated_sky[,, 4] = 1
  if (!is.null(band)) {
    attr(generated_sky, "L_band") = band
  }

  if (!is.na(filename)) {
    warn_precision_loss(filename)
    rayimage::ray_write_image(generated_sky, filename)
    return(invisible(generated_sky))
  } else {
    converted = rayimage::ray_read_image(
      generated_sky,
      assume_white = "D65",
      assume_colorspace = rayimage::CS_SRGB
    )
    if (!is.null(band)) {
      attr(converted, "L_band") = band
    }
    return(converted)
  }
}

#' Generate a location and time-specific sky dome (optionally with stars)
#'
#' @description Convenience wrapper around [generate_sky()] that:
#' 1. Computes the Sun's apparent position for `datetime`, `lat`, and `lon`
#'    (via Swiss Ephemeris / **swephR**).
#' 2. Renders the corresponding sky model.
#' 3. Optionally overlays a star field using [generate_stars()] or moon with [generate_moon_latlong()].
#'
#' @param datetime           POSIX-compatible date-time.
#' @param lat                Observer latitude (degrees N).
#' @param lon                Observer longitude (degrees E; west < 0).
#' @param filename            Default `NA`. Path to the `.exr` file to write.
#' @param albedo             Default `0.5`. Ground albedo, range 0 to 1.
#' @param turbidity          Default `3`. Atmospheric turbidity, range
#'   1.7 to 10 (*Hosek only*).
#' @param altitude           Default `0`. Observer altitude (m), range
#'   0 to 15000 (*Prague only*).
#' @param resolution         Default `2048`. Image height in pixels (width = 2 * height).
#' @param numbercores        Default `1`. CPU threads to use.
#' @param hosek              Default `TRUE`. `FALSE` selects the Prague model.
#' @param wide_spectrum      Default `FALSE`. 55-channel Prague coefficients (altitude = 0m only).
#' @param visibility         Default `50`. Meteorological range (km); *Prague only*.
#' @param stars              Default `FALSE`. If `TRUE`, composite a star field
#'   using [generate_stars()]. Automatically falls back to a pure star render
#'   when the black-sky condition is met (see *Details*).
#' @param star_width         Default `1`. Passed to [generate_stars()] to control stellar point-spread size.
#' @param planets            Default `FALSE`. If `TRUE`, composite bright planets via [generate_planets()].
#' @param moon               Default `FALSE`. If `TRUE`, overlay a moon render from [generate_moon_latlong()].
#' @param moon_atmosphere    Default `FALSE`. If `TRUE`, this generates atmospheric scattering from light from the moon.
#' @param moon_hosek         Default `TRUE`. Whether to use the faster (but less accurate) Hosek model for atmospheric scattering from the moon. Note
#' that the light scattered from the moon is much less intense than the sun, and thus small inaccuracies are much less likely to be noticable.
#' @param render_mode        Default `"all"`. One of `"all"`, `"atmosphere"`, or `"sun"`.
#'   Use `"all"` for atmosphere + solar disk, `"atmosphere"` for atmospheric radiance only, or `"sun"` for the solar disk only.
#' @param below_horizon      Default `TRUE`. Whether to sample atmospheric scattering below the horizon, which is non-zero when altitude > 0.
#' @param stars_exposure     Default `0`. Increases star exposure by `2^exposure`. Non-physical, this just controls adjustments for artistic effect.
#' @param verbose            Default `FALSE`. Whether to print progress bars/diagnostic info.
#' @param ...                Additional arguments passed to [generate_stars()], and when enabled, [generate_planets()] and [generate_moon_latlong()].
#'
#' @details
#' *Solar angles* - altitude (degrees above the horizon) and azimuth (degrees clockwise from
#' east, so 90 degrees = south) - are derived internally; you never have to supply
#' them directly.
#'
#' *Black-sky rule* - With the Prague model the sky radiance is defined only
#' down to -4.2 degrees, and with the Hosek model it is defined only about 0 degrees.
#' Below that the function skips the sky render and writes **only stars** when `stars = TRUE`.
#'
#' @return Either the raw data, or the data is invisibly returned if filename is given. The EXR is written to `filename`.
#' @export
#'
#' @examples
#' # Morning sunrise on spring solstice over Washington, DC with Prague model
#' if(run_documentation()) {
#' generate_sky_latlong(
#'   datetime    = as.POSIXct("2025-03-21 06:15:00",tz="EST"),
#'   lat         = 38.9072,
#'   lon         = -77.0369,
#'   numbercores = 2,
#'   hosek = FALSE
#' ) |>
#'   rayimage::plot_image()
#'}
#'if(run_documentation()) {
#' generate_sky_latlong(
#'   datetime    = as.POSIXct("2025-03-21 12:00:00",tz="EST"),
#'   lat         = 38.9072,
#'   lon         = -77.0369,
#'   numbercores = 2,
#' ) |>
#'   rayimage::plot_image()
#'}
#'if(run_documentation()) {
#' generate_sky_latlong(
#'   datetime    = as.POSIXct("2025-03-21 18:00:00",tz="EST"),
#'   lat         = 38.9072,
#'   lon         = -77.0369,
#'   numbercores = 2,
#' ) |>
#'   rayimage::plot_image()
#'}
#'if(run_documentation()) {
#' generate_sky_latlong(
#'   datetime    = as.POSIXct("2025-03-21 18:30:00",tz="EST"),
#'   lat         = 38.9072,
#'   lon         = -77.0369,
#'   numbercores = 2,
#'   hosek = FALSE,
#'   verbose=TRUE,
#' ) |>
#'   rayimage::render_exposure(exposure=2) |>
#'   rayimage::plot_image()
#' }
generate_sky_latlong = function(
  datetime,
  lat,
  lon,
  filename = NA,
  albedo = 0.5,
  turbidity = 3,
  altitude = 0,
  resolution = 2048,
  numbercores = 1,
  hosek = TRUE,
  wide_spectrum = FALSE,
  visibility = 50,
  stars = FALSE,
  star_width = 1.0,
  planets = FALSE,
  moon = FALSE,
  moon_atmosphere = FALSE,
  moon_hosek = TRUE,
  render_mode = "all",
  below_horizon = TRUE,
  verbose = FALSE,
  stars_exposure = 0,
  ...
) {
  render_mode = normalize_render_mode(render_mode)

  sunmoon_data = swe_dirs_topo_moon_sun(
    datetime = datetime,
    lat = lat,
    lon = lon,
    elev_m = altitude
  )
  sun_dir = sunmoon_data$sun_dir_topo
  elevation = asin(sun_dir[3]) * 180 / pi
  azimuth = (180 + atan2(sun_dir[1], sun_dir[2]) * 180 / pi + 360) %% 360

  if (verbose) {
    message(sprintf(
      "Sun: %0.1f elevation, %0.1f azimuth",
      elevation,
      azimuth
    ))
  }
  # Just add up all three
  sky_array = generate_sky(
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
    verbose = verbose,
    render_mode = render_mode,
    below_horizon = below_horizon
  )

  if (moon) {
    moon_array = generate_moon_latlong(
      filename = filename,
      datetime = datetime,
      lat = lat,
      lon = lon,
      albedo = albedo,
      turbidity = turbidity,
      altitude = altitude,
      resolution = resolution,
      numbercores = numbercores,
      hosek = moon_hosek,
      moon_atmosphere = moon_atmosphere,
      wide_spectrum = wide_spectrum,
      visibility = visibility,
      verbose = verbose,
      ...
    )
    sky_band = attr(sky_array, "L_band")
    moon_band = attr(moon_array, "L_band")
    sky_array = sky_array + moon_array
    if (!is.null(sky_band) && !is.null(moon_band)) {
      attr(sky_array, "L_band") = sky_band + moon_band
    }
  }
  if (stars) {
    stars_array = generate_stars(
      resolution = resolution,
      lon = lon,
      lat = lat,
      datetime = as.POSIXct(datetime),
      turbidity = turbidity,
      altitude = altitude,
      numbercores = numbercores,
      star_width = star_width,
      ...
    ) *
      (2^stars_exposure)
    sky_array = sky_array + stars_array
  }
  if (planets) {
    planets_array = generate_planets(
      resolution = resolution,
      datetime = as.POSIXct(datetime),
      lon = lon,
      lat = lat,
      turbidity = turbidity,
      altitude = altitude,
      numbercores = numbercores,
      planet_width = star_width,
      verbose = verbose,
      ...
    )
    sky_array = sky_array + planets_array
  }
  sky_array[,, 4] = 1
  if (!is.na(filename)) {
    warn_precision_loss(filename)
    rayimage::ray_write_image(sky_array, filename, clamp = FALSE)
    return(invisible(sky_array))
  } else {
    return(sky_array)
  }
}

#' Sample a direction from the Prague model.
#'
#' @description Evaluate the Prague spectral sky model at arbitrary spherical
#' directions without writing an image, returning radiance-only samples.
#'
#' @param phi                Horizontal angle of the sample, degrees. Vectorized. Range 0 to 360.
#' @param theta              Vertical angle of the sample, degrees. Vectorized. Range -90 to 90.
#' @param altitude 	         Default `0`, vectorized. Altitude of the viewer in meters. Range 0 to 15000.
#' @param elevation          Default `10`, vectorized. Solar elevation angle above/below the horizon (degrees). Range -4.2 to 90.
#' @param visibility         Default `50`, vectorized. Range 20 to 131.8. Meteorological range in kilometers for Prague model.
#' @param albedo             Default `0.5`, vectorized. Range 0 to 1. Ground albedo.
#' @param azimuth            Default `90`, single value. Solar azimuth (degrees). Defaults South.
#' @param numbercores        Default `1`. Number of threads to use in computation.
#' @param wide_spectrum      Default `FALSE`. Whether to use the wide-spectrum (55-channel, polarised) coefficients.
#' @param render_mode        Default `"all"`. One of `"all"`, `"atmosphere"`, or `"sun"`.
#'   Use `"all"` for atmosphere + solar disk, `"atmosphere"` for atmospheric radiance only, or `"sun"` for the solar disk only.
#'
#' @return 3-column RGB matrix.
#' @export
#' @examples
#' # Generate a basic atmosphere with the Prague model
#' if(run_documentation()) {
#' value_grid = expand.grid(
#'   phi = seq(0,360,by=30),
#'   theta = seq(0,90,by=10),
#'   altitude = c(0,10000)
#'  )
#' vals = calculate_sky_values(
#'   phi = value_grid$phi,
#'   theta = value_grid$theta,
#'   altitude = value_grid$altitude,
#'   elevation = 45,
#'   visibility = 120,
#'   albedo = 0
#'  )
#'  cbind(value_grid, vals)
#' }
calculate_sky_values = function(
  phi,
  theta,
  altitude = 0,
  elevation = 10,
  visibility = 50,
  albedo = 0.5,
  azimuth = 90,
  numbercores = 1,
  wide_spectrum = FALSE,
  render_mode = "all"
) {
  render_mode = normalize_render_mode(render_mode)
  sea_level = all(altitude == 0)
  stopifnot(all(phi <= 360 & phi >= 0))
  stopifnot(all(theta <= 90 & theta >= -90))
  stopifnot(all(altitude <= 15000 & altitude >= 0))
  stopifnot(all(visibility >= 20 & visibility <= 131.8))
  filesize = ""
  df_values = as.list(data.frame(
    phi = phi,
    theta = theta,
    altitude = altitude,
    elevation = elevation,
    visibility = visibility,
    albedo = albedo,
    azimuth = azimuth
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
  vals = calculate_raw_prague(
    df_values$phi,
    df_values$theta,
    df_values$elevation,
    df_values$albedo,
    df_values$altitude,
    df_values$visibility,
    df_values$azimuth,
    numbercores,
    coef_file,
    render_mode
  )
  colnames(vals) = c("r", "g", "b")
  return(vals)
}

#' Sample sun luminance at the disk center.
#'
#' @description Evaluate the Prague or Hosek sky model at the sun center and
#' integrate against the CIE Y curve to return a luminance value. This is
#' useful for relative attenuation (e.g., comparing zenith sun to a low sun).
#'
#' @param elevation          Default `10`. Solar elevation angle above/below
#'   the horizon (degrees). Range `[-4.2, 90]` for Prague, `[0, 90]` for Hosek.
#' @param azimuth            Default `90`. Solar azimuth (degrees).
#' @param albedo             Default `0.5`. Ground albedo, range 0 to 1.
#' @param turbidity          Default `3`. Atmospheric turbidity, range 1.7 to 10
#'   (*Hosek only*).
#' @param altitude           Default `0`. Observer altitude (m), range 0 to 15000
#'   (*Prague only*).
#' @param visibility         Default `50`. Meteorological range (km); *Prague only*.
#' @param hosek              Default `TRUE`. `FALSE` selects the Prague model.
#' @param wide_spectrum      Default `FALSE`. Use wide-spectrum (55-channel)
#'   coefficients for Prague at sea level only.
#' @param lambda_nm          Optional vector of wavelengths for Hosek sampling.
#'
#' @return Numeric scalar of sun luminance (CIE Y, relative scale).
#' @export
#' @examples
#' if(run_documentation()) {
#'   calculate_sun_brightness(elevation = 45, hosek = TRUE)
#' }
calculate_sun_brightness = function(
  elevation = 10,
  azimuth = 90,
  albedo = 0.5,
  turbidity = 3,
  altitude = 0,
  visibility = 50,
  hosek = TRUE,
  wide_spectrum = FALSE,
  lambda_nm = NULL
) {
  stopifnot(length(elevation) == 1)
  stopifnot(length(azimuth) == 1)
  stopifnot(length(albedo) == 1)
  stopifnot(length(turbidity) == 1)
  stopifnot(length(altitude) == 1)
  stopifnot(length(visibility) == 1)
  stopifnot(altitude >= 0)

  model = if (hosek) "hosek" else "prague"
  coef_file = ""

  if (!hosek) {
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
  }

  calculate_sun_brightness_rcpp(
    albedo = albedo,
    turbidity = turbidity,
    elevation = elevation,
    azimuth_deg = azimuth,
    model = model,
    prg_dataset = coef_file,
    altitude = altitude,
    visibility = visibility,
    lambda_nm = if (hosek) lambda_nm else NULL
  )
}
