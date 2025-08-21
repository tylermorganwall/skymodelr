#' Generate the atmosphere with the moon
#'
#' Note that this is just a scaled version of `generate_sky()`, scaled down by the luminance
#' of the moon as compared to the sun. This function takes the phase of the moon into account,
#' along with the increase in luminosity around a full moon (known as opposition surge).
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
#' @param verbose            Default `FALSE`. Whether to print progress bars/diagnostic info.
#' @param ...                Additional arguments passed to [generate_stars()]
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
generate_moon_latlong = function(
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
  verbose = FALSE,
  ...
) {
  moon_altitude_azimuth = suncalc::getMoonPosition(datetime, lat, lon)
  moon_elevation = moon_altitude_azimuth$altitude * 180 / pi
  moon_azimuth = 90 + moon_altitude_azimuth$azimuth * 180 / pi
  moon_phase_info = suncalc::getMoonIllumination(datetime) # 0.5 = full
  moon_phase = ((moon_phase_info$phase - 0.5) * 2 * pi) * 180 / pi
  #sun is 10,000 footcandles/100,000 lux at noon
  sun_brightness_ftcndl = 10000
  # m = -12.7 at angle == 0, which is correct
  m = -12.73 + 0.026 * abs(moon_phase) + 4 * 10^(-9) * moon_phase^4
  #opposition_surge_multiplier is 1.35 at max, 1.0 for continuity
  opposition_surge_multiplier = max(1.35 - 0.05 * moon_phase, 1)
  moon_ftcndl = opposition_surge_multiplier * 10^(-0.4 * (m + 16.57))

  scale_sun_to_moon = moon_ftcndl / sun_brightness_ftcndl
  if (verbose) {
    message(sprintf(
      "Moon: %0.1f° elevation, %0.1f° azimuth, %0.1f phase, %f lux",
      moon_elevation,
      moon_azimuth,
      moon_phase,
      moon_brightness_ftcndl * 10
    ))
  }
  moon_exr = generate_sky(
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
    square_projection = square_projection,
    verbose = verbose,
    render_solar_disk = TRUE
  )
  moon_exr[,, 1:3] = moon_exr[,, 1:3] * scale_sun_to_moon

  if (!is.na(outfile)) {
    rayimage::ray_write_image(moon_exr, outfile, clamp = FALSE)
    return(invisible(moon_exr))
  } else {
    return(moon_exr)
  }
}
