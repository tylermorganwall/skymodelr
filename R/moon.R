#' Generate the atmosphere with the moon
#'
#' Note that this is just a scaled version of `generate_sky()`, scaled down by the luminance
#' of the moon as compared to the sun. This function takes the phase of the moon into account,
#' along with the increase in luminosity around a full moon (known as opposition surge).
#'
#' @param filename           Default `NA`. Path to the `.exr` file to write.
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
#' @param verbose            Default `FALSE`. Whether to print progress bars/diagnostic info.
#' @param ...                Additional arguments passed to [generate_stars()]
#'
#' @return Either the raw data, or the data is invisibly returned if filename is given. The EXR is written to `filename`.
#' @export
#'
#' @examples
#' # Moonlit sky (Hosek), mid-evening in DC
#' if(run_documentation()) {
#' generate_moon_latlong(
#'   datetime   = as.POSIXct("2025-09-05 19:30:00",tz="America/New_York"),
#'   lat        = 38.9072,
#'   lon        = -77.0369,
#'   resolution = 400,
#'   turbidity  = 3,
#'   verbose    = TRUE
#' ) |>
#'   rayimage::render_exposure(15) |>
#'   rayimage::plot_image()
#' }
generate_moon_latlong = function(
  filename = NA,
  datetime = "2025-07-29 18:00:00",
  lat = 38.9072,
  lon = -77.0369,
  albedo = 0.5,
  turbidity = 3,
  altitude = 0,
  resolution = 2048,
  numbercores = 1,
  atmospheric_scattering = FALSE,
  hosek = TRUE,
  wide_spectrum = FALSE,
  visibility = 50,
  verbose = FALSE,
  ...
) {
  moon_altitude_azimuth = suncalc::getMoonPosition(datetime, lat, lon)
  moon_elevation = moon_altitude_azimuth$altitude * 180 / pi
  moon_azimuth = 90 + moon_altitude_azimuth$azimuth * 180 / pi
  if (moon_azimuth < 0) {
    moon_azimuth = moon_azimuth + 360
  }
  moon_sun_data = swe_dirs_topo_moon_sun(datetime, lat, lon, altitude)

  mag_to_lux = function(m) {
    10^((-14.18 - m) / 2.5)
  }
  moon_lux = mag_to_lux(moon_sun_data$moon_brightness_magnitude)
  sun_lux = mag_to_lux(moon_sun_data$sun_brightness_magnitude)

  moon_info_list = generate_moon_image_latlong(
    datetime,
    lat,
    lon,
    altitude,
    width = 401,
    height = 401
  )
  moon_luminance_array = moon_info_list$moon_luminance_array
  moon_angular_diameter_deg = moon_info_list$moon_angular_diameter_deg
  moon_angular_diameter_rad = moon_angular_diameter_deg * pi / 180

  scale_sun_to_moon = moon_lux / sun_lux
  if (verbose) {
    message(sprintf(
      "Moon: %0.1f° elevation, %0.1f° azimuth, %0.1f phase, %f lux",
      moon_elevation,
      moon_azimuth,
      0,
      # moon_phase,
      moon_lux
    ))
  }
  if (atmospheric_scattering) {
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
      verbose = verbose,
      render_solar_disk = FALSE
    )
  } else {
    moon_exr = rayimage::ray_read_image(array(
      0,
      dim = c(resolution, resolution * 2, 4)
    ))
    moon_exr[,, 4] = 1
  }
  moon_exr[,, 1:3] = moon_exr[,, 1:3] * scale_sun_to_moon

  moon_dir_vec = c(
    cospi(moon_azimuth / 180) * cospi(moon_elevation / 180),
    sinpi(moon_elevation / 180),
    sinpi(moon_azimuth / 180) * cospi(moon_elevation / 180)
  )

  moon_diameter_pixels = resolution / 180 * moon_angular_diameter_deg
  resize_moon_dim = ceiling(moon_diameter_pixels)
  if (resize_moon_dim %% 2 == 0) {
    resize_moon_dim = resize_moon_dim + 1
  }

  resized_moon_luminance_array = rayimage::render_resized(
    moon_luminance_array,
    dims = c(resize_moon_dim, resize_moon_dim)
  )

  nTheta = resolution
  nPhi = 2 * resolution

  r = moon_angular_diameter_rad / 2
  cos_r = cos(r)

  phi_c = moon_elevation * pi / 180
  j_min = max(1, floor(((pi / 2) - (phi_c + r)) / pi * (nTheta - 1)) + 1)
  j_max = min(nTheta, ceiling(((pi / 2) - (phi_c - r)) / pi * (nTheta - 1)) + 1)

  thetas = seq(pi / 2, -pi / 2, length.out = nTheta)
  phis = seq(0, 2 * pi, length.out = nPhi)

  cos_phi = cos(phis)
  sin_phi = sin(phis)
  cos_theta = cos(thetas)
  sin_theta = sin(thetas)

  # Tangent frame at moon center (ẑ = moon_dir_vec)
  lam_c = moon_azimuth * pi / 180
  z_hat = moon_dir_vec
  e_hat = c(-sin(lam_c), 0, cos(lam_c))
  e_hat = e_hat / sqrt(sum(e_hat * e_hat)) # local east
  n_hat = c(
    z_hat[2] * e_hat[3] - z_hat[3] * e_hat[2],
    z_hat[3] * e_hat[1] - z_hat[1] * e_hat[3],
    z_hat[1] * e_hat[2] - z_hat[2] * e_hat[1]
  )
  n_hat = n_hat / sqrt(sum(n_hat * n_hat)) # local north

  tan_r = tan(r) # r = moon_angular_diameter_rad/2

  for (j in seq(j_min, j_max)) {
    for (i in seq_len(nPhi)) {
      sample_dir_vec = c(
        cos_phi[i] * cos_theta[j],
        sin_theta[j],
        sin_phi[i] * cos_theta[j]
      )
      dot = sum(moon_dir_vec * sample_dir_vec)
      if (dot < cos_r) next

      # Camera (tangent) components of the sample direction
      vx = sum(sample_dir_vec * e_hat)
      vy = sum(sample_dir_vec * n_hat)
      vz = sum(sample_dir_vec * z_hat) # == dot

      # Inverse rectilinear mapping (square FOV 2*r):
      u_px = 201 + 200 * (vx / (tan_r * vz))
      v_px = 201 + 200 * (vy / (tan_r * vz))

      u_i = rayrender:::clamp(round(u_px), 1, 401)
      v_i = rayrender:::clamp(round(v_px), 1, 401)

      moon_exr[j, i, 1] = moon_luminance_array[v_i, u_i, 1]
      moon_exr[j, i, 2] = moon_luminance_array[v_i, u_i, 2]
      moon_exr[j, i, 3] = moon_luminance_array[v_i, u_i, 3]
    }
  }
  if (!is.na(filename)) {
    rayimage::ray_write_image(moon_exr, filename, clamp = FALSE)
    return(invisible(moon_exr))
  } else {
    return(moon_exr)
  }
}
