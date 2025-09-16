#' @keywords internal
normalize = function(v) {
  s = sqrt(sum(v * v))
  if (s == 0) stop("Zero-length vector")
  return(v / s)
}

#' @keywords internal
cross_prod = function(x, y) {
  return(c(
    x[2] * y[3] - x[3] * y[2],
    -(x[1] * y[3] - x[3] * y[1]),
    x[1] * y[2] - x[2] * y[1]
  ))
}

#Manual implementation:
altaz_to_enu = function(alt_rad, azN_rad) {
  east = sin(azN_rad) * cos(alt_rad)
  north = cos(azN_rad) * cos(alt_rad)
  up = sin(alt_rad)
  normalize(c(east, north, up)) # (E, N, U) matches engine XYZ
}

swe_dirs_topo_moon_sun = function(datetime, lat, lon, elev_m = 0) {
  swephR::swe_set_topo(lon, lat, elev_m)
  attr(datetime, "tzone") = "UTC"
  ts = as.POSIXlt(datetime, tz = "UTC")
  jd_ut = swephR::swe_utc_to_jd(
    ts$year + 1900,
    ts$mon + 1,
    ts$mday,
    ts$hour,
    ts$min,
    ts$sec,
    swephR::SE$GREG_CAL
  )$dret[1]

  flg_eq_topo = swephR::SE$FLG_SWIEPH +
    swephR::SE$FLG_EQUATORIAL +
    swephR::SE$FLG_TOPOCTR
  flg_eq_geo = swephR::SE$FLG_SWIEPH + swephR::SE$FLG_EQUATORIAL
  long_rad = lon * pi / 180
  lat_rad = lat * pi / 180

  enu_to_ecef = matrix(
    c(
      -sin(long_rad),
      cos(long_rad),
      0,
      -sin(lat_rad) * cos(long_rad),
      -sin(lat_rad) * sin(long_rad),
      cos(lat_rad),
      cos(lat_rad) * cos(long_rad),
      cos(lat_rad) * sin(long_rad),
      sin(lat_rad)
    ),
    ncol = 3,
    nrow = 3,
    byrow = TRUE
  )
  sun_t = swephR::swe_calc_ut(jd_ut, swephR::SE$SUN, flg_eq_topo)$xx
  moon_t = swephR::swe_calc_ut(jd_ut, swephR::SE$MOON, flg_eq_topo)$xx
  moon_size = swephR::swe_pheno_ut(jd_ut, swephR::SE$MOON, flg_eq_topo)$attr
  sun_size = swephR::swe_pheno_ut(jd_ut, swephR::SE$SUN, flg_eq_topo)$attr

  moon_diameter_degrees = moon_size[4]
  moon_brightness_magnitude = moon_size[5]

  sun_diameter_degrees = sun_size[4]
  sun_brightness_magnitude = sun_size[5]

  to_local_horizontal_coordinates = function(ra_dec) {
    aa = swephR::swe_azalt(
      jd_ut,
      swephR::SE$EQU2HOR,
      c(lon, lat, elev_m),
      0,
      0,
      xin = c(ra_dec[1], ra_dec[2], 1)
    )$xaz
    # aa[1] = (aa[1] + 180) %% 360 # south-based → north-based
    aa
  }
  to_rad = pi / 180

  local_horizon_sun = to_local_horizontal_coordinates(sun_t) * to_rad
  local_horizon_moon = to_local_horizontal_coordinates(moon_t) * to_rad

  azN_s = local_horizon_sun[1]
  alt_s = local_horizon_sun[2]
  azN_mt = local_horizon_moon[1]
  alt_mt = local_horizon_moon[2]

  list(
    sun_dir_topo = altaz_to_enu(alt_s, azN_s),
    sun_ecef_geo = enu_to_ecef %*% altaz_to_enu(alt_s, azN_s),
    moon_dir_topo = altaz_to_enu(alt_mt, azN_mt),
    moon_ecef_geo = enu_to_ecef %*% altaz_to_enu(alt_mt, azN_mt),
    local_up_geo = enu_to_ecef %*% c(0, 0, 1),
    moon_diameter_degrees = moon_diameter_degrees,
    moon_brightness_magnitude = moon_brightness_magnitude,
    sun_diameter_degrees = sun_diameter_degrees,
    sun_brightness_magnitude = sun_brightness_magnitude
  ) |>
    lapply(as.numeric)
}


sweph_info_to_stars_df = function(sweph_info) {
  sweph_info_df = do.call(rbind, lapply(sweph_info, as.data.frame))
  planets = rownames(sweph_info_df)
  data.frame(
    planet = planets,
    ra_rad = sweph_info_df$ra_rad,
    dec_rad = sweph_info_df$dec_rad,
    v_mag = sweph_info_df$v_mag,
    r = 1,
    g = 1,
    b = 1
  )
}

swe_dirs_topo_planets_df = function(datetime, lat, lon, elev_m = 0) {
  swephR::swe_set_topo(lon, lat, elev_m)
  attr(datetime, "tzone") = "UTC"
  ts = as.POSIXlt(datetime, tz = "UTC")
  jd_ut = swephR::swe_utc_to_jd(
    ts$year + 1900,
    ts$mon + 1,
    ts$mday,
    ts$hour,
    ts$min,
    ts$sec,
    swephR::SE$GREG_CAL
  )$dret[1]

  flg_eq_topo = swephR::SE$FLG_SWIEPH +
    swephR::SE$FLG_EQUATORIAL +
    swephR::SE$FLG_TOPOCTR
  flg_eq_geo = swephR::SE$FLG_SWIEPH + swephR::SE$FLG_EQUATORIAL
  long_rad = lon * pi / 180
  lat_rad = lat * pi / 180

  enu_to_ecef = matrix(
    c(
      -sin(long_rad),
      cos(long_rad),
      0,
      -sin(lat_rad) * cos(long_rad),
      -sin(lat_rad) * sin(long_rad),
      cos(lat_rad),
      cos(lat_rad) * cos(long_rad),
      cos(lat_rad) * sin(long_rad),
      sin(lat_rad)
    ),
    ncol = 3,
    nrow = 3,
    byrow = TRUE
  )

  to_rad = pi / 180

  to_local_horizontal_coordinates = function(ra_dec) {
    aa = swephR::swe_azalt(
      jd_ut,
      swephR::SE$EQU2HOR,
      c(lon, lat, elev_m),
      0,
      0,
      xin = c(ra_dec[1], ra_dec[2], 1)
    )$xaz
    aa
  }
  get_enu_magnitude = function(planet_flag) {
    position_rad = #to_local_horizontal_coordinates(
      swephR::swe_calc_ut(jd_ut, planet_flag, flg_eq_topo)$xx * to_rad
    # ) *
    # to_rad
    # position_enu = altaz_to_enu(position_hori[1], position_hori[2])
    info = swephR::swe_pheno_ut(jd_ut, planet_flag, flg_eq_topo)$attr
    list(ra_rad = position_rad[1], dec_rad = position_rad[2], v_mag = info[5])
  }
  planet_info = list(
    mercury = get_enu_magnitude(swephR::SE$MERCURY),
    venus = get_enu_magnitude(swephR::SE$VENUS),
    mars = get_enu_magnitude(swephR::SE$MARS),
    jupiter = get_enu_magnitude(swephR::SE$JUPITER),
    saturn = get_enu_magnitude(swephR::SE$SATURN),
    uranus = get_enu_magnitude(swephR::SE$URANUS),
    neptune = get_enu_magnitude(swephR::SE$NEPTUNE)
  )
  return(sweph_info_to_stars_df(planet_info))
}


build_from_z = function(x) {
  zz = normalize(x)
  a = if (abs(zz[1]) > 0.9999999) c(0, 1, 0) else c(1, 0, 0)
  yy = normalize(cross_prod(zz, a))
  xx = cross_prod(yy, zz)
  return(matrix(c(xx, yy, zz), ncol = 3, byrow = FALSE))
}

generate_moon_image_latlong = function(
  datetime,
  lat,
  lon,
  elev_m = 0,
  width = 400,
  height = 400
) {
  moon_sun_data = swe_dirs_topo_moon_sun(datetime, lat, lon, elev_m)
  dir_moon = normalize(moon_sun_data$moon_ecef_geo)
  # dir_sun = normalize(
  #   moon_sun_data$sun_ecef_geo - moon_sun_data$moon_ecef_geo
  # )
  dir_sun = -normalize(moon_sun_data$sun_ecef_geo)
  local_up = normalize(moon_sun_data$local_up_geo)

  moon_material = rayvertex::material_list(
    texture_location = "lroc_color_poles_1k.jpg",
    sigma = 10
  )

  rayvertex::sphere_mesh(material = moon_material, radius = 0.5) |>
    rayvertex::subdivide_mesh() |>
    rayvertex::rotate_mesh(c(0, -90, 6.68), order_rotation = c(3, 1, 2)) |>
    rayvertex::transform_mesh(rayvertex::lookat_transform(
      pos = -dir_moon,
      look = c(0, 0, 0),
      up = c(0, 0, 1)
    )) |>
    rayvertex::translate_mesh(dir_moon * 10) |>
    rayvertex::rasterize_scene(
      lookfrom = c(0, 0, 0),
      lookat = dir_moon * 10,
      camera_up = local_up,
      fov = 0,
      light_info = rayvertex::directional_light(
        direction = -dir_sun,
        intensity = 1
      ),
      transparent_background = TRUE,
      plot = FALSE,
      width = width,
      height = height
    ) -> moon_image
  mag_to_lux = function(m) {
    10^((-14.18 - m) / 2.5)
  }
  moon_mean = mean(rayimage::render_bw(moon_image)[,, 1:3])
  moon_image[,, 1:3] = (moon_image[,, 1:3] / moon_mean) *
    mag_to_lux(moon_sun_data$moon_brightness_magnitude)
  list(
    moon_luminance_array = moon_image,
    moon_angular_diameter_deg = moon_sun_data$moon_diameter_degrees
  )
}
