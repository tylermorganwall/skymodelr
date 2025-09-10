#Manual implementation:
altaz_to_enu = function(alt_rad, azN_rad) {
  east = sin(azN_rad) * cos(alt_rad)
  north = cos(azN_rad) * cos(alt_rad)
  up = sin(alt_rad)
  rayvertex:::normalize(c(east, north, up)) # (E, N, U) matches engine XYZ
}

swe_dirs_topo_moon_sun = function(datetime, lat, lon, elev_m = 0) {
  swephR::swe_set_topo(lon, lat, elev_m)
  attr(datetime, "tzone") = "UTC"
  ts = as.POSIXlt(datetime, tz = "UTC")
  print(ts)
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
  swephR::SE
  sun_t = swephR::swe_calc_ut(jd_ut, swephR::SE$SUN, flg_eq_topo)$xx
  moon_t = swephR::swe_calc_ut(jd_ut, swephR::SE$MOON, flg_eq_topo)$xx

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
    local_up_geo = enu_to_ecef %*% c(0, 0, 1)
  ) |>
    lapply(as.numeric)
}

build_from_z = function(x) {
  zz = rayvertex:::normalize(x)
  a = if (abs(zz[1]) > 0.9999999) c(0, 1, 0) else c(1, 0, 0)
  yy = rayvertex:::normalize(rayvertex:::cross_prod(zz, a))
  xx = rayvertex:::cross_prod(yy, zz)
  return(matrix(c(xx, yy, zz), ncol = 3, byrow = FALSE))
}

generate_moon_image_latlong = function(
  datetime,
  lat,
  lon,
  elev_m = 0
) {
  enu_vectors = swe_dirs_topo_moon_sun(datetime, lat, lon, elev_m)
  dir_moon = rayvertex:::normalize(enu_vectors$moon_ecef_geo)
  dir_sun = rayvertex:::normalize(
    enu_vectors$sun_ecef_geo - enu_vectors$moon_ecef_geo
  )
  local_up = rayvertex:::normalize(enu_vectors$local_up_geo)

  moon_material = rayvertex::material_list(
    texture_location = "lroc_color_poles_1k.jpg",
    sigma = 10
  )

  rayvertex::sphere_mesh(material = moon_material, radius = 0.5) |>
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
        direction = dir_sun,
        intensity = 1
      )
    )
}
