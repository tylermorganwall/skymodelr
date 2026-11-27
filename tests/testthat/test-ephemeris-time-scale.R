test_that("topocentric positions use the correct Swiss Ephemeris time scales", {
  datetime = as.POSIXct("2024-04-08 15:29:00", tz = "America/New_York")
  lat = 44.906907
  lon = -72.804280
  coords = c(lon, lat, 0)
  se = swephR::SE
  flags = se$FLG_SWIEPH + se$FLG_EQUATORIAL + se$FLG_TOPOCTR

  # Reference: the fixture is 19:29 UTC. Use swe_calc() with TT for
  # equatorial positions and UT1 only for the horizontal conversion.
  jd = swephR::swe_utc_to_jd(2024, 4, 8, 19, 29, 0, se$GREG_CAL)$dret
  swephR::swe_set_topo(lon, lat, 0)
  sunmoon = swe_dirs_topo_sunmoon_df(datetime, lat, lon)
  planets = swe_dirs_topo_planets_df(datetime, lat, lon)
  directions = swe_dirs_topo_moon_sun(datetime, lat, lon)
  bodies = c(sun = se$SUN, moon = se$MOON, venus = se$VENUS,
             mars = se$MARS, jupiter = se$JUPITER)

  for (body in names(bodies)) {
    eq = swephR::swe_calc(jd[1], bodies[[body]], flags)$xx
    positions = if (body %in% c("sun", "moon")) sunmoon else planets
    row = positions[positions$planet == body, , drop = FALSE]
    expect_equal(nrow(row), 1L, info = body)
    expect_lt(max(abs(c(row$ra_deg, row$dec_deg) - eq[1:2])), 1e-6)

    if (body %in% c("sun", "moon")) {
      horiz = swephR::swe_azalt(
        jd[2], se$EQU2HOR, coords, 0, 0, xin = c(eq[1:2], 1)
      )$xaz
      dir = directions[[paste0(body, "_dir_topo")]]
      az = atan2(dir[1], dir[2]) * 180 / pi
      el = asin(dir[3]) * 180 / pi
      az_error = abs((az - horiz[1] + 180) %% 360 - 180)
      expect_lt(az_error, 1e-6)
      expect_lt(abs(el - horiz[2]), 1e-6)
    }
  }
})
