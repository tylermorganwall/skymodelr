time <- as.POSIXct('2026-06-21 19:00:00', tz = 'America/New_York')
test_that('metadata agrees with the existing public Sun and calibration APIs', {
  files <- vapply(
    c(0, 1),
    function(altitude) {
      tryCatch(
        resolve_prague_coef_file(altitude = altitude, allow_download = FALSE),
        error = function(e) NA_character_
      )
    },
    character(1)
  )
  skip_if(
    anyNA(files) || !all(file.exists(files)),
    'Prague coefficient datasets are not installed'
  )
  info <- get_prague_sky_metadata(time, 40.7, -74)
  disk <- generate_sun_disk(time, 40.7, -74, resolution = 16)
  for (key in c('elevation_deg', 'azimuth_deg', 'angular_diameter_deg')) {
    expect_equal(info[[key]], disk[[key]])
  }
  calibration <- calculate_sky_values(
    0,
    60,
    elevation = info$elevation_deg,
    azimuth = info$azimuth_deg
  )
  expect_equal(info$rgb_gain, attr(calibration, 'prague_rgb_correction_gain'))
  expect_true(file.exists(info$filename))
  expect_match(info$filename, 'SkyModelDataset.dat', fixed = TRUE)
  expect_equal(
    get_prague_sky_metadata(
      time,
      40.7,
      -74,
      prague_rgb_correction = FALSE
    )$rgb_gain,
    c(R = 1, G = 1, B = 1)
  )
  expect_equal(
    get_prague_sky_metadata(
      time,
      40.7,
      -74,
      prague_rgb_correction_gain = c(1, 2, 3),
      prague_rgb_correction_strength = .5
    )$rgb_gain,
    c(R = 1, G = sqrt(2), B = sqrt(3))
  )
  expect_error(get_prague_sky_metadata(time, 100, -74), 'lat')
  expect_error(
    get_prague_sky_metadata(time, 40.7, -74, visibility = 10),
    'visibility'
  )
  expect_error(
    get_prague_sky_metadata(time, 40.7, -74, altitude = 15001),
    'altitude'
  )
  expect_error(get_prague_sky_metadata('time', 40.7, -74), 'POSIXct')
  expect_error(
    get_prague_sky_metadata(
      as.POSIXct('2026-06-21 00:00:00', tz = 'America/New_York'),
      40.7,
      -74
    ),
    'Sun elevations'
  )
  expect_identical(unserialize(serialize(info, NULL)), info)
})
