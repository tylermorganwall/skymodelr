test_that('legacy Moon clips after full-disk calibration for both atmosphere models', {
  elevation = 12
  tint_elevations = numeric()
  old_cores = getOption('cores')
  local_mocked_bindings(
    swe_dirs_topo_moon_sun = function(..., moon_extinction_kV = 0.172) {
      list(
        moon_dir_topo = c(0, -cospi(elevation / 180), sinpi(elevation / 180)),
        sun_dir_topo = c(0, 0, 1),
        moon_diameter_degrees = 20,
        moon_phase = 0,
        moon_brightness_lux = apply_airmass_extinction(
          0.3,
          elevation,
          moon_extinction_kV
        )
      )
    },
    generate_moon_image_latlong = function(...) {
      expect_equal(getOption('cores'), 2)
      list(
        moon_luminance_array = array(1, c(32, 32, 4)),
        moon_angular_diameter_deg = 20
      )
    },
    generate_sky = function(..., elevation, resolution, hosek, render_mode) {
      if (hosek && elevation < 0) {
        stop('invalid Hosek atmospheric elevation')
      }
      image = array(0, c(resolution, 2 * resolution, 4))
      if (render_mode == 'sun') {
        tint_elevations <<- c(tint_elevations, elevation)
        image[,, 1:3] = 1
      }
      image[,, 4] = 1
      attr(image, 'L_band') = image[,, 1] + image[,, 2] + image[,, 3]
      image
    },
    resolve_prague_coef_file = function(...) 'unused.dat',
    calculate_sun_radiance_band_rcpp = function(...) 1,
    lux_to_radiometric_irradiance = function(lux, efficacy) lux / 100
  )
  generate = function(degrees, hosek, atmosphere = FALSE, extinction = 0) {
    elevation <<- degrees
    generate_moon_latlong(
      as.POSIXct('2026-01-04', tz = 'UTC'),
      0,
      0,
      resolution = 180,
      number_cores = 2,
      hosek = hosek,
      moon_atmosphere = atmosphere,
      moon_extinction_kV = extinction,
      moon_texture_width = 32,
      moon_texture_height = 32
    )
  }
  power = function(image) {
    n = dim(image)[1]
    sum(
      rowSums(attr(image, 'L_band')) * cos(seq(pi / 2, -pi / 2, length.out = n))
    ) *
      (pi / n) *
      (2 * pi / (2 * n))
  }
  for (hosek in c(TRUE, FALSE)) {
    full = generate(12, hosek)
    half = generate(0, hosek)
    segment = generate(-5, hosek)
    set = generate(-11, hosek)
    expect_equal(power(full), 0.003, tolerance = 1e-10)
    expect_equal(power(half) / power(full), 0.5, tolerance = 0.03)
    expect_gt(power(segment), 0)
    expect_lt(power(segment), 0.3 * power(full))
    expect_equal(power(set), 0)
    for (image in list(full, half, segment, set)) {
      expect_true(all(is.finite(image)))
      expect_true(all(image[,, 4] == 1))
      expect_equal(max(abs(image[91:180, , 1:3])), 0)
      expect_equal(
        as.numeric(attr(image, 'L_band')),
        as.numeric(image[,, 1] + image[,, 2] + image[,, 3])
      )
    }
    atmospheric = generate(-1, hosek, atmosphere = TRUE)
    expect_gt(power(atmospheric), 0) # Model domain cannot switch off the disk.
    dark = generate(12, hosek, extinction = 1e6)
    expect_equal(power(dark), 0) # Zero target power cannot leave an unscaled texture.
    # The optional atmospheric scaling must not affect the independently calibrated disk.
    expect_equal(generate(12, hosek, atmosphere = TRUE)[,, 1:3], full[,, 1:3])
    elevation = -5
    overlay = generate_sky_latlong(
      as.POSIXct('2026-01-04', tz = 'UTC'),
      0,
      0,
      resolution = 180,
      number_cores = 2,
      moon = TRUE,
      moon_hosek = hosek,
      moon_extinction_kV = 0,
      moon_texture_width = 32,
      moon_texture_height = 32,
      render_mode = 'atmosphere',
      exr_metadata = FALSE
    )
    expect_gt(sum(overlay[,, 1:3]), 0)
    expect_true(all(overlay[,, 4] == 1))
  }
  expect_true(all(tint_elevations >= 0))
  expect_identical(getOption('cores'), old_cores)
})
