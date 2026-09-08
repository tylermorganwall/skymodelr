test_that("disk inputs are validated before model or raster work", {
  time = as.POSIXct("2026-01-28", tz = "UTC")
  defaults = list(datetime = time, lat = 0, lon = 0)
  for (f in list(generate_sun_disk, generate_moon_disk)) {
    for (field in c(
      "lat",
      "lon",
      "altitude",
      "resolution",
      "visibility",
      "albedo",
      "number_cores",
      "prague_rgb_correction_strength"
    )) {
      for (value in list(NA_real_, Inf, numeric(), c(1, 2), "1")) {
        expect_error(
          do.call(f, utils::modifyList(defaults, setNames(list(value), field))),
          field
        )
      }
    }
    for (update in list(
      list(lat = 91),
      list(lon = -181),
      list(altitude = 15001),
      list(resolution = 15),
      list(resolution = 16.5),
      list(visibility = 19),
      list(albedo = 1.1),
      list(number_cores = 1.5),
      list(wide_spectrum = NA),
      list(datetime = as.Date(time)),
      list(prague_rgb_correction = NA),
      list(prague_rgb_correction_gain = c(1, -1, 1))
    )) {
      expect_error(do.call(f, utils::modifyList(defaults, update)))
    }
  }
  for (update in list(
    list(earthshine = NA),
    list(earthshine_albedo = -1),
    list(solar_irradiance_w_m2 = 0),
    list(moon_extinction_kV = Inf)
  )) {
    expect_error(do.call(
      generate_moon_disk,
      utils::modifyList(defaults, update)
    ))
  }
})

test_that("Sun export provides ephemeris geometry and directional solar radiance", {
  seen = NULL
  ephem = list(
    sun_dir_topo = c(0, -cospi(1 / 6), sinpi(1 / 6)),
    sun_diameter_degrees = 0.53
  )
  local_mocked_bindings(
    swe_dirs_topo_moon_sun = function(...) ephem,
    calculate_sky_values = function(...) {
      seen <<- list(...)
      cbind(
        rep(3, length(seen$phi)),
        rep(2, length(seen$phi)),
        rep(-0.1, length(seen$phi))
      )
    }
  )
  disk = generate_sun_disk(
    as.POSIXct("2026-01-28", tz = "UTC"),
    0,
    0,
    resolution = 16
  )
  expect_identical(dim(disk$image), c(16L, 16L, 3L))
  expect_equal(disk$image[8, 8, ], c(3, 2, -0.1))
  expect_equal(disk$azimuth_deg, 0)
  expect_equal(disk$elevation_deg, 30)
  expect_equal(disk$angular_diameter_deg, 0.53)
  expect_identical(disk$projection, "rectilinear")
  expect_equal(seen$render_mode, "sun")
  expect_equal(seen$elevation, 30)
  expect_true(diff(range(seen$theta)) > 0.4)
  expect_true(all(seen$theta > 29.7 & seen$theta < 30.3))
  expect_gt(seen$theta[1], seen$theta[16])
  expect_lt(((seen$phi[1] + 180) %% 360) - 180, 0)
  expect_gt(((seen$phi[241] + 180) %% 360) - 180, 0)
  ephem$sun_dir_topo = c(0, -cospi(1 / 6), -sinpi(1 / 6))
  expect_equal(
    sum(
      generate_sun_disk(
        as.POSIXct("2026-01-28", tz = "UTC"),
        0,
        0,
        resolution = 16
      )$image
    ),
    0
  )
})

test_that("Moon export preserves phase, radiometry, coverage, and thread options", {
  seen = NULL
  old_cores = getOption("cores")
  diameter = 0.53
  local_mocked_bindings(
    swe_dirs_topo_moon_sun = function(...) {
      list(
        moon_dir_topo = c(0, 0, 1),
        moon_diameter_degrees = diameter,
        moon_brightness_lux_unattenuated = 0.3
      )
    },
    generate_moon_image_latlong = function(...) {
      seen <<- c(list(...), list(cores = getOption("cores")))
      patch = array(0, c(32, 32, 4))
      patch[9:24, 9:24, 4] = 1
      for (channel in 1:3) {
        patch[9:24, 9:24, channel] = rep(
          seq(0.1, 1, length.out = 16),
          each = 16
        )
      }
      patch[16, 16, 4] = 0.5
      patch[9, 9, 1] = -0.1
      list(moon_luminance_array = patch)
    },
    lux_to_radiometric_irradiance = function(lux, efficacy) lux / efficacy,
    compute_K_eff = function(...) 100,
    compute_spd_rgb_unit = function(...) c(0.4, 0.35, 0.25),
    calculate_sky_values = function(...) matrix(1, 1, 3)
  )
  disk = generate_moon_disk(
    as.POSIXct("2026-01-28", tz = "UTC"),
    0,
    0,
    resolution = 16,
    earthshine = FALSE,
    moon_extinction_kV = 0,
    altitude = 100,
    number_cores = 2
  )
  expect_identical(getOption("cores"), old_cores)
  expect_equal(seen$cores, 2)
  expect_equal(seen$elev_m, 100)
  expect_equal(seen$width, 32)
  expect_false(seen$earthshine)
  expect_identical(dim(disk$image), c(16L, 16L, 3L))
  expect_true(all(disk$image >= 0))
  expect_equal(disk$elevation_deg, 90)
  expect_equal(disk$image[8, 8, ] / disk$image[7, 8, ], rep(0.5, 3))
  expect_gt(mean(disk$image[, 13:16, 1]), 4 * mean(disk$image[, 1:4, 1]))
  uv = outer(
    2 * ((1:16 - 0.5) / 16) - 1,
    2 * ((1:16 - 0.5) / 16) - 1,
    function(x, y) x^2 + y^2
  )
  t2 = tan(diameter * pi / 360)^2
  weights = 4 * t2 / (1 + t2 * uv)^1.5 / 256
  weights[uv > 1] = 0
  expect_equal(
    vapply(1:3, function(c) sum(disk$image[,, c] * weights), numeric(1)),
    0.003 * c(0.4, 0.35, 0.25)
  )
  local_mocked_bindings(generate_moon_image_latlong = function(...) {
    stop("raster failure")
  })
  expect_error(
    generate_moon_disk(
      as.POSIXct("2026-01-28", tz = "UTC"),
      0,
      0,
      resolution = 16,
      number_cores = 2
    ),
    "raster failure"
  )
  expect_identical(getOption("cores"), old_cores)
})


test_that("Moon disks retain continuous radiance and tint across center moonset", {
  elevation = 0
  seen = NULL
  local_mocked_bindings(
    swe_dirs_topo_moon_sun = function(...) {
      list(
        moon_dir_topo = c(0, -cospi(elevation / 180), sinpi(elevation / 180)),
        moon_diameter_degrees = 0.55,
        moon_brightness_lux_unattenuated = 0.3,
        # Deliberately unusable: disk preparation must not use point-source brightness.
        moon_brightness_lux = 0
      )
    },
    generate_moon_image_latlong = function(...) {
      list(
        moon_luminance_array = array(1, c(16, 16, 4))
      )
    },
    calculate_sky_values = function(...) {
      seen <<- list(...)
      matrix(c(4, 2, 1), 1, 3)
    }
  )
  generate = function(degrees) {
    elevation <<- degrees
    generate_moon_disk(
      as.POSIXct("2026-01-04", tz = "UTC"),
      0,
      0,
      resolution = 16
    )
  }
  horizon = generate(0)
  expect_true(all(is.finite(horizon$image)))
  expect_gt(sum(horizon$image), 0)
  expect_equal(seen$theta, 0)
  expect_equal(seen$elevation, 0)
  above = generate(1e-8)
  below = generate(-1e-8)
  expect_equal(above$image, below$image, tolerance = 1e-6)
  half_set = generate(-0.14)
  expect_equal(half_set$elevation_deg, -0.14)
  expect_identical(half_set$image, horizon$image)
  expect_gt(sum(half_set$image[9:16, , ]), 0) # The texture itself is not horizon-masked.
  expect_equal(seen$theta, 0)
  expect_equal(seen$elevation, 0)
  fully_set = generate(-1)
  expect_identical(fully_set$image, horizon$image) # Consumers own visibility.
})

test_that("extended-source attenuation has a finite horizon limit", {
  k = 0.172
  expected = 0.3 * 10^(-0.4 * k * 40)
  for (elevation in c(0, -0.1, -1)) {
    expect_equal(apply_airmass_extinction(0.3, elevation, k), expected)
    expect_equal(
      apply_airmass_extinction(0.3, elevation, k, clip_horizon = TRUE),
      0
    )
    expect_equal(
      apply_airmass_extinction(0.3, elevation, k, clip_horizon = FALSE),
      expected
    )
    expect_equal(
      apply_airmass_extinction(0.3, elevation, 0, clip_horizon = FALSE),
      0.3
    )
  }
  expect_equal(
    apply_airmass_extinction(0.3, 30, k, clip_horizon = FALSE),
    apply_airmass_extinction(0.3, 30, k)
  )
  expect_equal(
    apply_airmass_extinction(0.3, 1e-8, k, clip_horizon = FALSE),
    expected,
    tolerance = 1e-6
  )
})
