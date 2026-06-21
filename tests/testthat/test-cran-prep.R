test_that("missing Prague coefficients fail clearly without prompting", {
  cache_dir = tempfile()
  dir.create(cache_dir)
  testthat::local_mocked_bindings(
    prague_coef_cache_dir = function() cache_dir,
    .package = "skymodelr"
  )

  testthat::expect_error(
    resolve_prague_coef_file(allow_download = FALSE),
    "SkyModelDatasetGround[.]dat.*download_sky_data",
    fixed = FALSE
  )
})

test_that("generate_planets passes latitude and longitude by name", {
  captured = new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    swe_dirs_topo_planets_df = function(datetime, lat, lon, elev_m = 0) {
      captured$lat = lat
      captured$lon = lon
      captured$elev_m = elev_m
      data.frame(
        ra_rad = 0,
        dec_rad = 0,
        v_mag = 0,
        r = 1,
        g = 1,
        b = 1
      )
    },
    make_starfield_rcpp = function(
      stars,
      resolution = 2048L,
      zero_point = 1.0,
      lon_deg = 0.0,
      lat_deg = 0.0,
      jd = 2451545.0,
      turbidity = 3.0,
      ozone_du = 300.0,
      altitude = 0.0,
      star_width = 1.0,
      use_rgb = TRUE,
      atmosphere_effects = TRUE,
      upper_hemisphere_only = TRUE,
      number_cores = 1L
    ) {
      captured$lon_deg = lon_deg
      captured$lat_deg = lat_deg
      array(0, dim = c(resolution, resolution * 2, 3))
    },
    .package = "skymodelr"
  )

  generate_planets(
    datetime = as.POSIXct("2025-03-21 02:20:00", tz = "UTC"),
    lon = -77.0369,
    lat = 38.9072,
    altitude = 123,
    resolution = 2,
    atmosphere_effects = FALSE
  )

  testthat::expect_equal(captured$lat, 38.9072)
  testthat::expect_equal(captured$lon, -77.0369)
  testthat::expect_equal(captured$elev_m, 123)
  testthat::expect_equal(captured$lat_deg, 38.9072)
  testthat::expect_equal(captured$lon_deg, -77.0369)
})

test_that("moon texture is bundled", {
  testthat::expect_true(file.exists(get_moon_texture()))
})
