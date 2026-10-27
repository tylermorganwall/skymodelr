prague_ground_dataset = function() {
  file.path(
    tools::R_user_dir("skymodelr", "data"),
    "SkyModelDatasetGround.dat"
  )
}

test_that("calculate_sky_radiance splits Prague atmosphere and sun", {
  coef_file = prague_ground_dataset()
  testthat::skip_if_not(file.exists(coef_file))

  args = list(
    phi = c(90, 180),
    theta = c(20, 20),
    lambda_nm = 550,
    altitude = 0,
    elevation = 20,
    visibility = 50,
    albedo = 0.2,
    azimuth = 90,
    number_cores = 1
  )

  L_all = do.call(calculate_sky_radiance, c(args, list(render_mode = "all")))
  L_atm = do.call(
    calculate_sky_radiance,
    c(args, list(render_mode = "atmosphere"))
  )
  L_sun = do.call(calculate_sky_radiance, c(args, list(render_mode = "sun")))

  testthat::expect_equal(L_all, L_atm + L_sun, tolerance = 1e-8)
  testthat::expect_true(all(is.finite(L_all)))
  testthat::expect_gt(L_sun[1], 0)
  testthat::expect_gte(L_sun[1], L_sun[2])
})

test_that("radiance rejects invalid inputs before loading coefficients", {
  expect_error(calculate_sky_radiance(90, 20, 550, backend = "bad"), "arg")
  for (name in c(
    "phi",
    "theta",
    "lambda_nm",
    "altitude",
    "elevation",
    "visibility",
    "albedo",
    "azimuth"
  )) {
    args = list(phi = 90, theta = 20, lambda_nm = 550, backend = "gpu")
    args[[name]] = NA_real_
    expect_error(do.call(calculate_sky_radiance, args), "finite numeric")
  }
  expect_error(calculate_sky_radiance("90", 20, 550), "finite numeric")
  expect_error(calculate_sky_radiance(1:3, 1:2, 550), "differing number")
})

test_that("unavailable GPU requests fail explicitly", {
  skip_if(prague_gpu_available(), "Metal device available")
  # No coefficient file is needed to diagnose an unavailable device.
  expect_error(
    calculate_raw_prague_radiance(
      90,
      20,
      550,
      20,
      0.2,
      0,
      50,
      90,
      backend = "gpu"
    ),
    "available Metal device"
  )
})

test_that("Metal radiance agrees with CPU over a grid and all render modes", {
  skip_on_cran()
  skip_if_not(prague_gpu_available(), "Metal device unavailable")
  skip_if_not(
    file.exists(prague_ground_dataset()),
    "Prague ground data missing"
  )
  grid = expand.grid(
    phi = c(0, 90, 180, 360),
    theta = c(-30, 0, 20, 60, 90),
    lambda_nm = c(320, 540 - 1e-8, 540, 550, 740),
    elevation = c(-4.2, 0, 20, 90),
    visibility = c(20, 80, 131.8),
    albedo = c(0, 0.3, 1)
  )
  results = list()
  for (mode in c("all", "atmosphere", "sun")) {
    args = c(as.list(grid), list(render_mode = mode))
    cpu = do.call(calculate_sky_radiance, c(args, list(backend = "cpu")))
    gpu = do.call(calculate_sky_radiance, c(args, list(backend = "gpu")))
    expect_length(gpu, nrow(grid))
    expect_true(all(is.finite(gpu)))
    # Float geometry and interpolation can differ at small radiances. Use a
    # per-sample bound so bright solar samples cannot hide atmospheric errors.
    expect_true(all(abs(gpu - cpu) <= 1e-6 + 2e-3 * abs(cpu)), info = mode)
    results[[mode]] = gpu
  }
  expect_equal(results$all, results$atmosphere + results$sun, tolerance = 1e-6)
  expect_gt(max(results$sun), 0)
  expect_error(calculate_sky_radiance(90, 20, 1, backend = "gpu"), "within")
})

test_that("Metal batches preserve row order across a batch boundary", {
  skip_on_cran()
  skip_if_not(prague_gpu_available(), "Metal device unavailable")
  skip_if_not(
    file.exists(prague_ground_dataset()),
    "Prague ground data missing"
  )
  n = 65539L
  args = list(
    phi = rep(c(90, 180, 270), length.out = n),
    theta = rep(c(20, 45, 80), length.out = n),
    lambda_nm = rep(c(420, 550, 680), length.out = n),
    elevation = 20,
    azimuth = rep(c(90, 180, 270), length.out = n),
    backend = "gpu"
  )
  actual = do.call(calculate_sky_radiance, args)
  reference = calculate_sky_radiance(
    phi = c(90, 180, 270),
    theta = c(20, 45, 80),
    lambda_nm = c(420, 550, 680),
    elevation = 20,
    azimuth = c(90, 180, 270),
    backend = "cpu"
  )
  expected = rep(reference, length.out = n)
  expect_length(actual, n)
  expect_true(all(abs(actual - expected) <= 1e-6 + 2e-3 * abs(expected)))
})

test_that("Metal coefficients refresh when the dataset changes", {
  skip_on_cran()
  skip_if_not(prague_gpu_available(), "Metal device unavailable")
  skip_if_not(
    file.exists(prague_ground_dataset()),
    "Prague ground data missing"
  )
  wide_file = file.path(
    prague_coef_cache_dir(),
    "PragueSkyModelDatasetGroundInfra.dat"
  )
  skip_if_not(file.exists(wide_file), "Prague wide-spectrum data missing")
  args = list(phi = c(90, 180), theta = c(20, 45), lambda_nm = c(550, 680))
  before = do.call(calculate_sky_radiance, c(args, list(backend = "gpu")))
  wide_args = list(
    phi = c(90, 180),
    theta = c(20, 45),
    lambda_nm = c(1000, 2000),
    wide_spectrum = TRUE
  )
  cpu = do.call(calculate_sky_radiance, c(wide_args, list(backend = "cpu")))
  gpu = do.call(calculate_sky_radiance, c(wide_args, list(backend = "gpu")))
  expect_true(all(abs(gpu - cpu) <= 1e-6 + 2e-3 * abs(cpu)))
  after = do.call(calculate_sky_radiance, c(args, list(backend = "gpu")))
  expect_identical(after, before)
})

test_that("Metal radiance supports grids of observer altitudes", {
  skip_on_cran()
  skip_if_not(prague_gpu_available(), "Metal device unavailable")
  full_file = file.path(prague_coef_cache_dir(), "SkyModelDataset.dat")
  skip_if_not(file.exists(full_file), "Prague full-altitude data missing")
  grid = expand.grid(
    phi = c(90, 180),
    theta = c(-20, 0, 20, 90),
    lambda_nm = c(320, 550, 740),
    altitude = c(0, 100, 5000, 15000)
  )
  args = c(as.list(grid), list(elevation = 20))
  cpu = do.call(calculate_sky_radiance, c(args, list(backend = "cpu")))
  gpu = do.call(calculate_sky_radiance, c(args, list(backend = "gpu")))
  expect_true(all(is.finite(gpu)))
  expect_true(all(abs(gpu - cpu) <= 1e-6 + 2e-3 * abs(cpu)))
})
