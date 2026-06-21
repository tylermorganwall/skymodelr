test_that("basic star generation works without downloads", {
  result = generate_stars(
    resolution = 8,
    lon = -77.0369,
    lat = 38.9072,
    datetime = as.POSIXct("2025-03-21 02:20:00", tz = "UTC"),
    atmosphere_effects = FALSE,
    upper_hemisphere_only = TRUE,
    number_cores = 1
  )

  expect_true(is.array(result))
  expect_equal(dim(result), c(8, 16, 4))
  expect_true(all(is.finite(result)))
})

test_that("sky data cache helpers are programmatic and path-safe", {
  cache_dir = tempfile()
  dir.create(cache_dir)
  testthat::local_mocked_bindings(
    prague_coef_cache_dir = function() cache_dir,
    .package = "skymodelr"
  )

  expect_equal(nrow(list_sky_data()), 0)

  cached_file = file.path(cache_dir, "SkyModelDatasetGround.dat")
  writeBin(charToRaw("cached"), cached_file)

  listed = list_sky_data()
  expect_s3_class(listed, "data.frame")
  expect_named(listed, c("file", "path", "size", "modified"))
  expect_equal(listed$file, "SkyModelDatasetGround.dat")
  expect_equal(listed$size, file.info(cached_file)$size)

  missing = NULL
  expect_message(
    missing <- clear_sky_data(files = "not-cached.dat", ask = FALSE),
    "No matching cached sky data files found.",
    fixed = TRUE
  )
  expect_identical(missing, character())
  expect_true(file.exists(cached_file))

  outside_dir = tempfile()
  dir.create(outside_dir)
  outside_file = file.path(outside_dir, "not-in-cache.dat")
  writeBin(charToRaw("outside"), outside_file)
  expect_message(
    clear_sky_data(files = outside_file, ask = FALSE),
    "No matching cached sky data files found.",
    fixed = TRUE
  )
  expect_true(file.exists(outside_file))
  expect_true(file.exists(cached_file))

  outside_same_name = file.path(outside_dir, basename(cached_file))
  writeBin(charToRaw("outside"), outside_same_name)
  removed = clear_sky_data(
    files = file.path("..", basename(cached_file)),
    ask = FALSE
  )
  expect_false(file.exists(cached_file))
  expect_true(file.exists(outside_file))
  expect_true(file.exists(outside_same_name))
  expect_equal(basename(removed), "SkyModelDatasetGround.dat")
})
