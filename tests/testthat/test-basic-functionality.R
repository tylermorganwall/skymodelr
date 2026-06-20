test_that("basic Hosek brightness calculation works without downloads", {
  result = calculate_sun_brightness(
    elevation = 45,
    hosek = TRUE,
    lambda_nm = seq(380, 720, by = 80)
  )

  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(is.finite(result))
  expect_gt(result, 0)
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

  outside_file = tempfile()
  writeBin(charToRaw("outside"), outside_file)
  expect_message(
    clear_sky_data(files = outside_file, ask = FALSE),
    "No matching cached sky data files found.",
    fixed = TRUE
  )
  expect_true(file.exists(outside_file))
  expect_true(file.exists(cached_file))

  removed = clear_sky_data(
    files = file.path("..", "SkyModelDatasetGround.dat"),
    ask = FALSE
  )
  expect_false(file.exists(cached_file))
  expect_equal(basename(removed), "SkyModelDatasetGround.dat")
})
