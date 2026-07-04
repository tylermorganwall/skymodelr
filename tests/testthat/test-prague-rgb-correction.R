testthat::skip_on_cran()

prague_rgb_correction_ground_dataset = function() {
  file.path(
    tools::R_user_dir("skymodelr", "data"),
    "SkyModelDatasetGround.dat"
  )
}

test_that("Prague RGB correction is multiplicative for generate_sky", {
  testthat::skip_if_not(file.exists(prague_rgb_correction_ground_dataset()))

  raw = generate_sky(
    resolution = 16,
    elevation = 60,
    azimuth = 315,
    albedo = 0.2,
    hosek = FALSE,
    visibility = 131.8,
    render_mode = "all",
    prague_rgb_correction = FALSE
  )

  corrected = generate_sky(
    resolution = 16,
    elevation = 60,
    azimuth = 315,
    albedo = 0.2,
    hosek = FALSE,
    visibility = 131.8,
    render_mode = "all",
    prague_rgb_correction = TRUE
  )
  default_corrected = generate_sky(
    resolution = 16,
    elevation = 60,
    azimuth = 315,
    albedo = 0.2,
    hosek = FALSE,
    visibility = 131.8,
    render_mode = "all"
  )

  gain = attr(corrected, "prague_rgb_correction_gain")

  expect_equal(default_corrected, corrected)
  expect_equal(as.numeric(corrected[,, 1]), as.numeric(raw[,, 1] * gain[["R"]]))
  expect_equal(as.numeric(corrected[,, 2]), as.numeric(raw[,, 2] * gain[["G"]]))
  expect_equal(as.numeric(corrected[,, 3]), as.numeric(raw[,, 3] * gain[["B"]]))
  expect_equal(as.numeric(corrected[,, 4]), as.numeric(raw[,, 4]))
  expect_equal(attr(corrected, "L_band"), attr(raw, "L_band"))

  expect_equal(attr(corrected, "prague_rgb_correction"), "constant")
  expect_equal(attr(corrected, "prague_rgb_correction_strength"), 1)
})

test_that("Prague RGB correction strength interpolates logarithmically", {
  testthat::skip_if_not(file.exists(prague_rgb_correction_ground_dataset()))

  raw = generate_sky(
    resolution = 16,
    elevation = 60,
    azimuth = 315,
    albedo = 0.2,
    hosek = FALSE,
    visibility = 131.8,
    render_mode = "all",
    prague_rgb_correction = FALSE
  )

  half = generate_sky(
    resolution = 16,
    elevation = 60,
    azimuth = 315,
    albedo = 0.2,
    hosek = FALSE,
    visibility = 131.8,
    render_mode = "all",
    prague_rgb_correction = TRUE,
    prague_rgb_correction_strength = 0.5
  )

  gain = attr(half, "prague_rgb_correction_gain")

  expect_equal(as.numeric(half[,, 1]), as.numeric(raw[,, 1] * gain[["R"]]))
  expect_equal(as.numeric(half[,, 2]), as.numeric(raw[,, 2] * gain[["G"]]))
  expect_equal(as.numeric(half[,, 3]), as.numeric(raw[,, 3] * gain[["B"]]))
  expect_equal(as.numeric(half[,, 4]), as.numeric(raw[,, 4]))
})

test_that("Prague RGB correction is ignored for Hosek generate_sky output", {
  raw = generate_sky(
    resolution = 16,
    elevation = 60,
    azimuth = 315,
    hosek = TRUE,
    prague_rgb_correction = FALSE
  )

  corrected = generate_sky(
    resolution = 16,
    elevation = 60,
    azimuth = 315,
    hosek = TRUE,
    prague_rgb_correction = TRUE
  )

  expect_equal(corrected, raw)
  expect_null(attr(corrected, "prague_rgb_correction"))
})

test_that("calculate_sky_values uses the same Prague RGB correction", {
  testthat::skip_if_not(file.exists(prague_rgb_correction_ground_dataset()))

  raw = calculate_sky_values(
    phi = c(90, 180),
    theta = c(45, 60),
    altitude = c(0, 1000),
    elevation = 60,
    visibility = 131.8,
    albedo = 0.2,
    prague_rgb_correction = FALSE
  )

  corrected = calculate_sky_values(
    phi = c(90, 180),
    theta = c(45, 60),
    altitude = c(0, 1000),
    elevation = 60,
    visibility = 131.8,
    albedo = 0.2,
    prague_rgb_correction = TRUE
  )

  gain = attr(corrected, "prague_rgb_correction_gain")

  expect_equal(corrected[, 1], raw[, 1] * gain[["R"]])
  expect_equal(corrected[, 2], raw[, 2] * gain[["G"]])
  expect_equal(corrected[, 3], raw[, 3] * gain[["B"]])
  expect_equal(attr(corrected, "prague_rgb_correction"), "constant")
})

test_that("invalid Prague RGB correction options error clearly", {
  expect_error(
    normalize_prague_rgb_correction(NA),
    "prague_rgb_correction"
  )

  expect_error(
    prepare_prague_rgb_gain(strength = -1),
    "prague_rgb_correction_strength"
  )

  expect_error(
    prepare_prague_rgb_gain(gain = c(R = 1, G = 1)),
    "prague_rgb_correction_gain"
  )

  expect_error(
    prepare_prague_rgb_gain(gain = c(R = 1, G = 1, B = NA)),
    "finite and positive"
  )
})
