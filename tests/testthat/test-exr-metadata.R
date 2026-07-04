skip_if_no_exr_metadata = function() {
  testthat::skip_if_not(
    length(find.package("libopenexr", quiet = TRUE)) > 0,
    "libopenexr is not installed"
  )
  testthat::skip_if_not(
    "metadata" %in% names(formals(libopenexr::write_exr)),
    "installed libopenexr does not support metadata"
  )
}

exr_metadata_prague_ground_dataset = function() {
  file.path(
    tools::R_user_dir("skymodelr", "data"),
    "SkyModelDatasetGround.dat"
  )
}

xy_from_xyz = function(xyz) {
  unname(c(xyz[1] / sum(xyz), xyz[2] / sum(xyz)))
}

expect_srgb_image = function(
  image,
  expected_white_xy,
  expect_srgb_name = TRUE
) {
  testthat::expect_s3_class(image, "rayimg")
  colorspace = attr(image, "colorspace", exact = TRUE)
  if (expect_srgb_name) {
    testthat::expect_equal(colorspace$name, "sRGB")
    testthat::expect_equal(colorspace$white_name, "D65")
  }
  testthat::expect_equal(
    unname(colorspace$primaries$r),
    unname(rayimage::CS_SRGB$primaries$r),
    tolerance = 1e-6
  )
  testthat::expect_equal(
    unname(colorspace$primaries$g),
    unname(rayimage::CS_SRGB$primaries$g),
    tolerance = 1e-6
  )
  testthat::expect_equal(
    unname(colorspace$primaries$b),
    unname(rayimage::CS_SRGB$primaries$b),
    tolerance = 1e-6
  )
  testthat::expect_equal(
    xy_from_xyz(attr(image, "white_current", exact = TRUE)),
    expected_white_xy,
    tolerance = 1e-6
  )
  testthat::expect_equal(
    attr(image, "exr", exact = TRUE)$envmap,
    "latlong"
  )
}

expect_srgb_metadata = function(metadata, expected_white_xy) {
  expected = sky_exr_metadata(
    tag_skymodelr_exr_metadata(
      array(0, dim = c(1, 1, 4)),
      adopted_white_xy = expected_white_xy
    )
  )
  testthat::expect_equal(
    unname(metadata$chromaticities$red),
    expected$chromaticities$red,
    tolerance = 1e-6
  )
  testthat::expect_equal(
    unname(metadata$chromaticities$green),
    expected$chromaticities$green,
    tolerance = 1e-6
  )
  testthat::expect_equal(
    unname(metadata$chromaticities$blue),
    expected$chromaticities$blue,
    tolerance = 1e-6
  )
  testthat::expect_equal(
    unname(metadata$chromaticities$white),
    expected$chromaticities$white,
    tolerance = 1e-6
  )
  testthat::expect_equal(
    unname(metadata$adoptedNeutral),
    expected_white_xy,
    tolerance = 1e-6
  )
  testthat::expect_equal(metadata$envmap, "latlong")
}

test_that("generate_sky tags returned Prague sky with sRGB colorspace and D60 adopted neutral", {
  testthat::skip_if_not(file.exists(exr_metadata_prague_ground_dataset()))

  sky = generate_sky(
    resolution = 16,
    elevation = 60,
    azimuth = 315,
    hosek = FALSE,
    visibility = 131.8,
    albedo = 0.2,
    render_mode = "all",
    exr_metadata = TRUE,
    exr_adopted_white = "D60"
  )

  expected_d60 = xy_to_xyz_y1(c(0.32168, 0.33767))
  cs = attr(sky, "colorspace", exact = TRUE)
  white_current = attr(sky, "white_current", exact = TRUE)
  exr = attr(sky, "exr", exact = TRUE)

  expect_s3_class(sky, "rayimg")
  expect_equal(cs$name, "sRGB")
  expect_equal(cs$white_name, "D65")
  expect_equal(cs$primaries$r, c(0.64, 0.33), tolerance = 1e-6)
  expect_equal(cs$primaries$g, c(0.30, 0.60), tolerance = 1e-6)
  expect_equal(cs$primaries$b, c(0.15, 0.06), tolerance = 1e-6)
  expect_equal(
    as.numeric(white_current),
    as.numeric(expected_d60),
    tolerance = 1e-8
  )

  expect_true(is.list(exr))
  expect_equal(exr$skymodelr_adopted_neutral, "D60")
  expect_equal(
    exr$skymodelr_rgb_encoding,
    "linear sRGB / Rec.709, D65 chromaticities"
  )
})

test_that("generate_sky can tag D65 adopted neutral", {
  testthat::skip_if_not(file.exists(exr_metadata_prague_ground_dataset()))

  sky = generate_sky(
    resolution = 16,
    elevation = 60,
    azimuth = 315,
    hosek = FALSE,
    exr_metadata = TRUE,
    exr_adopted_white = "D65"
  )

  expected_d65 = rayimage::CS_SRGB$white_xyz
  expect_equal(
    as.numeric(attr(sky, "white_current", exact = TRUE)),
    as.numeric(expected_d65),
    tolerance = 1e-8
  )
  expect_equal(attr(sky, "exr", exact = TRUE)$skymodelr_adopted_neutral, "D65")
})

test_that("generate_sky can disable EXR metadata tagging", {
  testthat::skip_if_not(file.exists(exr_metadata_prague_ground_dataset()))

  sky = generate_sky(
    resolution = 16,
    elevation = 60,
    azimuth = 315,
    hosek = FALSE,
    exr_metadata = FALSE
  )

  exr = attr(sky, "exr", exact = TRUE)
  expect_null(exr$skymodelr_adopted_neutral)
  expect_null(exr$skymodelr_rgb_encoding)
})

test_that("EXR metadata tagging does not change pixels", {
  testthat::skip_if_not(file.exists(exr_metadata_prague_ground_dataset()))

  untagged = generate_sky(
    resolution = 16,
    elevation = 60,
    azimuth = 315,
    hosek = FALSE,
    visibility = 131.8,
    albedo = 0.2,
    render_mode = "all",
    exr_metadata = FALSE
  )

  tagged = generate_sky(
    resolution = 16,
    elevation = 60,
    azimuth = 315,
    hosek = FALSE,
    visibility = 131.8,
    albedo = 0.2,
    render_mode = "all",
    exr_metadata = TRUE,
    exr_adopted_white = "D60"
  )

  expect_equal(as.numeric(tagged), as.numeric(untagged), tolerance = 0)
  expect_equal(attr(tagged, "L_band"), attr(untagged, "L_band"))
})

test_that("generate_sky EXRs carry sRGB/D60 metadata matching memory output", {
  skip_if_no_exr_metadata()

  cases = list(
    list(elevation = 15, render_mode = "atmosphere"),
    list(elevation = -1, render_mode = "all")
  )
  expected_d60_xy = .skymodelr_d60_xy

  for (case in cases) {
    path = tempfile(fileext = ".exr")
    memory_sky = generate_sky(
      resolution = 8,
      elevation = case$elevation,
      azimuth = 135,
      render_mode = case$render_mode
    )
    file_sky = generate_sky(
      filename = path,
      resolution = 8,
      elevation = case$elevation,
      azimuth = 135,
      render_mode = case$render_mode
    )

    expect_srgb_image(memory_sky, expected_d60_xy)
    expect_srgb_image(file_sky, expected_d60_xy)
    testthat::expect_equal(
      as.numeric(file_sky),
      as.numeric(memory_sky),
      tolerance = 1e-12
    )

    file_metadata = libopenexr::read_exr(path)$metadata
    read_sky = rayimage::ray_read_image(path, normalize = FALSE)
    expect_srgb_image(read_sky, expected_d60_xy, expect_srgb_name = FALSE)
    expect_srgb_metadata(file_metadata, expected_d60_xy)
    testthat::expect_null(
      attr(read_sky, "exr", exact = TRUE)$skymodelr_adopted_neutral
    )
  }
})

test_that("EXR round trip preserves D60 adopted neutral metadata", {
  testthat::skip_if_not(file.exists(exr_metadata_prague_ground_dataset()))
  skip_if_no_exr_metadata()

  filename = tempfile(fileext = ".exr")

  generate_sky(
    filename = filename,
    resolution = 16,
    elevation = 60,
    azimuth = 315,
    hosek = FALSE,
    visibility = 131.8,
    albedo = 0.2,
    render_mode = "all",
    exr_metadata = TRUE,
    exr_adopted_white = "D60"
  )

  img = rayimage::ray_read_image(filename, normalize = FALSE)

  cs = attr(img, "colorspace", exact = TRUE)
  white_current = attr(img, "white_current", exact = TRUE)
  exr = attr(img, "exr", exact = TRUE)

  expect_equal(cs$primaries$r, c(0.64, 0.33), tolerance = 1e-6)
  expect_equal(cs$primaries$g, c(0.30, 0.60), tolerance = 1e-6)
  expect_equal(cs$primaries$b, c(0.15, 0.06), tolerance = 1e-6)

  expected_d60 = xy_to_xyz_y1(c(0.32168, 0.33767))
  expect_equal(
    as.numeric(white_current),
    as.numeric(expected_d60),
    tolerance = 1e-6
  )

  expect_true(is.list(exr))
  expect_null(exr$skymodelr_adopted_neutral)
})

test_that("generate_sky_latlong EXRs match memory output metadata", {
  skip_if_no_exr_metadata()

  args = list(
    datetime = as.POSIXct("2025-03-21 12:00:00", tz = "America/New_York"),
    lat = 38.9072,
    lon = -77.0369,
    resolution = 8,
    number_cores = 1
  )
  expected_d60_xy = .skymodelr_d60_xy
  path = tempfile(fileext = ".exr")
  memory_sky = do.call(generate_sky_latlong, args)
  file_sky = do.call(
    generate_sky_latlong,
    c(args, list(filename = path))
  )

  expect_srgb_image(memory_sky, expected_d60_xy)
  expect_srgb_image(file_sky, expected_d60_xy)
  testthat::expect_equal(
    as.numeric(file_sky),
    as.numeric(memory_sky),
    tolerance = 1e-12
  )

  file_metadata = libopenexr::read_exr(path)$metadata
  read_sky = rayimage::ray_read_image(path, normalize = FALSE)
  expect_srgb_image(read_sky, expected_d60_xy, expect_srgb_name = FALSE)
  expect_srgb_metadata(file_metadata, expected_d60_xy)
  testthat::expect_null(
    attr(read_sky, "exr", exact = TRUE)$skymodelr_adopted_neutral
  )
})
