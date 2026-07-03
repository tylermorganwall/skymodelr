skip_if_no_exr_metadata = function() {
  testthat::skip_if_not(
    length(find.package("libopenexr", quiet = TRUE)) > 0,
    "libopenexr is not installed"
  )
  supports_metadata = getFromNamespace(
    "libopenexr_supports_metadata",
    "rayimage"
  )()
  testthat::skip_if_not(
    supports_metadata,
    "installed libopenexr does not support EXR metadata"
  )
}

xy_from_xyz = function(xyz) {
  unname(c(xyz[1] / sum(xyz), xyz[2] / sum(xyz)))
}

expect_srgb_d65_image = function(image, expect_srgb_name = TRUE) {
  testthat::expect_s3_class(image, "rayimg")
  colorspace = attr(image, "colorspace", exact = TRUE)
  if (expect_srgb_name) {
    testthat::expect_equal(colorspace$name, "sRGB")
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
    sky_exr_metadata()$adoptedNeutral,
    tolerance = 1e-6
  )
  testthat::expect_equal(
    attr(image, "exr", exact = TRUE)$envmap,
    "latlong"
  )
}

expect_srgb_d65_metadata = function(metadata) {
  expected = sky_exr_metadata()
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
    expected$adoptedNeutral,
    tolerance = 1e-6
  )
  testthat::expect_equal(metadata$envmap, "latlong")
}

test_that("generate_sky EXRs carry sRGB/D65 metadata matching memory output", {
  skip_if_no_exr_metadata()

  cases = list(
    list(elevation = 15, render_mode = "atmosphere"),
    list(elevation = -1, render_mode = "all")
  )

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

    expect_srgb_d65_image(memory_sky)
    expect_srgb_d65_image(file_sky)
    testthat::expect_equal(
      as.numeric(file_sky),
      as.numeric(memory_sky),
      tolerance = 1e-12
    )

    file_metadata = libopenexr::read_exr(path)$metadata
    read_sky = rayimage::ray_read_image(path)
    expect_srgb_d65_image(read_sky, expect_srgb_name = FALSE)
    expect_srgb_d65_metadata(file_metadata)
    testthat::expect_equal(
      attr(read_sky, "exr", exact = TRUE)$envmap,
      "latlong"
    )
  }
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
  path = tempfile(fileext = ".exr")
  memory_sky = do.call(generate_sky_latlong, args)
  file_sky = do.call(
    generate_sky_latlong,
    c(args, list(filename = path))
  )

  expect_srgb_d65_image(memory_sky)
  expect_srgb_d65_image(file_sky)
  testthat::expect_equal(
    as.numeric(file_sky),
    as.numeric(memory_sky),
    tolerance = 1e-12
  )

  file_metadata = libopenexr::read_exr(path)$metadata
  read_sky = rayimage::ray_read_image(path)
  expect_srgb_d65_image(read_sky, expect_srgb_name = FALSE)
  expect_srgb_d65_metadata(file_metadata)
  testthat::expect_equal(
    attr(read_sky, "exr", exact = TRUE)$envmap,
    "latlong"
  )
})
