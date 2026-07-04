#' Attach skymodelr RGB color metadata
#'
#' @param image Image array or rayimage image.
#' @keywords internal
as_sky_image = function(image) {
  l_band = attr(image, "L_band", exact = TRUE)
  exr_metadata = attr(image, "exr", exact = TRUE)
  prague_rgb_correction = attr(image, "prague_rgb_correction", exact = TRUE)
  prague_rgb_correction_gain = attr(
    image,
    "prague_rgb_correction_gain",
    exact = TRUE
  )
  prague_rgb_correction_strength = attr(
    image,
    "prague_rgb_correction_strength",
    exact = TRUE
  )
  image = rayimage::ray_read_image(
    image,
    assume_white = "D65",
    assume_colorspace = rayimage::CS_SRGB
  )
  if (!is.null(l_band)) {
    attr(image, "L_band") = l_band
  }
  if (!is.null(prague_rgb_correction)) {
    attr(image, "prague_rgb_correction") = prague_rgb_correction
  }
  if (!is.null(prague_rgb_correction_gain)) {
    attr(image, "prague_rgb_correction_gain") = prague_rgb_correction_gain
  }
  if (!is.null(prague_rgb_correction_strength)) {
    attr(image, "prague_rgb_correction_strength") =
      prague_rgb_correction_strength
  }
  if (!is.list(exr_metadata)) {
    exr_metadata = list()
  }
  attr(image, "exr") = utils::modifyList(
    exr_metadata,
    list(envmap = "latlong")
  )
  image
}

#' Build EXR metadata for skymodelr RGB output
#'
#' @keywords internal
sky_exr_metadata = function() {
  colorspace = rayimage::CS_SRGB
  colorspace_white_xyz = colorspace$white_xyz
  colorspace_white_xy = unname(c(
    colorspace_white_xyz[1] / sum(colorspace_white_xyz),
    colorspace_white_xyz[2] / sum(colorspace_white_xyz)
  ))
  d65_white_xyz = c(0.95047, 1, 1.08883)
  d65_white_xy = unname(c(
    d65_white_xyz[1] / sum(d65_white_xyz),
    d65_white_xyz[2] / sum(d65_white_xyz)
  ))

  list(
    chromaticities = list(
      red = unname(colorspace$primaries$r),
      green = unname(colorspace$primaries$g),
      blue = unname(colorspace$primaries$b),
      white = colorspace_white_xy
    ),
    adoptedNeutral = d65_white_xy,
    envmap = "latlong"
  )
}

#' Write a skymodelr image with EXR color metadata
#'
#' @param image Image array or rayimage image.
#' @param filename Destination image path.
#' @param ... Additional arguments passed to [rayimage::ray_write_image()].
#' @keywords internal
write_sky_image = function(image, filename, ...) {
  dots = list(...)
  if (tolower(tools::file_ext(filename)) == "exr") {
    metadata = dots$metadata
    if (!is.null(metadata) && !is.list(metadata)) {
      stop("EXR metadata must be a list.", call. = FALSE)
    }
    dots$metadata = if (is.null(metadata)) {
      sky_exr_metadata()
    } else {
      utils::modifyList(sky_exr_metadata(), metadata)
    }
  }

  do.call(
    rayimage::ray_write_image,
    c(
      list(image = image, filename = filename),
      dots
    )
  )
}
