#' Attach skymodelr RGB color metadata
#'
#' @param image Image array or rayimage image.
#' @keywords internal
as_sky_image = function(image) {
  l_band = attr(image, "L_band", exact = TRUE)
  exr_metadata = attr(image, "exr", exact = TRUE)
  colorspace = attr(image, "colorspace", exact = TRUE)
  white_current = attr(image, "white_current", exact = TRUE)
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
    normalize = FALSE,
    assume_white = if (is.null(white_current)) "D65" else white_current,
    assume_colorspace = rayimage::CS_SRGB
  )
  if (!is.null(l_band)) {
    attr(image, "L_band") = l_band
  }
  if (!is.null(colorspace)) {
    attr(image, "colorspace") = colorspace
  }
  if (!is.null(white_current)) {
    attr(image, "white_current") = white_current
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

.skymodelr_d60_xy = c(0.32168, 0.33767)

.skymodelr_exr_rgb_encoding_label =
  "linear sRGB / Rec.709, D65 chromaticities"
.skymodelr_exr_adopted_neutral_label = "D60"

#' Convert xy chromaticity to XYZ with Y equal to 1
#'
#' @param xy Finite xy chromaticity vector.
#' @keywords internal
xy_to_xyz_y1 = function(xy) {
  if (!is.numeric(xy) || length(xy) != 2L || any(!is.finite(xy))) {
    stop("xy must be a finite numeric vector of length 2.", call. = FALSE)
  }

  x = xy[1]
  y = xy[2]

  if (y <= 0) {
    stop("xy y-coordinate must be positive.", call. = FALSE)
  }

  c(X = x / y, Y = 1, Z = (1 - x - y) / y)
}

#' Get skymodelr EXR adopted white
#'
#' @param exr_adopted_white Adopted white name or XYZ vector.
#' @keywords internal
get_skymodelr_adopted_white = function(exr_adopted_white) {
  if (is.character(exr_adopted_white) && length(exr_adopted_white) == 1L) {
    white_name = toupper(exr_adopted_white)

    if (white_name == "D60") {
      return(list(
        name = "D60",
        xy = .skymodelr_d60_xy,
        xyz = xy_to_xyz_y1(.skymodelr_d60_xy)
      ))
    }

    if (white_name == "D65") {
      xy = c(0.3127, 0.3290)
      return(list(
        name = "D65",
        xy = xy,
        xyz = xy_to_xyz_y1(xy)
      ))
    }

    stop(
      "exr_adopted_white must be \"D60\", \"D65\", or numeric XYZ with Y = 1.",
      call. = FALSE
    )
  }

  if (is.numeric(exr_adopted_white) && length(exr_adopted_white) == 3L) {
    if (any(!is.finite(exr_adopted_white)) || exr_adopted_white[2] <= 0) {
      stop(
        "Numeric exr_adopted_white must be finite XYZ with positive Y.",
        call. = FALSE
      )
    }

    xyz = as.numeric(exr_adopted_white)
    xyz = xyz / xyz[2]
    sum_xyz = sum(xyz)

    if (!is.finite(sum_xyz) || sum_xyz <= 0) {
      stop(
        "Numeric exr_adopted_white must have positive XYZ sum.",
        call. = FALSE
      )
    }

    xy = c(xyz[1] / sum_xyz, xyz[2] / sum_xyz)

    return(list(
      name = "custom",
      xy = xy,
      xyz = xyz
    ))
  }

  stop(
    "exr_adopted_white must be \"D60\", \"D65\", or numeric XYZ with Y = 1.",
    call. = FALSE
  )
}

#' Tag skymodelr EXR metadata
#'
#' @param sky Sky image array.
#' @param adopted_white_xy Default `.skymodelr_d60_xy`. Adopted white xy.
#' @param adopted_white_name Default `.skymodelr_exr_adopted_neutral_label`.
#'   Adopted white name.
#' @param rgb_colorspace Default `rayimage::CS_SRGB`. RGB colorspace.
#' @param model_name Default `NULL`. Sky model name.
#' @param prague_rgb_correction Default `NULL`. Prague RGB correction label.
#' @keywords internal
tag_skymodelr_exr_metadata = function(
  sky,
  adopted_white_xy = .skymodelr_d60_xy,
  adopted_white_name = .skymodelr_exr_adopted_neutral_label,
  rgb_colorspace = rayimage::CS_SRGB,
  model_name = NULL,
  prague_rgb_correction = NULL
) {
  l_band = attr(sky, "L_band", exact = TRUE)
  exr = attr(sky, "exr", exact = TRUE)
  prague_rgb_correction_attr = attr(
    sky,
    "prague_rgb_correction",
    exact = TRUE
  )
  prague_rgb_correction_gain = attr(
    sky,
    "prague_rgb_correction_gain",
    exact = TRUE
  )
  prague_rgb_correction_strength = attr(
    sky,
    "prague_rgb_correction_strength",
    exact = TRUE
  )

  out = rayimage::ray_read_image(
    sky,
    normalize = FALSE,
    assume_colorspace = rgb_colorspace,
    assume_white = xy_to_xyz_y1(adopted_white_xy)
  )

  if (!is.null(l_band)) {
    attr(out, "L_band") = l_band
  }
  if (!is.null(prague_rgb_correction_attr)) {
    attr(out, "prague_rgb_correction") = prague_rgb_correction_attr
  }
  if (!is.null(prague_rgb_correction_gain)) {
    attr(out, "prague_rgb_correction_gain") = prague_rgb_correction_gain
  }
  if (!is.null(prague_rgb_correction_strength)) {
    attr(out, "prague_rgb_correction_strength") =
      prague_rgb_correction_strength
  }

  if (!is.list(exr)) {
    exr = list()
  }

  exr$skymodelr_rgb_encoding = .skymodelr_exr_rgb_encoding_label
  exr$skymodelr_adopted_neutral = adopted_white_name
  exr$skymodelr_adopted_neutral_xy = sprintf(
    "%.8f %.8f",
    adopted_white_xy[1],
    adopted_white_xy[2]
  )

  if (!is.null(model_name)) {
    exr$skymodelr_model = as.character(model_name)
  }

  if (!is.null(prague_rgb_correction)) {
    exr$skymodelr_prague_rgb_correction =
      as.character(prague_rgb_correction)
  }

  attr(out, "exr") = exr

  if (inherits(sky, "sky_image") && !inherits(out, "sky_image")) {
    class(out) = unique(c("sky_image", class(out)))
  }

  out
}

#' Tag generated sky EXR metadata
#'
#' @param sky Sky image array.
#' @param hosek Whether the Hosek-Wilkie model was used.
#' @param exr_metadata Whether to attach skymodelr EXR metadata.
#' @param exr_adopted_white Adopted white name or XYZ vector.
#' @keywords internal
tag_generated_sky_exr_metadata = function(
  sky,
  hosek,
  exr_metadata,
  exr_adopted_white
) {
  if (!isTRUE(exr_metadata)) {
    return(sky)
  }

  adopted_white = get_skymodelr_adopted_white(exr_adopted_white)

  model_name = if (isTRUE(hosek)) {
    "Hosek-Wilkie"
  } else {
    "Prague"
  }

  prague_correction_label = NULL

  if (!isTRUE(hosek)) {
    prague_correction_label = attr(
      sky,
      "prague_rgb_correction",
      exact = TRUE
    )

    if (is.null(prague_correction_label)) {
      prague_correction_label = "none"
    }
  }

  tag_skymodelr_exr_metadata(
    sky,
    adopted_white_xy = adopted_white$xy,
    adopted_white_name = adopted_white$name,
    rgb_colorspace = rayimage::CS_SRGB,
    model_name = model_name,
    prague_rgb_correction = prague_correction_label
  )
}

#' Build EXR metadata for skymodelr RGB output
#'
#' @param image Default `NULL`. Image to read colorspace and white metadata from.
#' @keywords internal
sky_exr_metadata = function(image = NULL) {
  colorspace = attr(image, "colorspace", exact = TRUE)
  white_current = attr(image, "white_current", exact = TRUE)
  exr_metadata = attr(image, "exr", exact = TRUE)

  if (is.null(colorspace)) {
    colorspace = rayimage::CS_SRGB
  }
  if (is.null(white_current)) {
    white_current = colorspace$white_xyz
  }
  if (!is.list(exr_metadata)) {
    exr_metadata = list()
  }

  colorspace_white_xyz = colorspace$white_xyz
  colorspace_white_xy = unname(c(
    colorspace_white_xyz[1] / sum(colorspace_white_xyz),
    colorspace_white_xyz[2] / sum(colorspace_white_xyz)
  ))
  adopted_white_xy = unname(c(
    white_current[1] / sum(white_current),
    white_current[2] / sum(white_current)
  ))

  utils::modifyList(
    exr_metadata,
    list(
      chromaticities = list(
        red = unname(colorspace$primaries$r),
        green = unname(colorspace$primaries$g),
        blue = unname(colorspace$primaries$b),
        white = colorspace_white_xy
      ),
      adoptedNeutral = adopted_white_xy,
      envmap = "latlong"
    )
  )
}

#' Keep EXR metadata fields supported by libopenexr
#'
#' @param metadata EXR metadata list.
#' @keywords internal
filter_supported_exr_metadata = function(metadata) {
  metadata[intersect(
    names(metadata),
    c("chromaticities", "adoptedNeutral", "whiteLuminance", "envmap")
  )]
}

#' Write a skymodelr image with EXR color metadata
#'
#' @param image Image array or rayimage image.
#' @param filename Destination image path.
#' @param ... Additional arguments passed to [rayimage::ray_write_image()].
#' @keywords internal
write_sky_image = function(image, filename, ...) {
  dots = list(...)
  image_to_write = image
  if (tolower(tools::file_ext(filename)) == "exr") {
    metadata = dots$metadata
    if (!is.null(metadata) && !is.list(metadata)) {
      stop("EXR metadata must be a list.", call. = FALSE)
    }
    dots$metadata = filter_supported_exr_metadata(
      if (is.null(metadata)) {
        sky_exr_metadata(image)
      } else {
        utils::modifyList(sky_exr_metadata(image), metadata)
      }
    )
    attr(image_to_write, "exr") = filter_supported_exr_metadata(
      attr(image_to_write, "exr", exact = TRUE)
    )
  }

  do.call(
    rayimage::ray_write_image,
    c(
      list(image = image_to_write, filename = filename),
      dots
    )
  )
}
