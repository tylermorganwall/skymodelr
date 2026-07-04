.prague_rgb_correction_gain = c(
  R = 0.94438727,
  G = 1.02157200,
  B = 0.95012063
)

#' Normalize Prague RGB correction option
#'
#' @param prague_rgb_correction Prague RGB correction option.
#' @keywords internal
normalize_prague_rgb_correction = function(prague_rgb_correction) {
  if (is.logical(prague_rgb_correction)) {
    if (length(prague_rgb_correction) != 1L || is.na(prague_rgb_correction)) {
      stop(
        "prague_rgb_correction must be TRUE, FALSE, \"constant\", or \"none\".",
        call. = FALSE
      )
    }

    return(if (isTRUE(prague_rgb_correction)) "constant" else "none")
  }

  prague_rgb_correction = tolower(as.character(prague_rgb_correction))
  match.arg(prague_rgb_correction, c("none", "constant"))
}

#' Validate Prague RGB correction gain
#'
#' @param gain RGB correction gain.
#' @keywords internal
validate_prague_rgb_gain = function(gain) {
  if (!is.numeric(gain) || length(gain) != 3L) {
    stop(
      "prague_rgb_correction_gain must be a numeric vector of length 3.",
      call. = FALSE
    )
  }

  if (is.null(names(gain))) {
    names(gain) = c("R", "G", "B")
  }

  if (!all(c("R", "G", "B") %in% names(gain))) {
    stop(
      "prague_rgb_correction_gain must have names R, G, and B.",
      call. = FALSE
    )
  }

  gain = gain[c("R", "G", "B")]

  if (any(!is.finite(gain)) || any(gain <= 0)) {
    stop(
      "prague_rgb_correction_gain values must be finite and positive.",
      call. = FALSE
    )
  }

  gain
}

#' Prepare Prague RGB correction gain
#'
#' @param gain Default `"auto"`. RGB correction gain.
#' @param strength Default `1`. Correction strength.
#' @keywords internal
prepare_prague_rgb_gain = function(
  gain = "auto",
  strength = 1
) {
  if (identical(gain, "auto")) {
    gain = .prague_rgb_correction_gain
  }

  gain = validate_prague_rgb_gain(gain)

  if (
    !is.numeric(strength) ||
      length(strength) != 1L ||
      !is.finite(strength) ||
      strength < 0
  ) {
    stop(
      "prague_rgb_correction_strength must be a finite non-negative scalar.",
      call. = FALSE
    )
  }

  exp(log(gain) * strength)
}

#' Apply Prague RGB correction gain
#'
#' @param x RGB array or matrix.
#' @param gain RGB correction gain.
#' @keywords internal
apply_prague_rgb_gain = function(x, gain) {
  gain = validate_prague_rgb_gain(gain)

  out = x
  d = dim(out)

  if (length(d) == 3L) {
    if (d[3] < 3L) {
      stop("RGB array must have at least 3 channels.", call. = FALSE)
    }

    out[,, 1] = out[,, 1] * gain[["R"]]
    out[,, 2] = out[,, 2] * gain[["G"]]
    out[,, 3] = out[,, 3] * gain[["B"]]
    return(out)
  }

  if (is.matrix(out)) {
    if (ncol(out) < 3L) {
      stop("RGB matrix must have at least 3 columns.", call. = FALSE)
    }

    out[, 1] = out[, 1] * gain[["R"]]
    out[, 2] = out[, 2] * gain[["G"]]
    out[, 3] = out[, 3] * gain[["B"]]
    return(out)
  }

  stop("Expected an RGB array or RGB matrix.", call. = FALSE)
}
