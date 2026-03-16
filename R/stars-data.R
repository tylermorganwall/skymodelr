#' Star catalog used for sky rendering
#'
#' @description
#' A data frame of stellar positions and photometric/color information used by
#' [generate_stars()] to render star fields.
#'
#' @format
#' A tibble with 9,110 rows and 8 variables:
#' \describe{
#'   \item{bsc_number}{Numeric identifier from the source catalog.}
#'   \item{ra_rad}{Right ascension in radians.}
#'   \item{dec_rad}{Declination in radians.}
#'   \item{v_mag}{Apparent visual magnitude (V band).}
#'   \item{spec}{Spectral type string.}
#'   \item{r}{Relative red channel weight derived from spectral type.}
#'   \item{g}{Relative green channel weight derived from spectral type.}
#'   \item{b}{Relative blue channel weight derived from spectral type.}
#' }
#'
#' @source
#' Bright Star Catalogue, 5th Revised Edition (BSC5), stored in `data-raw/BSC5`.
#'
#' @usage data(stars)
#' @keywords datasets internal
"stars"
