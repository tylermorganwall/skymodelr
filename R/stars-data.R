#' Star catalog used for sky rendering
#'
#' @description
#' A data frame of stellar positions and photometric/color information used by
#' [generate_stars()] to render star fields.
#'
#' @format
#' A data.frame with 9,110 rows and 8 variables:
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
#' Derived from the Bright Star Catalogue, 5th Revised Edition (BSC5), also
#' known as the Yale Bright Star Catalogue. The source catalogue is commonly
#' identified as Hoffleit, D. and Warren, W. H. Jr. (1991), "The Bright Star
#' Catalogue, 5th Revised Ed.", Yale University Observatory, and is distributed
#' as machine-readable catalogue V/50 by the Centre de Donnees astronomiques de
#' Strasbourg (CDS): <https://cdsarc.cds.unistra.fr/viz-bin/cat/V/50>.
#'
#' @usage data(stars)
#' @keywords datasets internal
"stars"
