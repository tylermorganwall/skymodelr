#' Write a Hosek–Wilkie sky dome to EXR
#'
#' @param outfile            Default `""`. Path to `.exr` file.
#' @param albedo             Default `0.5`. 0–1 ground albedo.
#' @param turbidity          Default `3`. 1.7–10 atmospheric turbidity.
#' @param elevation          Default `10`. Solar elevation above the horizon (°).
#' @param resolution         Default `2048`. Linear resolution of the square map (≥ 1).
#' @param numbercores        Default `1`. Threads for \code{RcppThread}.
#' @param square_projection  Default `FALSE`. If \code{TRUE} use equal‑area square mapping,
#'                           else latitude–longitude square mapping.
#'
#' @return Invisible `NULL`.  The EXR is written to `outfile`.
#' @export
#' @examples
#' makesky("sky.exr", turbidity = 2.4, elevation = 6, numbercores = 1)
makesky = function(
  outfile,
  albedo = 0.5,
  turbidity = 3,
  elevation = 10,
  resolution = 2048,
  numbercores = 1,
  square_projection = FALSE
) {
  makesky_rcpp(
    outfile,
    albedo,
    turbidity,
    elevation,
    resolution,
    numbercores,
    square_projection
  )
  invisible(NULL)
}
