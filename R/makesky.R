#' Write a Hosek–Wilkie sky dome to EXR
#'
#' @param outfile            Default `""`. Path to `.exr` file.
#' @param albedo             Default `0.5`. 0–1 ground albedo.
#' @param turbidity          Default `3`. 1.7–10 atmospheric turbidity.
#' @param elevation          Default `10`. Solar elevation above the horizon (°).
#' @param azimuth            Default `90`. Solar azimuth (°). Defaults South.
#' @param altitude 	         Default `0`. Altitude of the viewer in meters. Valid [0,15000]. Only valid for the
#' Prague model.
#' @param resolution         Default `2048`. Linear resolution of the square map (≥ 1).
#' @param numbercores        Default `1`. Threads for \code{RcppThread}.
#' @param square_projection  Default `FALSE`. If \code{TRUE} use equal‑area square mapping,
#'                           else latitude–longitude square mapping.
#' @param model Default `"hosek"`. Use `"prague"` to enable the Prague 2021-22 spectral sky model.
#' @param prg_dataset Default `""`. Full path to the Prague binary dataset (`*.dat`) when `model = "prague"`.
#' @param visibility Default `50`. Meteorological range in kilometres for Prague model.
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
	azimuth = 90,
	altitude = 0,
	resolution = 2048,
	numbercores = 1,
	model = "hosek",
	prg_dataset = "",
	visibility = 50,
	square_projection = FALSE
) {
	stopifnot(altitude >= 0 && altitude <= 15000)
	makesky_rcpp(
		outfile,
		albedo,
		turbidity,
		elevation,
		azimuth,
		resolution,
		numbercores,
		square_projection,
		visibility = visibility,
		prg_dataset = prg_dataset,
		altitude,
		model = model
	)
	invisible(NULL)
}
