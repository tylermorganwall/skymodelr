# ---------------------------------------------------------------------------
# Helpers to convert POSIXct → Julian Date (UTC noon 1 Jan 2000 = 2451545).
#' @noRd
jd_utc = function(time_utc) {
	unclass(time_utc) / 86400 + 2440587.5 # Unix epoch to JD
}

# ---------------------------------------------------------------------------
#' Write a star‑field EXR aligned with `makesky()`
#'
#' @param outfile    `.exr` output path.
#' @param stars      Tibble from `read_bsc5()`.  `NULL` ➜ packaged version.
#' @param resolution Map half‑width (image is 2×`resolution` × `resolution`).
#' @param zero_point Exposure scale (larger = brighter stars).
#' @param lon_deg Observer longitude (east +)
#' @param lat_deg Observer latitude (deg).
#' @param time_utc   `POSIXct` UTC time stamp used to compute sidereal time.
#' @param ncores     Threads for `RcppThread`.
#'
#' @export
makestars = function(
	outfile = "stars.exr",
	stars = NULL,
	resolution = 2048,
	zero_point = 1,
	lon_deg = 0,
	lat_deg = 0,
	time_utc = as.POSIXct("2000-01-01 00:00:00", tz = "UTC"),
	turbidity = 3.0,
	ozone_du = 300.0,
	altitude = 0.0,
	color = TRUE,
	ncores = parallel::detectCores()
) {
	if (!inherits(time_utc, "POSIXct")) {
		stop("time_utc must be POSIXct in UTC")
	}
	attr(time_utc, "tzone") <- "UTC"
	jd = jd_utc(time_utc)

	make_starfield_rcpp(
		outfile = outfile,
		stars = skymodelr::stars,
		resolution = resolution,
		zero_point = zero_point,
		lon_deg = lon_deg,
		lat_deg = lat_deg,
		jd = jd,
		use_rgb = color,
		turbidity = turbidity,
		ozone_du = ozone_du,
		altitude = altitude,
		ncores = ncores
	)
	invisible(NULL)
}
