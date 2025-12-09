# ---------------------------------------------------------------------------
# Helpers to convert POSIXct → Julian Date (UTC noon 1 Jan 2000 = 2451545).
#' Convert POSIXct timestamps to Julian Date
#'
#' @description Transform a POSIXct time value into an astronomical Julian Date
#' (days since noon UTC 1 Jan 4713 BCE).
#' @param time_utc POSIXct vector in UTC.
#' @keywords internal
jd_utc = function(time_utc) {
	unclass(time_utc) / 86400 + 2440587.5 # Unix epoch to JD
}

# ---------------------------------------------------------------------------
#' Generate a star‑field array aligned with `generate_sky()`
#'
#' @description Render a star map for a given observer location, time, and atmospheric
#' conditions so it can be composited with [generate_sky()]. Returns a
#' `(resolution, 2 * resolution, 4)` array with an opaque alpha channel. An
#' image file is written only when `filename` is supplied.
#'
#' @param datetime Default `as.POSIXct("2000-01-01 00:00:00", tz = "UTC")`. `POSIXct` timestamp used to compute local sidereal time.
#' @param lon Default `0`. Observer longitude in degrees (east positive).
#' @param lat Default `0`. Observer latitude in degrees.
#' @param filename Default `NA`. Path to an image file to write. If `NA`, the
#'   image array is returned without writing.
#' @param resolution Default `2048`. Map half-width; the output image is `2 * resolution` × `resolution`.
#' @param turbidity Default `3.0`. Atmospheric turbidity controlling aerosol optical depth for extinction/reddening.
#' @param ozone_du Default `300.0`. Column ozone in Dobson Units used in atmospheric absorption.
#' @param altitude Default `0.0`. Observer altitude above mean sea level in meters.
#' @param color Default `TRUE`. If `TRUE`, render RGB star colors; if `FALSE`, render monochrome luminance.
#' @param star_width Default `1`. Approximate stellar point-spread size in pixels (controls apparent star sharpness).
#' @param upper_hemisphere_only Default `TRUE`. If `TRUE`, pixels below the local horizon are suppressed to match `generate_sky()`’s visible hemisphere.
#' @param atmosphere_effects Default `TRUE`. If `TRUE`, apply atmospheric extinction and color shift using `turbidity`, `ozone_du`, and `altitude`.
#' @param numbercores Default `1`. Number of CPU threads to use.
#'
#' @return Either the image array, or the array is invisibly returned if a file
#'   is written. The array has dimensions `(resolution, 2 * resolution, 4)`.
#' @note Writing to non-EXR formats will introduce precision loss because HDR
#'   data are quantised to the destination format, and low dynamic range outputs like PNG
#'   and JPEG files will not represent the true luminosity values encoded in the array.
#'
#' @export
#' @examples
#' # Note: exposure has been increased for all examples (via white_point) for
#' # ease of visibility in documentation
#'
#' # Basic star field over Washington, DC at a fixed time
#' if(run_documentation()) {
#' generate_stars(
#'   resolution = 400,
#'   lon        = -77.0369,
#'   lat        = 38.9072,
#'   datetime   = as.POSIXct("2025-03-21 02:20:00", tz = "EST"),
#'   color      = TRUE,
#'   star_width = 1,
#'   atmosphere_effects   = TRUE,
#'   upper_hemisphere_only = TRUE,
#'   numbercores = 2
#' ) |>
#'   rayimage::plot_image()
#'}
#'if(run_documentation()) {
#' # Monochrome stars, no atmospheric extinction/reddening, full sphere
#' generate_stars(
#'   resolution = 400,
#'   lon        = -122.4194,
#'   lat        = 37.7749,
#'   datetime   = as.POSIXct("2025-06-01 08:00:00", tz = "UTC"),
#'   color      = FALSE,
#'   star_width = 1,
#'   upper_hemisphere_only = FALSE,
#'   atmosphere_effects    = FALSE
#' ) |>
#'   rayimage::plot_image()
#'}
#'if(run_documentation()) {
#' # Sharper stars (smaller PSF) with ozone/turbidity and altitude
#' generate_stars(
#'   resolution = 400,
#'   lon        = 10,
#'   lat        = 45,
#'   datetime   = as.POSIXct("2025-12-01 22:00:00", tz = "UTC"),
#'   star_width = 0.5,
#'   turbidity  = 3.5,
#'   ozone_du   = 320,
#'   altitude   = 1000,
#'   color      = TRUE
#' ) |>
#'   rayimage::plot_image()
#'}
generate_stars = function(
	lon = 0,
	lat = 0,
	datetime = as.POSIXct("2000-01-01 00:00:00", tz = "UTC"),
	filename = NA,
	resolution = 2048,
	turbidity = 3.0,
	ozone_du = 300.0,
	altitude = 0.0,
	color = TRUE,
	star_width = 1,
	upper_hemisphere_only = TRUE,
	atmosphere_effects = TRUE,
	numbercores = 1
) {
	if (!inherits(datetime, "POSIXct")) {
		stop("datetime must be POSIXct in UTC")
	}
	attr(datetime, "tzone") = "UTC"
	jd = jd_utc(datetime)

	star_rgb = make_starfield_rcpp(
		stars = skymodelr::stars,
		resolution = resolution,
		lon_deg = lon,
		lat_deg = lat,
		jd = jd,
		use_rgb = color,
		turbidity = turbidity,
		ozone_du = ozone_du,
		altitude = altitude,
		star_width = star_width,
		atmosphere_effects = atmosphere_effects,
		upper_hemisphere_only = upper_hemisphere_only,
		numbercores = numbercores
	)
	star_array = array(0, dim = c(resolution, resolution * 2, 4))
	star_array[,, 1:3] = star_rgb
	star_array[,, 4] = 1
	star_array = rayimage::ray_read_image(
		star_array,
		assume_white = "D65",
		assume_colorspace = rayimage::CS_SRGB
	)
	if (!is.na(filename)) {
		warn_precision_loss(filename)
		rayimage::ray_write_image(star_array, filename)
		return(invisible(star_array))
	} else {
		return(star_array)
	}
}
