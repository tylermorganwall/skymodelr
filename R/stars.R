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

spec_to_bv = function(spec) {
	if (is.na(spec) || !nzchar(spec)) {
		return(NA_real_)
	}
	s = toupper(trimws(spec))
	letter = substr(s, 1, 1)
	base = c(
		O = -0.33,
		B = -0.30,
		A = 0.00,
		F = 0.30,
		G = 0.58,
		K = 0.81,
		M = 1.40
	)
	if (!letter %in% names(base)) {
		return(NA_real_)
	}
	digit = suppressWarnings(as.numeric(substr(s, 2, 2)))
	if (!is.finite(digit)) {
		digit = 0
	}
	digit = max(0, min(9, digit))
	next_letter = switch(
		letter,
		O = "B",
		B = "A",
		A = "F",
		F = "G",
		G = "K",
		K = "M",
		M = NA
	)
	next_base = if (!is.na(next_letter)) {
		base[[next_letter]]
	} else {
		base[["M"]] + 0.35
	}
	b_v = base[[letter]] + (digit / 10) * (next_base - base[[letter]])
	max(-0.4, min(2.0, b_v))
}

stars_with_bv = local({
	cache = NULL
	function() {
		if (!is.null(cache)) {
			return(cache)
		}
		star_df = skymodelr::stars
		if (!("b_v" %in% names(star_df))) {
			star_df$b_v = vapply(star_df$spec, spec_to_bv, numeric(1))
		}
		cache = star_df
		cache
	}
})

# ---------------------------------------------------------------------------
#' Generate a star‑field array aligned with `generate_sky()`
#'
#' @description Render a star map for a given observer location, time, and atmospheric
#' conditions so it can be composited with [generate_sky()]. Returns a
#' `(resolution, 2 * resolution, 4)` array with an opaque alpha channel. An
#' image file is written only when `filename` is supplied.
#'
#' @param datetime `POSIXct` timestamp used to compute local sidereal time.
#' @param lon Observer longitude in degrees (east positive).
#' @param lat Observer latitude in degrees.
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
#' @param number_cores Default `1`. Number of CPU threads to use.
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
#'   number_cores = 2
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
	lon,
	lat,
	datetime,
	filename = NA,
	resolution = 2048,
	turbidity = 3.0,
	ozone_du = 300.0,
	altitude = 0.0,
	color = TRUE,
	star_width = 1,
	upper_hemisphere_only = TRUE,
	atmosphere_effects = TRUE,
	number_cores = 1
) {
	if (!inherits(datetime, "POSIXct")) {
		stop("datetime must be POSIXct in UTC")
	}
	attr(datetime, "tzone") = "UTC"
	jd = jd_utc(datetime)

	star_rgb = make_starfield_rcpp(
		stars = stars_with_bv(),
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
		number_cores = number_cores
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
