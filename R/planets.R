#' Generate a bright-planet image array
#'
#' @description Build a planetary luminance map aligned with the sky dome for
#' compositing within [generate_sky_latlong()].
#'
#' @param datetime POSIXct timestamp in UTC used for ephemerides.
#' @param lon Observer longitude in degrees (east positive).
#' @param lat Observer latitude in degrees.
#' @param filename Default `NA`. Destination image path to write. When `NA`, the
#'   image array is returned without writing.
#' @param resolution Default `2048`. Map half-width (image is `2 * resolution` × `resolution`).
#' @param turbidity Atmospheric turbidity for extinction modelling.
#' @param ozone_du Column ozone (Dobson Units) for colour shifts.
#' @param altitude Observer altitude in metres.
#' @param color Render RGB (`TRUE`) stars or monochrome (`FALSE`).
#' @param planet_width Approximate point-spread size for planets in pixels.
#' @param upper_hemisphere_only If `TRUE`, mask pixels below the horizon.
#' @param atmosphere_effects If `TRUE`, apply atmospheric extinction.
#' @param numbercores CPU threads used for rendering.
#' @param verbose Emit diagnostic output when `TRUE`.
#'
#' @return Either the image array, or the array is invisibly returned if a file
#'   is written. The array has dimensions `(resolution, 2 * resolution, 4)`.
#' @note Writing to non-EXR formats will introduce precision loss because HDR
#'   data are quantised to the destination format, and low dynamic range outputs like PNG
#'   and JPEG files will not represent the true luminosity values encoded in the array.
#'
#' @export
#' @examples
#' # Basic star field over Washington, DC at a fixed time
#' if(run_documentation()) {
#' generate_planets(
#'   datetime   = as.POSIXct("2025-03-21 02:20:00", tz = "EST"),
#'   lon        = -77.0369,
#'   lat        = 38.9072,
#'   resolution = 400,
#'   color      = TRUE,
#'   planet_width = 1,
#'   zero_point = 10,
#'   atmosphere_effects   = TRUE,
#'   upper_hemisphere_only = TRUE,
#'   numbercores = 2
#' ) |>
#'   rayimage::plot_image()
#'}
generate_planets = function(
	datetime = as.POSIXct("2000-01-01 00:00:00", tz = "UTC"),
	lon = 0,
	lat = 0,
	filename = NA,
	resolution = 2048,
	turbidity = 3.0,
	ozone_du = 300.0,
	altitude = 0.0,
	color = FALSE,
	planet_width = 1,
	upper_hemisphere_only = TRUE,
	atmosphere_effects = TRUE,
	numbercores = 1,
	verbose = FALSE
) {
	if (!inherits(datetime, "POSIXct")) {
		stop("datetime must be POSIXct in UTC")
	}
	attr(datetime, "tzone") = "UTC"
	jd = jd_utc(datetime)

	planet_temp = swe_dirs_topo_planets_df(datetime, lon, lat)
	if (verbose) {
		print(planet_temp)
	}
	planet_rgb = make_starfield_rcpp(
		stars = planet_temp,
		resolution = resolution,
		lon_deg = lon,
		lat_deg = lat,
		jd = jd,
		use_rgb = color,
		turbidity = turbidity,
		ozone_du = ozone_du,
		altitude = altitude,
		star_width = planet_width,
		atmosphere_effects = atmosphere_effects,
		upper_hemisphere_only = upper_hemisphere_only,
		numbercores = numbercores
	)
	planet_array = array(0, dim = c(resolution, resolution * 2, 4))
	planet_array[,, 1:3] = planet_rgb
	planet_array[,, 4] = 1
	planet_array = rayimage::ray_read_image(
		planet_array,
		assume_white = "D65",
		assume_colorspace = rayimage::CS_SRGB
	)
	if (!is.na(filename)) {
		warn_precision_loss(filename)
		rayimage::ray_write_image(planet_array, filename)
		return(invisible(planet_array))
	} else {
		return(planet_array)
	}
}
