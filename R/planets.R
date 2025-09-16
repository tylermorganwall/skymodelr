#' Render bright planets into an EXR map
#'
#' @description Build a planetary luminance map aligned with the sky dome for
#' compositing within [generate_sky_latlong()].
#'
#' @param datetime POSIXct timestamp in UTC used for ephemerides.
#' @param lon Observer longitude in degrees (east positive).
#' @param lat Observer latitude in degrees.
#' @param filename Default `NA`. Destination `.exr` path to write.
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
	planet_tmp = tempfile(fileext = ".exr")
	if (!inherits(datetime, "POSIXct")) {
		stop("datetime must be POSIXct in UTC")
	}
	attr(datetime, "tzone") = "UTC"
	jd = jd_utc(datetime)

	planet_temp = swe_dirs_topo_planets_df(datetime, lon, lat)
	if (verbose) {
		print(planet_temp)
	}
	make_starfield_rcpp(
		outfile = planet_tmp,
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
		numbercores = numbercores,
		precision_multiplier = 1000000
	)
	planet_exr = rayimage::ray_read_image(planet_tmp) / 1000000
	if (!is.na(filename)) {
		file.copy(filename, planet_tmp, overwrite = TRUE)
		return(invisible(planet_exr))
	} else {
		return(planet_exr)
	}
}
