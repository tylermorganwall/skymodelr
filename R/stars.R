# ---------------------------------------------------------------------------
# Helpers to convert POSIXct → Julian Date (UTC noon 1 Jan 2000 = 2451545).
#' @noRd
jd_utc = function(time_utc) {
  unclass(time_utc) / 86400 + 2440587.5 # Unix epoch to JD
}

# ---------------------------------------------------------------------------
#' Write a star‑field EXR aligned with `generate_sky()`
#'
#' @param outfile Default `NA`. Path to a `.exr` file to write. If `NA`, no file is written and the image data are returned.
#' @param resolution Default `2048`. Map half-width; the output image is `2 * resolution` × `resolution`.
#' @param lon Default `0`. Observer longitude in degrees (east positive).
#' @param lat Default `0`. Observer latitude in degrees.
#' @param datetime Default `as.POSIXct("2000-01-01 00:00:00", tz = "UTC")`. `POSIXct` timestamp used to compute local sidereal time.
#' @param turbidity Default `3.0`. Atmospheric turbidity controlling aerosol optical depth for extinction/reddening.
#' @param ozone_du Default `300.0`. Column ozone in Dobson Units used in atmospheric absorption.
#' @param altitude Default `0.0`. Observer altitude above mean sea level in meters.
#' @param color Default `TRUE`. If `TRUE`, render RGB star colors; if `FALSE`, render monochrome luminance.
#' @param star_width Default `1`. Approximate stellar point-spread size in pixels (controls apparent star sharpness).
#' @param upper_hemisphere_only Default `TRUE`. If `TRUE`, pixels below the local horizon are suppressed to match `generate_sky()`’s visible hemisphere.
#' @param atmosphere_effects Default `TRUE`. If `TRUE`, apply atmospheric extinction and color shift using `turbidity`, `ozone_du`, and `altitude`.
#' @param numbercores Default `1`. Number of CPU threads to use.
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
#'   zero_point = 10,
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
#'   zero_point = 10,
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
  outfile = NA,
  resolution = 2048,
  zero_point = 1,
  lon = 0,
  lat = 0,
  datetime = as.POSIXct("2000-01-01 00:00:00", tz = "UTC"),
  turbidity = 3.0,
  ozone_du = 300.0,
  altitude = 0.0,
  color = TRUE,
  star_width = 1,
  upper_hemisphere_only = TRUE,
  atmosphere_effects = TRUE,
  numbercores = 1
) {
  starfield_tmp = tempfile(fileext = ".exr")
  if (!inherits(datetime, "POSIXct")) {
    stop("datetime must be POSIXct in UTC")
  }
  attr(datetime, "tzone") = "UTC"
  jd = jd_utc(datetime)

  make_starfield_rcpp(
    outfile = starfield_tmp,
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
  star_exr = rayimage::ray_read_image(starfield_tmp)
  if (!is.na(outfile)) {
    file.copy(outfile, starfield_tmp, overwrite = TRUE)
    return(invisible(star_exr))
  } else {
    return(star_exr)
  }
}
