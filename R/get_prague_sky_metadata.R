#' Prepare Location and Time Metadata for Native Prague Sky Queries
#'
#' Resolve the full-altitude Prague dataset and the Sun's topocentric position
#' without loading coefficients or generating an image. This lets rendering
#' libraries evaluate the model natively while using skymodelr's ephemeris and
#' RGB calibration through a public interface.
#'
#' @param datetime A single `POSIXct` instant. Its time zone is respected.
#' @param lat Observer latitude in degrees, between -90 and 90.
#' @param lon Observer longitude in degrees, between -180 and 180.
#' @param altitude Default `0`. Reference altitude in meters, from 0 to 15000.
#' @param visibility Default `50`. Meteorological visibility in kilometers,
#'   from 20 to 131.8.
#' @param albedo Default `0.5`. Ground albedo, from 0 to 1.
#' @param prague_rgb_correction Default `TRUE`. Apply the Prague RGB calibration,
#'   as in [calculate_sky_values()].
#' @param prague_rgb_correction_strength Default `1`. Nonnegative calibration
#'   strength.
#' @param prague_rgb_correction_gain Default `"auto"`. Automatic calibration or
#'   three positive RGB gains, as in [calculate_sky_values()].
#'
#' @return A list with `filename` (the installed full-altitude coefficient file),
#'   `elevation_deg`, `azimuth_deg` (clockwise from north),
#'   `angular_diameter_deg`, `rgb_gain` (the applied R, G, B gains), and the
#'   reference `altitude`, `visibility`, and `albedo`. No native pointers are
#'   included; the result can be serialized. The coefficient path is local to
#'   this computer and should be resolved again on another machine.
#'
#' @details Install the data explicitly with
#'   `download_sky_data(sea_level = FALSE)`. This function never downloads files
#'   or prompts. It supports the visible-spectrum, full-altitude dataset and
#'   checks the Prague solar-elevation range of -4.2 to 90 degrees. It does not
#'   evaluate radiance, attenuation, or atmospheric refraction.
#' @export
#' @examplesIf interactive()
#' info = get_prague_sky_metadata(
#'   as.POSIXct("2026-06-21 20:00:00", tz = "America/New_York"),
#'   lat = 40.7, lon = -74
#' )
#' info[c("elevation_deg", "azimuth_deg", "rgb_gain")]
get_prague_sky_metadata <- function(
  datetime,
  lat,
  lon,
  altitude = 0,
  visibility = 50,
  albedo = 0.5,
  prague_rgb_correction = TRUE,
  prague_rgb_correction_strength = 1,
  prague_rgb_correction_gain = "auto"
) {
  celestial_disk_settings(
    datetime,
    lat,
    lon,
    altitude,
    16,
    visibility,
    albedo,
    1,
    FALSE,
    prague_rgb_correction,
    prague_rgb_correction_strength,
    prague_rgb_correction_gain
  )
  filename <- resolve_prague_coef_file(altitude = 1, allow_download = FALSE)
  ephemeris <- swe_dirs_topo_moon_sun(datetime, lat, lon, elev_m = altitude)
  disk <- celestial_disk_result(
    NULL,
    ephemeris$sun_dir_topo,
    ephemeris$sun_diameter_degrees
  )
  if (disk$elevation_deg < -4.2 || disk$elevation_deg > 90) {
    stop(
      "The Prague model supports Sun elevations from -4.2 to 90 degrees.",
      call. = FALSE
    )
  }
  gain <- if (
    normalize_prague_rgb_correction(prague_rgb_correction) == "constant"
  ) {
    prepare_prague_rgb_gain(
      prague_rgb_correction_gain,
      prague_rgb_correction_strength
    )
  } else {
    c(R = 1, G = 1, B = 1)
  }
  c(
    list(filename = filename),
    disk[c("elevation_deg", "azimuth_deg", "angular_diameter_deg")],
    list(
      rgb_gain = gain,
      altitude = altitude,
      visibility = visibility,
      albedo = albedo
    )
  )
}
