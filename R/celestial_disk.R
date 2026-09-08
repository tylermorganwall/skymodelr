#' Generate Sun and Moon Disk Radiance Images
#'
#' Generate a celestial disk directly at texture resolution, independently of
#' an environment map. The Sun uses the Prague solar radiance profile. The Moon
#' uses the existing lunar surface texture, topocentric phase orientation,
#' earthshine, photometry, and atmospheric extinction.
#'
#' @param datetime A single `POSIXct` instant. Its time zone is respected.
#' @param lat Observer latitude in degrees, between -90 and 90.
#' @param lon Observer longitude in degrees, between -180 and 180.
#' @param altitude Default `0`. Observer altitude in meters, from 0 to 15000.
#' @param resolution Default `256` for the Sun and `1024` for the Moon. Target
#'   disk width and height in pixels, an integer of at least 16. The Moon is
#'   cropped from a padded raster without resampling; edge coverage can add a
#'   few pixels to the returned dimensions.
#' @param visibility Default `50`. Prague visibility in kilometers, from 20 to
#'   131.8.
#' @param albedo Default `0.5`. Ground albedo, from 0 to 1.
#' @param number_cores Default `1`. Positive integer number of worker threads.
#' @param wide_spectrum Default `FALSE`. Use the extended Prague spectral data.
#' @param prague_rgb_correction Default `TRUE`. RGB correction passed to
#'   [calculate_sky_values()].
#' @param prague_rgb_correction_strength Default `1`. Nonnegative correction
#'   strength passed to [calculate_sky_values()].
#' @param prague_rgb_correction_gain Default `"auto"`. Automatic correction or
#'   three positive RGB gains, as in [calculate_sky_values()].
#' @param earthshine Default `TRUE`. Include Earth-reflected light on the Moon.
#' @param earthshine_albedo Default `0.19`. Nonnegative Earth albedo used for
#'   earthshine.
#' @param solar_irradiance_w_m2 Default `1300`. Positive solar irradiance in
#'   watts per square meter used to illuminate the lunar surface.
#' @param moon_extinction_kV Default `0.172`. Nonnegative lunar atmospheric
#'   extinction coefficient in magnitudes per airmass.
#'
#' @return A list containing:
#' * `image`: a numeric array with dimensions `c(height, width, 3)` containing
#'   linear sRGB radiance, with no exposure, tone mapping, or alpha channel.
#'   Antialiased lunar coverage is already applied: do not multiply by alpha.
#'   Out-of-gamut solar RGB values may be negative and should be preserved.
#' * `azimuth_deg`: disk-center azimuth, clockwise from true north, in [0, 360).
#' * `elevation_deg`: disk-center elevation above the local horizon, in degrees.
#' * `angular_diameter_deg`: topocentric apparent angular diameter in degrees.
#' * `projection`: `"rectilinear"`, using the disk mapping described below.
#'
#' @details
#' Prague data must be installed with [download_sky_data()]. The atmosphere is
#' evaluated for one observer altitude and time; this is a static light image.
#' Pair these images with an atmosphere-only sky and disable the corresponding
#' rasterized celestial bodies in that sky to avoid counting their light twice.
#'
#' The disk is inscribed in the image rectangle. For column `j`, row `i`, width
#' `w`, and height `h`, pixel centers have coordinates
#' `x = 2 * (j - 0.5) / w - 1` and `y = 1 - 2 * (i - 0.5) / h`.
#' Sample only within `x^2 + y^2 <= 1`. The corresponding ray direction is
#' `normalize(forward + tan(radius) * (x * right + y * up))`, where
#' `radius = angular_diameter_deg * pi / 360`. Image-up is the projection of
#' local vertical into the disk plane; image-right points toward increasing
#' azimuth. At the zenith or nadir, north replaces the vertical reference axis.
#' There is no image reprojection or implicit rotation in the return value.
#'
#' The solar model has a fixed angular profile, mapped here to the apparent
#' ephemeris diameter without changing radiance. Lunar radiance is normalized
#' over disk solid angle to the existing phase-dependent, atmosphere-attenuated
#' irradiance, with the solar spectral RGB distribution and Prague atmospheric
#' tint. Lunar attenuation and tint are evaluated at the disk center. At or
#' below zero center elevation, their finite horizon values are retained so the
#' upper limb does not disappear prematurely. This continuation preserves a full
#' disk texture; it is not a model of below-horizon visibility or horizon
#' depression at altitude. Geometric horizon clipping is the consumer's
#' responsibility; sample the full disk and discard directions below the 
#' horizon when appropriate.
#'
#' To save an EXR without changing radiance, tag `image` with
#' `rayimage::ray_read_image(image, normalize = FALSE, source_linear = TRUE,
#' assume_colorspace = rayimage::CS_SRGB)` and write with `clamp = FALSE`.
#'
#' @export
#' @examplesIf interactive()
#' time = as.POSIXct("2026-01-28 21:00:00", tz = "Pacific/Auckland")
#' moon = generate_moon_disk(time, lat = -36.87593, lon = 174.7647)
#' moon[c("azimuth_deg", "elevation_deg", "angular_diameter_deg")]
generate_sun_disk = function(
  datetime,
  lat,
  lon,
  altitude = 0,
  resolution = 256,
  visibility = 50,
  albedo = 0.5,
  number_cores = 1,
  wide_spectrum = FALSE,
  prague_rgb_correction = TRUE,
  prague_rgb_correction_strength = 1,
  prague_rgb_correction_gain = "auto"
) {
  settings = celestial_disk_settings(
    datetime,
    lat,
    lon,
    altitude,
    resolution,
    visibility,
    albedo,
    number_cores,
    wide_spectrum,
    prague_rgb_correction,
    prague_rgb_correction_strength,
    prague_rgb_correction_gain
  )
  ephemeris = swe_dirs_topo_moon_sun(datetime, lat, lon, elev_m = altitude)
  d = ephemeris$sun_dir_topo
  image = celestial_sun_pixels(c(d[1], d[3], -d[2]), resolution, settings)
  celestial_disk_result(image, d, ephemeris$sun_diameter_degrees)
}

#' @rdname generate_sun_disk
#' @export
generate_moon_disk = function(
  datetime,
  lat,
  lon,
  altitude = 0,
  resolution = 1024,
  visibility = 50,
  albedo = 0.5,
  number_cores = 1,
  wide_spectrum = FALSE,
  prague_rgb_correction = TRUE,
  prague_rgb_correction_strength = 1,
  prague_rgb_correction_gain = "auto",
  earthshine = TRUE,
  earthshine_albedo = 0.19,
  solar_irradiance_w_m2 = 1300,
  moon_extinction_kV = 0.172
) {
  settings = celestial_disk_settings(
    datetime,
    lat,
    lon,
    altitude,
    resolution,
    visibility,
    albedo,
    number_cores,
    wide_spectrum,
    prague_rgb_correction,
    prague_rgb_correction_strength,
    prague_rgb_correction_gain
  )
  validate_disk_logical(earthshine, "earthshine")
  validate_disk_scalar(earthshine_albedo, "earthshine_albedo", 0, Inf)
  validate_disk_scalar(solar_irradiance_w_m2, "solar_irradiance_w_m2", 0, Inf)
  if (solar_irradiance_w_m2 == 0) {
    stop("solar_irradiance_w_m2 must be positive.", call. = FALSE)
  }
  validate_disk_scalar(moon_extinction_kV, "moon_extinction_kV", 0, Inf)
  args = list(
    earthshine = earthshine,
    earthshine_albedo = earthshine_albedo,
    solar_irradiance_w_m2 = solar_irradiance_w_m2,
    moon_extinction_kV = moon_extinction_kV
  )
  ephemeris = swe_dirs_topo_moon_sun(
    datetime,
    lat,
    lon,
    elev_m = altitude,
    moon_extinction_kV = moon_extinction_kV
  )
  image = celestial_moon_pixels(
    datetime,
    lat,
    lon,
    resolution,
    args,
    ephemeris,
    settings
  )
  celestial_disk_result(
    image,
    ephemeris$moon_dir_topo,
    ephemeris$moon_diameter_degrees
  )
}

#' @keywords internal
validate_disk_scalar = function(value, name, lower, upper, integer = FALSE) {
  if (
    !is.numeric(value) ||
      length(value) != 1L ||
      !is.finite(value) ||
      value < lower ||
      value > upper ||
      (integer && value != floor(value))
  ) {
    stop(
      name,
      " must be a finite ",
      if (integer) "integer" else "number",
      " between ",
      lower,
      " and ",
      upper,
      ".",
      call. = FALSE
    )
  }
}

#' @keywords internal
validate_disk_logical = function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
}

#' @keywords internal
celestial_disk_settings = function(
  datetime,
  lat,
  lon,
  altitude,
  resolution,
  visibility,
  albedo,
  number_cores,
  wide_spectrum,
  prague_rgb_correction,
  prague_rgb_correction_strength,
  prague_rgb_correction_gain
) {
  if (
    !inherits(datetime, "POSIXct") ||
      length(datetime) != 1L ||
      !is.finite(as.numeric(datetime))
  ) {
    stop("datetime must be a single finite POSIXct instant.", call. = FALSE)
  }
  validate_disk_scalar(lat, "lat", -90, 90)
  validate_disk_scalar(lon, "lon", -180, 180)
  validate_disk_scalar(altitude, "altitude", 0, 15000)
  validate_disk_scalar(
    resolution,
    "resolution",
    16,
    .Machine$integer.max / 2,
    TRUE
  )
  validate_disk_scalar(visibility, "visibility", 20, 131.8)
  validate_disk_scalar(albedo, "albedo", 0, 1)
  validate_disk_scalar(
    number_cores,
    "number_cores",
    1,
    .Machine$integer.max,
    TRUE
  )
  validate_disk_logical(wide_spectrum, "wide_spectrum")
  normalize_prague_rgb_correction(prague_rgb_correction)
  prepare_prague_rgb_gain(
    prague_rgb_correction_gain,
    prague_rgb_correction_strength
  )
  list(
    altitude = altitude,
    visibility = visibility,
    albedo = albedo,
    number_cores = number_cores,
    wide_spectrum = wide_spectrum,
    prague_rgb_correction = prague_rgb_correction,
    prague_rgb_correction_strength = prague_rgb_correction_strength,
    prague_rgb_correction_gain = prague_rgb_correction_gain
  )
}

#' @keywords internal
celestial_disk_result = function(image, direction, diameter) {
  if (any(!is.finite(image))) {
    stop(
      "Nonfinite celestial radiance returned by the sky model.",
      call. = FALSE
    )
  }
  # swe_dirs_topo_moon_sun retains Swiss Ephemeris's south-to-west azimuth
  # convention. Convert to the public north-to-east convention used by the sky.
  list(
    image = image,
    azimuth_deg = (180 + atan2(direction[1], direction[2]) * 180 / pi) %% 360,
    elevation_deg = asin(pmax(-1, pmin(1, direction[3]))) * 180 / pi,
    angular_diameter_deg = diameter,
    projection = "rectilinear"
  )
}

#' @keywords internal
celestial_frame = function(direction) {
  forward = direction / sqrt(sum(direction^2))
  vertical = if (abs(forward[2]) < 0.999999) c(0, 1, 0) else c(0, 0, 1)
  right = c(
    forward[2] * vertical[3] - forward[3] * vertical[2],
    forward[3] * vertical[1] - forward[1] * vertical[3],
    forward[1] * vertical[2] - forward[2] * vertical[1]
  )
  right = right / sqrt(sum(right^2))
  up = c(
    right[2] * forward[3] - right[3] * forward[2],
    right[3] * forward[1] - right[1] * forward[3],
    right[1] * forward[2] - right[2] * forward[1]
  )
  list(forward = forward, right = right, up = up)
}

#' @keywords internal
celestial_sun_pixels = function(direction, n, settings) {
  elevation = asin(direction[2]) * 180 / pi
  if (elevation < -4.2) {
    return(array(0, c(n, n, 3)))
  }
  azimuth = (atan2(-direction[1], direction[3]) * 180 / pi) %% 360
  frame = celestial_frame(direction)
  grid = expand.grid(
    y = 1 - 2 * ((seq_len(n) - 0.5) / n),
    x = 2 * ((seq_len(n) - 0.5) / n) - 1
  )
  # SUN_RADIUS in src/PragueSkyModel/PragueSkyModel.cpp. Sample its profile
  # here; consumers map the resulting image to the ephemeris angular diameter.
  radius = 0.004654793
  dirs = outer(grid$x * tan(radius), frame$right) +
    outer(grid$y * tan(radius), frame$up) +
    matrix(frame$forward, n * n, 3, byrow = TRUE)
  dirs = dirs / sqrt(rowSums(dirs^2))
  phi = (atan2(-dirs[, 1], dirs[, 3]) * 180 / pi) %% 360
  theta = asin(pmax(-1, pmin(1, dirs[, 2]))) * 180 / pi
  values = do.call(
    calculate_sky_values,
    c(
      list(
        phi = phi,
        theta = theta,
        elevation = elevation,
        azimuth = azimuth,
        render_mode = "sun"
      ),
      settings
    )
  )
  array(as.numeric(values), c(n, n, 3))
}

#' @keywords internal
celestial_moon_pixels = function(
  datetime,
  lat,
  lon,
  resolution,
  args,
  ephemeris,
  settings
) {
  # skymodelr rasterizes a padded patch using rayvertex. Scope its thread count;
  # preparation happens before the path tracer's worker pool starts.
  old_options = options(cores = settings$number_cores)
  on.exit(options(old_options), add = TRUE)
  patch = do.call(
    generate_moon_image_latlong,
    c(
      list(
        datetime = datetime,
        lat = lat,
        lon = lon,
        elev_m = settings$altitude,
        width = 2 * resolution,
        height = 2 * resolution
      ),
      args
    )
  )$moon_luminance_array
  mask = patch[,, 4] > 0
  indices = which(mask, arr.ind = TRUE)
  if (!nrow(indices)) {
    stop("skymodelr returned an empty Moon texture.", call. = FALSE)
  }
  rows = seq(min(indices[, 1]), max(indices[, 1]))
  cols = seq(min(indices[, 2]), max(indices[, 2]))
  patch = patch[rows, cols, , drop = FALSE]
  pixels = patch[,, 1:3, drop = FALSE]
  # The rasterizer's reconstruction filter can undershoot at the silhouette.
  # Remove negative coverage/radiance before the physical normalization.
  coverage = pmax(0, pmin(1, patch[,, 4]))
  for (c in 1:3) {
    pixels[,, c] = pmax(0, pixels[,, c]) * coverage
  }

  direction = c(
    ephemeris$moon_dir_topo[1],
    ephemeris$moon_dir_topo[3],
    -ephemeris$moon_dir_topo[2]
  )
  elevation = asin(pmax(-1, pmin(1, direction[2]))) * 180 / pi
  # Calibrate a complete disk. Its center may be below zero while the upper
  # limb remains visible; horizon visibility belongs to the consuming renderer.
  illuminance = apply_airmass_extinction(
    ephemeris$moon_brightness_lux_unattenuated,
    elevation,
    kV = args$moon_extinction_kV,
    clip_horizon = FALSE
  )
  irradiance = lux_to_radiometric_irradiance(
    illuminance,
    compute_K_eff("BB5778")
  )
  rgb = compute_spd_rgb_unit("BB5778")
  # Use the same finite horizon limit for color, avoiding a tint discontinuity
  # when the center crosses zero. This remains a disk-wide approximation.
  tint_elevation = max(0, elevation)
  azimuth = (atan2(-direction[1], direction[3]) * 180 / pi) %% 360
  tint = do.call(
    calculate_sky_values,
    c(
      list(
        phi = azimuth,
        theta = tint_elevation,
        elevation = tint_elevation,
        azimuth = azimuth,
        render_mode = "sun"
      ),
      settings
    )
  )
  if (all(is.finite(tint)) && sum(tint) > 0) {
    rgb = rgb * pmax(as.numeric(tint), 0)
  }
  rgb = rgb / sum(rgb)
  n = dim(pixels)
  xx = matrix(
    rep(2 * ((seq_len(n[2]) - 0.5) / n[2]) - 1, each = n[1]),
    n[1],
    n[2]
  )
  yy = matrix(rep(2 * ((seq_len(n[1]) - 0.5) / n[1]) - 1, n[2]), n[1], n[2])
  r2 = xx^2 + yy^2
  t2 = tan(ephemeris$moon_diameter_degrees * pi / 360)^2
  domega = 4 * t2 / (1 + t2 * r2)^1.5 / (n[1] * n[2])
  domega[r2 > 1] = 0
  for (c in 1:3) {
    integral = sum(pixels[,, c] * domega)
    pixels[,, c] = if (integral > 0) {
      pixels[,, c] * irradiance * rgb[c] / integral
    } else {
      0
    }
  }
  pixels
}
