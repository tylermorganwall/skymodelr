# Generate Sun and Moon Disk Radiance Images

Generate a celestial disk directly at texture resolution, independently
of an environment map. The Sun uses the Prague solar radiance profile.
The Moon uses the existing lunar surface texture, topocentric phase
orientation, earthshine, photometry, and atmospheric extinction.

## Usage

``` r
generate_sun_disk(
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
  prague_rgb_correction_gain = "auto",
  atmospheric_attenuation = TRUE
)

generate_moon_disk(
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
  moon_extinction_kV = 0.172,
  atmospheric_attenuation = TRUE
)
```

## Arguments

- datetime:

  A single `POSIXct` instant. Its time zone is respected.

- lat:

  Observer latitude in degrees, between -90 and 90.

- lon:

  Observer longitude in degrees, between -180 and 180.

- altitude:

  Default `0`. Observer altitude in meters, from 0 to 15000.

- resolution:

  Default `256` for the Sun and `1024` for the Moon. Target disk width
  and height in pixels, an integer of at least 16. The Moon is cropped
  from a padded raster without resampling; edge coverage can add a few
  pixels to the returned dimensions.

- visibility:

  Default `50`. Prague visibility in kilometers, from 20 to 131.8.

- albedo:

  Default `0.5`. Ground albedo, from 0 to 1.

- number_cores:

  Default `1`. Positive integer number of worker threads.

- wide_spectrum:

  Default `FALSE`. Use the extended Prague spectral data.

- prague_rgb_correction:

  Default `TRUE`. RGB correction passed to
  [`calculate_sky_values()`](https://tylermorganwall.github.io/skymodelr/reference/calculate_sky_values.md).

- prague_rgb_correction_strength:

  Default `1`. Nonnegative correction strength passed to
  [`calculate_sky_values()`](https://tylermorganwall.github.io/skymodelr/reference/calculate_sky_values.md).

- prague_rgb_correction_gain:

  Default `"auto"`. Automatic correction or three positive RGB gains, as
  in
  [`calculate_sky_values()`](https://tylermorganwall.github.io/skymodelr/reference/calculate_sky_values.md).

- atmospheric_attenuation:

  Default `TRUE`. Include atmospheric extinction and tint for the
  specified observer. Set `FALSE` to return radiance before it enters
  Earth's atmosphere; the consuming renderer must apply atmospheric
  transport and horizon visibility. Lunar phase, surface detail, and
  earthshine are retained. `moon_extinction_kV` is ignored when this is
  `FALSE`.

- earthshine:

  Default `TRUE`. Include Earth-reflected light on the Moon.

- earthshine_albedo:

  Default `0.19`. Nonnegative Earth albedo used for earthshine.

- solar_irradiance_w_m2:

  Default `1300`. Positive solar irradiance in watts per square meter
  used to illuminate the lunar surface.

- moon_extinction_kV:

  Default `0.172`. Nonnegative lunar atmospheric extinction coefficient
  in magnitudes per airmass.

## Value

A list containing:

- `image`: a numeric array with dimensions `c(height, width, 3)`
  containing linear sRGB radiance, with no exposure, tone mapping, or
  alpha channel. Antialiased lunar coverage is already applied: do not
  multiply by alpha. Out-of-gamut solar RGB values may be negative and
  should be preserved.

- `azimuth_deg`: disk-center azimuth, clockwise from true north, in \[0,
  360).

- `elevation_deg`: disk-center elevation above the local horizon, in
  degrees.

- `angular_diameter_deg`: topocentric apparent angular diameter in
  degrees.

- `atmospheric_attenuation`: whether atmospheric filtering is included.

- `projection`: `"rectilinear"`, using the disk mapping described below.

## Details

Prague data must be installed with
[`download_sky_data()`](https://tylermorganwall.github.io/skymodelr/reference/download_sky_data.md).
The atmosphere is evaluated for one observer altitude and time when
`atmospheric_attenuation` is `TRUE`. With `FALSE`, solar radiance comes
from the model's intrinsic solar spectrum and lunar radiance uses
unattenuated photometry without the atmospheric tint. Solar texture
values are then independent of altitude and solar elevation; ephemeris
placement and lunar phase still use the requested location and time.
Neither mode clips the texture to a geometric horizon. Pair these images
with an atmosphere-only sky and disable the corresponding rasterized
celestial bodies in that sky to avoid counting their light twice.

The disk is inscribed in the image rectangle. For column `j`, row `i`,
width `w`, and height `h`, pixel centers have coordinates
`x = 2 * (j - 0.5) / w - 1` and `y = 1 - 2 * (i - 0.5) / h`. Sample only
within `x^2 + y^2 <= 1`. The corresponding ray direction is
`normalize(forward + tan(radius) * (x * right + y * up))`, where
`radius = angular_diameter_deg * pi / 360`. Image-up is the projection
of local vertical into the disk plane; image-right points toward
increasing azimuth. At the zenith or nadir, north replaces the vertical
reference axis. There is no image reprojection or implicit rotation in
the return value.

The solar model has a fixed angular profile, mapped here to the apparent
ephemeris diameter without changing radiance. Lunar radiance is
normalized over disk solid angle to the existing phase-dependent
irradiance, with the solar spectral RGB distribution. With
`atmospheric_attenuation = TRUE`, lunar attenuation and Prague tint are
evaluated at the disk center. At or below zero center elevation, their
finite horizon values are retained so the upper limb does not disappear
prematurely. This continuation preserves a full disk texture; it is not
a model of below-horizon visibility or horizon depression at altitude.
Geometric horizon clipping is the consumer's responsibility; sample the
full disk and discard directions below the horizon when appropriate.

To save an EXR without changing radiance, tag `image` with
`rayimage::ray_read_image(image, normalize = FALSE, source_linear = TRUE, assume_colorspace = rayimage::CS_SRGB)`
and write with `clamp = FALSE`.

## Examples

``` r
if (FALSE) { # interactive()
time = as.POSIXct("2026-01-28 21:00:00", tz = "Pacific/Auckland")
moon = generate_moon_disk(time, lat = -36.87593, lon = 174.7647)
moon[c("azimuth_deg", "elevation_deg", "angular_diameter_deg")]
# Let a renderer evaluate Earth's atmosphere at each interaction instead.
intrinsic_moon = generate_moon_disk(time, lat = -36.87593, lon = 174.7647,
  atmospheric_attenuation = FALSE)
}
```
