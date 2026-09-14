# Prepare Location and Time Metadata for Native Prague Sky Queries

Resolve the full-altitude Prague dataset and the Sun's topocentric
position without loading coefficients or generating an image. This lets
rendering libraries evaluate the model natively while using skymodelr's
ephemeris and RGB calibration through a public interface.

## Usage

``` r
get_prague_sky_metadata(
  datetime,
  lat,
  lon,
  altitude = 0,
  visibility = 50,
  albedo = 0.5,
  prague_rgb_correction = TRUE,
  prague_rgb_correction_strength = 1,
  prague_rgb_correction_gain = "auto"
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

  Default `0`. Reference altitude in meters, from 0 to 15000.

- visibility:

  Default `50`. Meteorological visibility in kilometers, from 20 to
  131.8.

- albedo:

  Default `0.5`. Ground albedo, from 0 to 1.

- prague_rgb_correction:

  Default `TRUE`. Apply the Prague RGB calibration, as in
  [`calculate_sky_values()`](https://tylermorganwall.github.io/skymodelr/reference/calculate_sky_values.md).

- prague_rgb_correction_strength:

  Default `1`. Nonnegative calibration strength.

- prague_rgb_correction_gain:

  Default `"auto"`. Automatic calibration or three positive RGB gains,
  as in
  [`calculate_sky_values()`](https://tylermorganwall.github.io/skymodelr/reference/calculate_sky_values.md).

## Value

A list with `filename` (the installed full-altitude coefficient file),
`elevation_deg`, `azimuth_deg` (clockwise from north),
`angular_diameter_deg`, `rgb_gain` (the applied R, G, B gains), and the
reference `altitude`, `visibility`, and `albedo`. No native pointers are
included; the result can be serialized. The coefficient path is local to
this computer and should be resolved again on another machine.

## Details

Install the data explicitly with `download_sky_data(sea_level = FALSE)`.
This function never downloads files or prompts. It supports the
visible-spectrum, full-altitude dataset and returns the actual Sun
position even below the Prague model's minimum elevation of -4.2
degrees. Native consumers should use zero solar radiance below that
limit while retaining transmission for other light sources. This
function does not evaluate radiance, attenuation, or atmospheric
refraction.

## Examples

``` r
if (FALSE) { # interactive()
info = get_prague_sky_metadata(
  as.POSIXct("2026-06-21 20:00:00", tz = "America/New_York"),
  lat = 40.7, lon = -74
)
info[c("elevation_deg", "azimuth_deg", "rgb_gain")]
}
```
