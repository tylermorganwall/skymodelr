# Generate a star‑field array aligned with `generate_sky()`

Render a star map for a given observer location, time, and atmospheric
conditions so it can be composited with
[`generate_sky()`](https://tylermorganwall.github.io/skymodelr/reference/generate_sky.md).
Returns a `(resolution, 2 * resolution, 4)` array with an opaque alpha
channel. An image file is written only when `filename` is supplied.

## Usage

``` r
generate_stars(
  lon,
  lat,
  datetime,
  filename = NA,
  resolution = 2048,
  turbidity = 3,
  ozone_du = 300,
  altitude = 0,
  color = TRUE,
  star_width = 1,
  upper_hemisphere_only = TRUE,
  atmosphere_effects = TRUE,
  number_cores = 1
)
```

## Arguments

- lon:

  Observer longitude in degrees (east positive).

- lat:

  Observer latitude in degrees.

- datetime:

  `POSIXct` timestamp used to compute local sidereal time.

- filename:

  Default `NA`. Path to an image file to write. If `NA`, the image array
  is returned without writing.

- resolution:

  Default `2048`. Map half-width; the output image is `2 * resolution` ×
  `resolution`.

- turbidity:

  Default `3.0`. Atmospheric turbidity controlling aerosol optical depth
  for extinction/reddening.

- ozone_du:

  Default `300.0`. Column ozone in Dobson Units used in atmospheric
  absorption.

- altitude:

  Default `0.0`. Observer altitude above mean sea level in meters.

- color:

  Default `TRUE`. If `TRUE`, render RGB star colors; if `FALSE`, render
  monochrome luminance.

- star_width:

  Default `1`. Approximate stellar point-spread size in pixels (controls
  apparent star sharpness).

- upper_hemisphere_only:

  Default `TRUE`. If `TRUE`, pixels below the local horizon are
  suppressed to match
  [`generate_sky()`](https://tylermorganwall.github.io/skymodelr/reference/generate_sky.md)’s
  visible hemisphere.

- atmosphere_effects:

  Default `TRUE`. If `TRUE`, apply atmospheric extinction and color
  shift using `turbidity`, `ozone_du`, and `altitude`.

- number_cores:

  Default `1`. Number of CPU threads to use.

## Value

Either the image array, or the array is invisibly returned if a file is
written. The array has dimensions `(resolution, 2 * resolution, 4)`.

## Note

Writing to non-EXR formats will introduce precision loss because HDR
data are quantised to the destination format, and low dynamic range
outputs like PNG and JPEG files will not represent the true luminosity
values encoded in the array.

## Examples

``` r
# Note: exposure has been increased for all examples (via white_point) for
# ease of visibility in documentation

# Basic star field over Washington, DC at a fixed time
if(run_documentation()) {
generate_stars(
  resolution = 400,
  lon        = -77.0369,
  lat        = 38.9072,
  datetime   = as.POSIXct("2025-03-21 02:20:00", tz = "EST"),
  color      = TRUE,
  star_width = 1,
  atmosphere_effects   = TRUE,
  upper_hemisphere_only = TRUE,
  number_cores = 2
) |>
  rayimage::plot_image()
}

if(run_documentation()) {
# Monochrome stars, no atmospheric extinction/reddening, full sphere
generate_stars(
  resolution = 400,
  lon        = -122.4194,
  lat        = 37.7749,
  datetime   = as.POSIXct("2025-06-01 08:00:00", tz = "UTC"),
  color      = FALSE,
  star_width = 1,
  upper_hemisphere_only = FALSE,
  atmosphere_effects    = FALSE
) |>
  rayimage::plot_image()
}

if(run_documentation()) {
# Sharper stars (smaller PSF) with ozone/turbidity and altitude
generate_stars(
  resolution = 400,
  lon        = 10,
  lat        = 45,
  datetime   = as.POSIXct("2025-12-01 22:00:00", tz = "UTC"),
  star_width = 0.5,
  turbidity  = 3.5,
  ozone_du   = 320,
  altitude   = 1000,
  color      = TRUE
) |>
  rayimage::plot_image()
}
```
