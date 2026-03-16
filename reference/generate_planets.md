# Generate a bright-planet image array

Build a planetary luminance map aligned with the sky dome for
compositing within
[`generate_sky_latlong()`](https://tylermorganwall.github.io/skymodelr/reference/generate_sky_latlong.md).

## Usage

``` r
generate_planets(
  datetime,
  lon,
  lat,
  filename = NA,
  resolution = 2048,
  turbidity = 3,
  ozone_du = 300,
  altitude = 0,
  color = FALSE,
  planet_width = 1,
  upper_hemisphere_only = TRUE,
  atmosphere_effects = TRUE,
  number_cores = 1,
  verbose = FALSE
)
```

## Arguments

- datetime:

  POSIXct timestamp used for ephemerides.

- lon:

  Observer longitude in degrees (east positive).

- lat:

  Observer latitude in degrees.

- filename:

  Default `NA`. Destination image path to write. When `NA`, the image
  array is returned without writing.

- resolution:

  Default `2048`. Map half-width (image is `2 * resolution` ×
  `resolution`).

- turbidity:

  Atmospheric turbidity for extinction modelling.

- ozone_du:

  Column ozone (Dobson Units) for colour shifts.

- altitude:

  Observer altitude in metres.

- color:

  Render RGB (`TRUE`) stars or monochrome (`FALSE`).

- planet_width:

  Approximate point-spread size for planets in pixels.

- upper_hemisphere_only:

  If `TRUE`, mask pixels below the horizon.

- atmosphere_effects:

  If `TRUE`, apply atmospheric extinction.

- number_cores:

  CPU threads used for rendering.

- verbose:

  Emit diagnostic output when `TRUE`.

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
# Basic star field over Washington, DC at a fixed time
if(run_documentation()) {
generate_planets(
  datetime   = as.POSIXct("2025-03-21 02:20:00", tz = "EST"),
  lon        = -77.0369,
  lat        = 38.9072,
  resolution = 400,
  color      = TRUE,
  planet_width = 1,
  atmosphere_effects   = TRUE,
  upper_hemisphere_only = TRUE,
  number_cores = 2
) |>
  rayimage::plot_image()
}
```
