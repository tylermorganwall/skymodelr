# Generate a Hosek-Wilkie sky dome array

Evaluate either the Hosek-Wilkie or Prague analytic sky models and
return a high-dynamic-range image array for the given solar
configuration. An image file is written only when `filename` is
supplied.

## Usage

``` r
generate_sky(
  filename = NA,
  albedo = 0.1,
  turbidity = 3,
  elevation = 10,
  azimuth = 90,
  altitude = 0,
  resolution = 2048,
  number_cores = 1,
  hosek = TRUE,
  wide_spectrum = FALSE,
  visibility = 50,
  verbose = FALSE,
  render_mode = "all",
  below_horizon = TRUE
)
```

## Arguments

- filename:

  Default `NA`. Path to an image file to write. If not given, the array
  is returned instead.

- albedo:

  Default `0.1`. 0.0-1.0 ground albedo. Grass has an albedo of about
  0.09, while a landscape covered in snow will have an albedo of 1.0.

- turbidity:

  Default `3`. 1.7-10 atmospheric turbidity. Only valid for Hosek model.

- elevation:

  Default `10`. Solar elevation above the horizon (degrees).

- azimuth:

  Default `90`, sun directly east. Solar azimuth (degrees). The left
  edge of the image faces north and the middle faces south.

- altitude:

  Default `0`. Altitude of the viewer in meters. Valid range: 0
  to 15000. Only valid for the Prague model.

- resolution:

  Default `2048`. Height of the image. Width is twice this number.

- number_cores:

  Default `1`. Number of threads to use in computation.

- hosek:

  Default `TRUE`. Set to `FALSE` to enable the Prague 2021-22 spectral
  sky model.

- wide_spectrum:

  Default `FALSE`. Use the 55-channel Prague coefficients (sea level
  only).

- visibility:

  Default `50`. Meteorological range in kilometres for Prague model.

- verbose:

  Default `FALSE`. Whether to print progress bars/diagnostic info.

- render_mode:

  Default `"all"`. One of `"all"`, `"atmosphere"`, or `"sun"`. Use
  `"all"` for atmosphere + solar disk, `"atmosphere"` for atmospheric
  radiance only, or `"sun"` for the solar disk only.

- below_horizon:

  Default `TRUE`. Whether to sample atmospheric scattering below the
  horizon, which is non-zero when altitude \> 0.

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
if(run_documentation()) {
# Hosek model (default): clear morning, Sun SE, with solar disk
generate_sky(
  resolution = 400,
  elevation  = 15,
  azimuth    = 135,
  turbidity  = 3,
  render_mode = "all"
) |>
  rayimage::plot_image()
}

if(run_documentation()) {
# Same view but hazier and without the solar disk
generate_sky(
  resolution = 400,
  elevation  = 15,
  azimuth    = 135,
  turbidity  = 6,
  render_mode = "atmosphere"
) |>
  rayimage::plot_image()
}

# Prague model (may prompt to download coefficients on first use)
# \donttest{
if(run_documentation()) {
generate_sky(
  resolution = 400,
  hosek      = FALSE,
  altitude   = 0,
  visibility = 80,
  albedo     = 0.2,
  elevation  = 5,
  azimuth    = 220,
  number_cores = 2
) |>
  rayimage::plot_image()
}
#>  Coefficient file for this setting not yet present: this is a large file (107MB), download? [y/n] 
#> Error in check_coef_file("SkyModelDatasetGround.dat"): Input not recognized.
# }
```
