# Sample a direction from the Prague model.

Evaluate the Prague spectral sky model at arbitrary spherical directions
without writing an image, returning radiance-only samples.

## Usage

``` r
calculate_sky_values(
  phi,
  theta,
  altitude = 0,
  elevation = 10,
  visibility = 50,
  albedo = 0.5,
  azimuth = 90,
  number_cores = 1,
  wide_spectrum = FALSE,
  render_mode = "all"
)
```

## Arguments

- phi:

  Horizontal angle of the sample, degrees. Vectorized. Range 0 to 360.

- theta:

  Vertical angle of the sample, degrees. Vectorized. Range -90 to 90.

- altitude:

  Default `0`, vectorized. Altitude of the viewer in meters. Range 0 to
  15000.

- elevation:

  Default `10`, vectorized. Solar elevation angle above/below the
  horizon (degrees). Range -4.2 to 90.

- visibility:

  Default `50`, vectorized. Range 20 to 131.8. Meteorological range in
  kilometers for Prague model.

- albedo:

  Default `0.5`, vectorized. Range 0 to 1. Ground albedo.

- azimuth:

  Default `90`, single value. Solar azimuth (degrees). Defaults South.

- number_cores:

  Default `1`. Number of threads to use in computation.

- wide_spectrum:

  Default `FALSE`. Whether to use the wide-spectrum (55-channel,
  polarised) coefficients.

- render_mode:

  Default `"all"`. One of `"all"`, `"atmosphere"`, or `"sun"`. Use
  `"all"` for atmosphere + solar disk, `"atmosphere"` for atmospheric
  radiance only, or `"sun"` for the solar disk only.

## Value

3-column RGB matrix.

## Examples

``` r
# Generate a basic atmosphere with the Prague model
if(run_documentation()) {
value_grid = expand.grid(
  phi = seq(0,360,by=30),
  theta = seq(0,90,by=10),
  altitude = c(0,10000)
 )
vals = calculate_sky_values(
  phi = value_grid$phi,
  theta = value_grid$theta,
  altitude = value_grid$altitude,
  elevation = 45,
  visibility = 120,
  albedo = 0
 )
 cbind(value_grid, vals)
}
#>  Coefficient file for this setting not yet present: this is a large file (2.4GB), download? [y/n] 
#> Error in check_coef_file("SkyModelDataset.dat"): Input not recognized.
```
