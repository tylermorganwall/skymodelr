# Sample Prague sky radiance at a chosen wavelength.

Evaluate the Prague spectral sky model at arbitrary spherical directions
and return radiance at a user-specified wavelength. Use
`render_mode = "all"` for atmosphere + solar disk, `"atmosphere"` for
atmospheric radiance only, or `"sun"` for the solar disk only.

## Usage

``` r
calculate_sky_radiance(
  phi,
  theta,
  lambda_nm,
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

  Azimuthal angle of the sample, degrees. Vectorized. Range 0 to 360.

- theta:

  Vertical angle of the sample, degrees. Vectorized. Range -90 to 90.

- lambda_nm:

  Wavelength in nanometers. Vectorized. Must lie within the loaded
  Prague dataset range.

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

Numeric vector of radiance values at the requested wavelength(s).

## Examples

``` r
if(run_documentation()) {
lambda_vals = calculate_sky_radiance(
  phi = c(90, 90),
  theta = c(45, 45),
  lambda_nm = 550,
  altitude = c(0, 10000),
  elevation = 20,
  visibility = 80,
  albedo = 0.1
)
cbind(altitude = c(0, 10000), radiance = lambda_vals)
}
#>  Coefficient file for this setting not yet present: this is a large file (2.4GB), download? [y/n] 
#> Error in check_coef_file("SkyModelDataset.dat"): Input not recognized.
```
