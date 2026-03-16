# Sample sun luminance at the disk center.

Evaluate the Prague or Hosek sky model at the sun center and integrate
against the CIE Y curve to return a luminance value. This is useful for
relative attenuation (e.g., comparing zenith sun to a low sun).

## Usage

``` r
calculate_sun_brightness(
  elevation = 10,
  azimuth = 90,
  albedo = 0.5,
  turbidity = 3,
  altitude = 0,
  visibility = 50,
  hosek = TRUE,
  wide_spectrum = FALSE,
  lambda_nm = NULL
)
```

## Arguments

- elevation:

  Default `10`. Solar elevation angle above/below the horizon (degrees).
  Range `[-4.2, 90]` for Prague, `[0, 90]` for Hosek.

- azimuth:

  Default `90`. Solar azimuth (degrees).

- albedo:

  Default `0.5`. Ground albedo, range 0 to 1.

- turbidity:

  Default `3`. Atmospheric turbidity, range 1.7 to 10 (*Hosek only*).

- altitude:

  Default `0`. Observer altitude (m), range 0 to 15000 (*Prague only*).

- visibility:

  Default `50`. Meteorological range (km); *Prague only*.

- hosek:

  Default `TRUE`. `FALSE` selects the Prague model.

- wide_spectrum:

  Default `FALSE`. Use wide-spectrum (55-channel) coefficients for
  Prague at sea level only.

- lambda_nm:

  Optional vector of wavelengths for Hosek sampling.

## Value

Numeric scalar of sun luminance (CIE Y, relative scale).

## Examples

``` r
if(run_documentation()) {
  calculate_sun_brightness(elevation = 45, hosek = TRUE)
}
#> [1] 2237321
```
