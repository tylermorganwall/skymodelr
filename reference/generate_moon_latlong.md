# Generate the atmosphere with the moon

Note that this is just a scaled version of
[`generate_sky()`](https://tylermorganwall.github.io/skymodelr/reference/generate_sky.md),
scaled down by the luminance of the moon as compared to the sun. This
function takes the phase of the moon into account, along with the
increase in luminosity around a full moon (known as opposition surge).
Moonlight attenuation uses the Rozenberg/Krisciunas-Schaefer airmass
approximation with configurable `moon_extinction_kV`.

## Usage

``` r
generate_moon_latlong(
  datetime,
  lat,
  lon,
  filename = NA,
  albedo = 0.5,
  turbidity = 3,
  altitude = 0,
  resolution = 2048,
  number_cores = 1,
  moon_atmosphere = FALSE,
  earthshine = TRUE,
  earthshine_albedo = 0.19,
  solar_irradiance_w_m2 = 1300,
  moon_extinction_kV = 0.172,
  hosek = TRUE,
  wide_spectrum = FALSE,
  visibility = 50,
  moon_texture_width = 801,
  moon_texture_height = 801,
  verbose = FALSE
)
```

## Arguments

- datetime:

  POSIX-compatible date-time.

- lat:

  Observer latitude (degrees N).

- lon:

  Observer longitude (degrees E; west \< 0).

- filename:

  Default `NA`. Path to the image file to write.

- albedo:

  Default `0.5`. Ground albedo, range 0 to 1.

- turbidity:

  Default `3`. Atmospheric turbidity, range 1.7 to 10 (*Hosek only*).

- altitude:

  Default `0`. Observer altitude (m), range 0 to 15000 (*Prague only*).

- resolution:

  Default `2048`. Image height in pixels (width = 2 \* height).

- number_cores:

  Default `1`. CPU threads to use.

- moon_atmosphere:

  Default `FALSE`. If `TRUE`, this generates atmospheric scattering from
  light from the moon.

- earthshine:

  Default `TRUE`. If `FALSE`, disable earthshine contribution on the
  moon.

- earthshine_albedo:

  Default `0.19`. Effective Earth albedo term used in the earthshine
  irradiance approximation.

- solar_irradiance_w_m2:

  Default `1300`. Reference solar irradiance at 1 AU (W/m^2) used to
  normalize earthshine emission intensity.

- moon_extinction_kV:

  Default `0.172`. V-band atmospheric extinction coefficient for direct
  moonlight attenuation.

- hosek:

  Default `TRUE`. `FALSE` selects the Prague model.

- wide_spectrum:

  Default `FALSE`. 55-channel Prague coefficients (altitude = 0m only).

- visibility:

  Default `50`. Meteorological range (km); *Prague only*.

- moon_texture_width:

  Default `801`. Internal moon-texture render width.

- moon_texture_height:

  Default `801`. Internal moon-texture render height.

- verbose:

  Default `FALSE`. Whether to print progress bars/diagnostic info.

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
# Moonlit sky (Hosek), mid-evening in DC
if(run_documentation()) {
generate_moon_latlong(
  datetime   = as.POSIXct("2025-09-05 19:30:00",tz="America/New_York"),
  lat        = 38.9072,
  lon        = -77.0369,
  resolution = 400,
  turbidity  = 3,
  verbose    = TRUE
) |>
  rayimage::render_exposure(15) |>
  rayimage::plot_image()
}
#> phi=0.3991, earth_phase=2.7425, E_em=0.002508 W/m^2, emission_intensity=1.929e-06
#> Moon: 8.1 elevation, 120.6 azimuth, 22.809 phase, 0.064058 lux
```
