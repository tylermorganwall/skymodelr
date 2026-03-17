# Generate a location and time-specific sky dome (optionally with stars)

Convenience wrapper around
[`generate_sky()`](https://tylermorganwall.github.io/skymodelr/reference/generate_sky.md)
that:

1.  Computes the Sun's apparent position for `datetime`, `lat`, and
    `lon` (via Swiss Ephemeris / **swephR**).

2.  Renders the corresponding sky model.

3.  Optionally overlays a star field using
    [`generate_stars()`](https://tylermorganwall.github.io/skymodelr/reference/generate_stars.md)
    or moon with
    [`generate_moon_latlong()`](https://tylermorganwall.github.io/skymodelr/reference/generate_moon_latlong.md).

## Usage

``` r
generate_sky_latlong(
  datetime,
  lat,
  lon,
  filename = NA,
  albedo = 0.5,
  turbidity = 3,
  altitude = 0,
  resolution = 2048,
  number_cores = 1,
  hosek = TRUE,
  wide_spectrum = FALSE,
  visibility = 50,
  stars = FALSE,
  star_width = 1,
  planets = FALSE,
  moon = FALSE,
  moon_atmosphere = FALSE,
  moon_hosek = TRUE,
  render_mode = "all",
  below_horizon = TRUE,
  verbose = FALSE,
  stars_exposure = 0,
  ...
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

  Default `NA`. Path to the `.exr` file to write.

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

- hosek:

  Default `TRUE`. `FALSE` selects the Prague model.

- wide_spectrum:

  Default `FALSE`. 55-channel Prague coefficients (altitude = 0m only).

- visibility:

  Default `50`. Meteorological range (km); *Prague only*.

- stars:

  Default `FALSE`. If `TRUE`, composite a star field using
  [`generate_stars()`](https://tylermorganwall.github.io/skymodelr/reference/generate_stars.md).
  Automatically falls back to a pure star render when the black-sky
  condition is met (see *Details*).

- star_width:

  Default `1`. Passed to
  [`generate_stars()`](https://tylermorganwall.github.io/skymodelr/reference/generate_stars.md)
  to control stellar point-spread size.

- planets:

  Default `FALSE`. If `TRUE`, composite bright planets via
  [`generate_planets()`](https://tylermorganwall.github.io/skymodelr/reference/generate_planets.md).

- moon:

  Default `FALSE`. If `TRUE`, overlay a moon render from
  [`generate_moon_latlong()`](https://tylermorganwall.github.io/skymodelr/reference/generate_moon_latlong.md).

- moon_atmosphere:

  Default `FALSE`. If `TRUE`, this generates atmospheric scattering from
  light from the moon.

- moon_hosek:

  Default `TRUE`. Whether to use the faster (but less accurate) Hosek
  model for atmospheric scattering from the moon. Note that the light
  scattered from the moon is much less intense than the sun, and thus
  small inaccuracies are much less likely to be noticable.

- render_mode:

  Default `"all"`. One of `"all"`, `"atmosphere"`, or `"sun"`. Use
  `"all"` for atmosphere + solar disk, `"atmosphere"` for atmospheric
  radiance only, or `"sun"` for the solar disk only.

- below_horizon:

  Default `TRUE`. Whether to sample atmospheric scattering below the
  horizon, which is non-zero when altitude \> 0.

- verbose:

  Default `FALSE`. Whether to print progress bars/diagnostic info.

- stars_exposure:

  Default `0`. Increases star exposure by `2^exposure`. Non-physical,
  this just controls adjustments for artistic effect.

- ...:

  Additional **named** arguments forwarded to
  [`generate_stars()`](https://tylermorganwall.github.io/skymodelr/reference/generate_stars.md),
  and when enabled,
  [`generate_planets()`](https://tylermorganwall.github.io/skymodelr/reference/generate_planets.md)
  and
  [`generate_moon_latlong()`](https://tylermorganwall.github.io/skymodelr/reference/generate_moon_latlong.md).
  Unmatched names are ignored.

## Value

Either the raw data, or the data is invisibly returned if filename is
given. The EXR is written to `filename`.

## Details

*Solar angles* - altitude (degrees above the horizon) and azimuth
(degrees clockwise from east, so 90 degrees = south) - are derived
internally; you never have to supply them directly.

*Black-sky rule* - With the Prague model the sky radiance is defined
only down to -4.2 degrees, and with the Hosek model it is defined only
about 0 degrees. Below that the function skips the sky render and writes
**only stars** when `stars = TRUE`.

## Examples

``` r
# Morning sunrise on spring solstice over Washington, DC with Prague model
if(run_documentation()) {
generate_sky_latlong(
  datetime    = as.POSIXct("2025-03-21 06:15:00",tz="EST"),
  lat         = 38.9072,
  lon         = -77.0369,
  number_cores = 2,
  hosek = FALSE
) |>
  rayimage::plot_image()
}
#>  Coefficient file for this setting not yet present: this is a large file (107MB), download? [y/n] 
#> Error in check_coef_file("SkyModelDatasetGround.dat"): Input not recognized.
if(run_documentation()) {
generate_sky_latlong(
  datetime    = as.POSIXct("2025-03-21 12:00:00",tz="EST"),
  lat         = 38.9072,
  lon         = -77.0369,
  number_cores = 2,
) |>
  rayimage::plot_image()
}

if(run_documentation()) {
generate_sky_latlong(
  datetime    = as.POSIXct("2025-03-21 18:00:00",tz="EST"),
  lat         = 38.9072,
  lon         = -77.0369,
  number_cores = 2,
) |>
  rayimage::plot_image()
}

if(run_documentation()) {
generate_sky_latlong(
  datetime    = as.POSIXct("2025-03-21 18:30:00",tz="EST"),
  lat         = 38.9072,
  lon         = -77.0369,
  number_cores = 2,
  hosek = FALSE,
  verbose=TRUE,
) |>
  rayimage::render_exposure(exposure=2) |>
  rayimage::plot_image()
}
#> Sun: -2.7 elevation, 273.0 azimuth
#>  Coefficient file for this setting not yet present: this is a large file (107MB), download? [y/n] 
#> Error in check_coef_file("SkyModelDatasetGround.dat"): Input not recognized.
```
