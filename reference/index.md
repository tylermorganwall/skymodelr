# Package index

## Generate Sky EXRs

Functions for generating different sky models

- [`generate_moon_latlong()`](https://tylermorganwall.github.io/skymodelr/reference/generate_moon_latlong.md)
  : Generate the atmosphere with the moon

- [`generate_planets()`](https://tylermorganwall.github.io/skymodelr/reference/generate_planets.md)
  : Generate a bright-planet image array

- [`generate_sky()`](https://tylermorganwall.github.io/skymodelr/reference/generate_sky.md)
  : Generate a Hosek-Wilkie sky dome array

- [`generate_sky_latlong()`](https://tylermorganwall.github.io/skymodelr/reference/generate_sky_latlong.md)
  : Generate a location and time-specific sky dome (optionally with
  stars)

- [`generate_stars()`](https://tylermorganwall.github.io/skymodelr/reference/generate_stars.md)
  :

  Generate a star‑field array aligned with
  [`generate_sky()`](https://tylermorganwall.github.io/skymodelr/reference/generate_sky.md)

## Downloading

Download model data

- [`download_sky_data()`](https://tylermorganwall.github.io/skymodelr/reference/download_sky_data.md)
  : Download Prague Sky Model Coefficient Data

## Sample Values

Functions to sample specific values from the atmospheric model

- [`calculate_sky_values()`](https://tylermorganwall.github.io/skymodelr/reference/calculate_sky_values.md)
  : Sample a direction from the Prague model.
- [`calculate_sun_brightness()`](https://tylermorganwall.github.io/skymodelr/reference/calculate_sun_brightness.md)
  : Sample sun luminance at the disk center.

## pkgdown utils

Internal but must be exported to meet CRAN guidelines

- [`run_documentation()`](https://tylermorganwall.github.io/skymodelr/reference/run_documentation.md)
  : Run Documentation
