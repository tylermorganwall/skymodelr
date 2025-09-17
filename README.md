
# skymodelr

<img src="man/figures/startrailscolor.jpg"></img>

## Overview

**skymodelr** generates physically‑plausible sky domes and night skies
as high‑dynamic‑range EXR images/arrays directly from R. It implements
the Hosek–Wilkie analytic sky model and (optionally) the 2021–22 Prague
spectral sky model (below‑horizon sun, altitude, and wide‑spectrum
support). It also includes tools to add the moon as well as accurate
visible **star fields** aligned to observer location/time. Outputs are
lat‑long environment maps (2:1 equirectangular) that you can feed into
renderers (such as *rayrender*).

`generate_sky_latlong()` composes a full sky environment using the
functions `generate_sky()`, `generate_moon_latlong()`, and
`generate_stars()` for a specific latitude, longitude, and time. By
default `generate_sky_latlong()` only includes the sun’s contribution,
but you can also include stars and the moon by setting `moon = TRUE` and
`stars = TRUE`.

## Installation

``` r
# Latest version from GitHub
remotes::install_github("tylermorganwall/skymodelr")
```

If you plan to use the Prague spectral model, download its coefficient
dataset(s) once (see `download_sky_data()` below).

## Functions

- `generate_sky()` — Write/return an EXR sky dome using either the
  Hosek–Wilkie (default) or Prague model.

  - `filename = NA` to return the HDR array in‑memory (otherwise a `.exr`
    is written).
  - `albedo = 0.5` ground reflectance (0–1).
  - `turbidity = 3` atmospheric turbidity (1.7–10; Hosek only).
  - `elevation = 10`, `azimuth = 90` solar position (degrees).
  - `altitude = 0` observer altitude in meters (Prague only).
  - `resolution = 2048` image height (width is `2 * resolution`).
  - `numbercores = 1` threads.
  - `hosek = TRUE` set `FALSE` to use the Prague model; Prague options
    `wide_spectrum`, `visibility`.
  - `render_solar_disk = TRUE` set `FALSE` to render the atmospheric
    scattering effects without the solar disk.

- `generate_sky_latlong()` — Produce a complete equirectangular sky
  array/EXR. Accepts date/time and observer location, and (optionally)
  adds stars and a moon‑lit atmosphere.

  - Core args: `filename = NA`, `datetime`, `lat`, `lon`, `albedo`,
    `turbidity`, `resolution`, `numbercores`.
  - Model selection: `hosek = TRUE` (Hosek–Wilkie) or set
    `hosek = FALSE` to use the Prague spectral model; Prague extras:
    `wide_spectrum`, `visibility`, `altitude`.
  - Composition: `stars = FALSE`, `star_width = 1`, `moon = FALSE`.

- `generate_moon_latlong()` — Produce a moon‑lit atmosphere by scaling a
  sky dome to the moon’s luminance (phase + opposition surge). Computes
  the moon’s position from time/location. Arguments mirror
  `generate_sky_latlong()`.

- `generate_stars()` — Generate a star‑field EXR aligned with the sky
  dome:

  - `filename = NA`, `resolution = 2048`.
  - `lon`, `lat` observer longitude/latitude (deg) and `datetime` (used
    for local sidereal time).
  - Optional extinction/appearance controls: `turbidity`, `ozone_du`,
    `altitude`, `star_width`, `atmosphere_effects`,
    `upper_hemisphere_only`, `numbercores`.

- `calculate_sky_values()` — Sample radiance from the Prague model for
  given sky directions (`phi`, `theta`) and conditions (`elevation`,
  `altitude`, `visibility`, `albedo`, `azimuth`).

- `download_sky_data(sea_level = TRUE, wide_spectrum = FALSE)` —
  Download Prague model coefficient data:

  - Sea‑level, standard spectrum: `SkyModelDatasetGround.dat` (~107 MB)
  - Sea‑level, wide spectrum: `PragueSkyModelDatasetGroundInfra.dat`
    (~574 MB)
  - Full‑altitude dataset: `SkyModelDataset.dat` (~2.4 GB)

## Usage

Morning in DC:

``` r
library(skymodelr)
library(rayimage)

env = generate_sky_latlong(
  filename    = NA,
  datetime   = as.POSIXct("2025-03-21 06:15:00",tz="EST"),
  lat        = 38.9072,
  lon        = -77.0369,
  resolution = 800
)
rayimage::plot_image(env)
```

![](man/figures/full_sky-1.png)<!-- -->

Afternoon in DC:

``` r
env = generate_sky_latlong(
  filename    = NA,
  datetime   = as.POSIXct("2025-03-21 12:15:00",tz="EST"),
  lat        = 38.9072,
  lon        = -77.0369,
  resolution = 800
)
rayimage::plot_image(env)
```

![](man/figures/full_sky_afternoon-1.png)<!-- -->

Evening in DC:

``` r
env = generate_sky_latlong(
  filename    = NA,
  datetime   = as.POSIXct("2025-03-21 18:15:00",tz="EST"),
  lat        = 38.9072,
  lon        = -77.0369,
  resolution = 800
)
rayimage::plot_image(env)
```

![](man/figures/full_sky_evening-1.png)<!-- -->

Evening in DC (Prague model):

``` r
env = generate_sky_latlong(
  filename    = NA,
  datetime   = as.POSIXct("2025-03-21 18:15:00",tz="EST"),
  lat        = 38.9072,
  lon        = -77.0369,
  resolution = 800,
  hosek = FALSE
)
rayimage::plot_image(env)
```

![](man/figures/full_sky_evening_prague-1.png)<!-- -->

Full sun + moon + stars:

``` r
env = generate_sky_latlong(
  filename    = NA,
  datetime   = as.POSIXct("2025-03-21 18:00:00",tz="EST"),
  lat        = 38.9072,
  lon        = -77.0369,
  resolution = 800,
  stars = TRUE,
  moon = TRUE
)
rayimage::plot_image(env)
```

![](man/figures/full_sky_night-1.png)<!-- -->

## Custom sun position, multicore, no solar disk

``` r
sky = generate_sky(
  filename = NA,
  albedo = 0.2,
  turbidity = 3,
  elevation = 15,
  azimuth = 135,
  resolution = 800,
  numbercores = 2,
  hosek = TRUE,
  render_solar_disk = FALSE
)
rayimage::plot_image(sky)
```

![](man/figures/sky_basic-1.png)<!-- -->

## Moon‑lit atmosphere

``` r
moon_sky = generate_moon_latlong(
  filename   = NA,
  datetime  = as.POSIXct("2025-03-21 02:15:00",tz="EST"),
  lat       = 38.9072,
  lon       = -77.0369,
  albedo    = 0.2,
  turbidity = 3,
  resolution = 800,
  hosek = TRUE
)
#Increase exposure
moon_sky |> 
  rayimage::render_exposure(17) |> 
  rayimage::plot_image()
```

![](man/figures/sky_moon-1.png)<!-- -->

## Star field aligned to time and place

``` r
stars = generate_stars(
  filename   = NA,
  datetime  = as.POSIXct("2025-03-21 02:15:00",tz="EST"),
  lat       = 38.9072,
  lon       = -77.0369,
  resolution = 800,
  atmosphere_effects = TRUE,
  upper_hemisphere_only = TRUE
)
stars |> 
  render_exposure(5) |> 
  rayimage::plot_image()
```

![](man/figures/stars_basic-1.png)<!-- -->

Now render the entire sphere:

``` r
stars_full = generate_stars(
  filename   = NA,
  datetime  = as.POSIXct("2025-03-21 02:15:00",tz="EST"),
  lat       = 38.9072,
  lon       = -77.0369,
  resolution = 800,
  atmosphere_effects = FALSE,
  upper_hemisphere_only = FALSE
)
stars_full |> 
  render_exposure(5) |> 
  rayimage::plot_image()
```

![](man/figures/stars_basic_full-1.png)<!-- -->

## Use the Prague spectral model

``` r
# Download once (choose variant via args):
coef_path = download_sky_data(sea_level = TRUE, wide_spectrum = FALSE)

# Render with the Prague model:
sky_prague = generate_sky(
  albedo = 0.3,
  elevation = 10,
  azimuth = 90,
  altitude = 0,
  hosek = FALSE,
  wide_spectrum = FALSE,
  visibility = 50,
  resolution = 2048,
  numbercores = 4
)
plot_image(sky_prague)
```

## Acknowledgements

- Hosek–Wilkie sky model (analytic daylight model).
- Prague 2021–22 spectral sky model (wide‑spectrum, below‑horizon sun,
  altitude).
- Star positions from the Yale Bright Star Catalog.
