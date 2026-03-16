# Generate the moon into a lat-long image array patch

Produce a rasterised moon texture aligned for composition into the sky
dome. The Earth phase is inferred from the Sun-Moon geometry, and
earthshine is implemented as an emissive term scaled relative to solar
irradiance.

## Usage

``` r
generate_moon_image_latlong(
  datetime,
  lat,
  lon,
  elev_m = 0,
  width = 800,
  height = 800,
  earthshine = TRUE,
  earthshine_albedo = 0.19,
  solar_irradiance_w_m2 = 1300,
  moon_extinction_kV = 0.172,
  verbose = FALSE
)
```

## Arguments

- datetime:

  POSIXct observation time (UTC).

- lat:

  Latitude in degrees.

- lon:

  Longitude in degrees.

- elev_m:

  Observer elevation above sea level in metres.

- width:

  Output width in pixels.

- height:

  Output height in pixels.

- earthshine:

  Default `TRUE`. If `FALSE`, skip the earthshine emissive term.

- earthshine_albedo:

  Default `0.19`. Effective Earth albedo term used in the earthshine
  irradiance approximation.

- solar_irradiance_w_m2:

  Default `1300`. Reference solar irradiance at 1 AU (W/m^2) used to
  normalize earthshine emission intensity.

- moon_extinction_kV:

  Default `0.172`. V-band atmospheric extinction coefficient for direct
  moonlight attenuation.

- verbose:

  Default `FALSE`. If `TRUE`, print earthshine diagnostics.
