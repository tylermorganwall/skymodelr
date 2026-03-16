# Compute topocentric moon/sun directions

Query Swiss Ephemeris for topocentric sun and moon direction vectors and
magnitudes for a given observer location.

## Usage

``` r
swe_dirs_topo_moon_sun(
  datetime,
  lat,
  lon,
  elev_m = 0,
  moon_extinction_kV = 0.172
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

- moon_extinction_kV:

  Default `0.172`. V-band atmospheric extinction coefficient used in the
  Rozenberg/Krisciunas-Schaefer airmass attenuation model for moonlight.
