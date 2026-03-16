# Compute sun/moon positions via Swiss Ephemeris

Produce a data frame of apparent positions and magnitudes for the
observer.

## Usage

``` r
swe_dirs_topo_sunmoon_df(datetime, lat, lon, elev_m = 0)
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
