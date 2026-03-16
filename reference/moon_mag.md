# Moon V-band magnitude from phase angle

Moon V-band magnitude from phase angle

## Usage

``` r
moon_mag(phase_deg, dist_km = 384400, mean_km = 384400)
```

## Arguments

- phase_deg:

  Phase angle ψ in degrees (0 = full Moon, 180 = new).

- dist_km:

  Default `384400`. Current geocentric distance (km).

- mean_km:

  Default `384400`. Mean distance (km).

## Value

Apparent V magnitude `m`.
