# Star catalog used for sky rendering

A data frame of stellar positions and photometric/color information used
by
[`generate_stars()`](https://tylermorganwall.github.io/skymodelr/reference/generate_stars.md)
to render star fields.

## Usage

``` r
data(stars)
```

## Format

A tibble with 9,110 rows and 8 variables:

- bsc_number:

  Numeric identifier from the source catalog.

- ra_rad:

  Right ascension in radians.

- dec_rad:

  Declination in radians.

- v_mag:

  Apparent visual magnitude (V band).

- spec:

  Spectral type string.

- r:

  Relative red channel weight derived from spectral type.

- g:

  Relative green channel weight derived from spectral type.

- b:

  Relative blue channel weight derived from spectral type.

## Source

Bright Star Catalogue, 5th Revised Edition (BSC5), stored in
`data-raw/BSC5`.
