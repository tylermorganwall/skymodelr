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

A data.frame with 9,110 rows and 8 variables:

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

Derived from the Bright Star Catalogue, 5th Revised Edition (BSC5), also
known as the Yale Bright Star Catalogue. The source catalogue is
commonly identified as Hoffleit, D. and Warren, W. H. Jr. (1991), "The
Bright Star Catalogue, 5th Revised Ed.", Yale University Observatory,
and is distributed as machine-readable catalogue V/50 by the Centre de
Donnees astronomiques de Strasbourg (CDS):
<https://cdsarc.cds.unistra.fr/viz-bin/cat/V/50>.
