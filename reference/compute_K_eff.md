# Compute effective luminous efficacy for a chosen SPD

Compute effective luminous efficacy for a chosen SPD

## Usage

``` r
compute_K_eff(
  spd_type = c("D65", "BB5778"),
  wavelength_grid = NULL,
  V_lambda = NULL
)
```

## Arguments

- spd_type:

  SPD shape. Either `"D65"` or `"BB5778"`.

- wavelength_grid:

  Optional numeric vector of wavelengths (nm).

- V_lambda:

  Optional numeric vector of photopic V(lambda) samples.

## Value

Effective luminous efficacy in lm/W.
