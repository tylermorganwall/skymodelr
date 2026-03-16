# Compute luminous efficacy for a blackbody SPD

Compute luminous efficacy for a blackbody SPD

## Usage

``` r
luminous_efficacy_blackbody(
  temperature_K,
  wavelengths_nm = NULL,
  V_lambda = NULL
)
```

## Arguments

- temperature_K:

  Numeric vector of blackbody temperatures in Kelvin.

- wavelengths_nm:

  Optional wavelength grid in nm.

- V_lambda:

  Optional photopic V(lambda) samples matching the grid.

## Value

Numeric vector of K_eff values in lm/W.
