# Compute radiometric star amplitude and RGB unit vector

Compute radiometric star amplitude and RGB unit vector

## Usage

``` r
star_radiometric_amplitude(mV, bv, wavelengths_nm = NULL, V_lambda = NULL)
```

## Arguments

- mV:

  Numeric vector of V magnitudes.

- bv:

  Numeric vector of B-V values.

- wavelengths_nm:

  Optional wavelength grid in nm.

- V_lambda:

  Optional photopic V(lambda) samples matching the grid.

## Value

List with E_v (lux), E_e (W/m^2), K_eff (lm/W), and rgb_unit.
