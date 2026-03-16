# Compute a unit RGB vector for a chosen SPD

Compute a unit RGB vector for a chosen SPD

## Usage

``` r
compute_spd_rgb_unit(spd_type = c("D65", "BB5778"))
```

## Arguments

- spd_type:

  SPD shape. Either `"D65"` or `"BB5778"`.

## Value

Numeric vector length 3, normalized so sum equals 1.
