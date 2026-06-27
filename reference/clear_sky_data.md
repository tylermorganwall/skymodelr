# Clear cached sky data files

Removes files cached by
[`download_sky_data()`](https://tylermorganwall.github.io/skymodelr/reference/download_sky_data.md)
under `tools::R_user_dir("skymodelr", "data")`.

## Usage

``` r
clear_sky_data(files = NULL, ask = interactive())
```

## Arguments

- files:

  Default `NULL`. Character vector of cached file basenames to remove.

- ask:

  Default [`interactive()`](https://rdrr.io/r/base/interactive.html).
  Whether to ask for confirmation before deleting files.

## Value

Invisibly, the paths successfully removed.

## Examples

``` r
if (FALSE) { # interactive()
clear_sky_data()
}
```
