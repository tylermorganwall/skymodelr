# List cached sky data files

Lists files cached by
[`download_sky_data()`](https://tylermorganwall.github.io/skymodelr/reference/download_sky_data.md)
under `tools::R_user_dir("skymodelr", "data")`.

## Usage

``` r
list_sky_data()
```

## Value

A data frame with columns `file`, `path`, `size`, and `modified`.

## Examples

``` r
list_sky_data()
#>                                   file
#> 1 PragueSkyModelDatasetGroundInfra.dat
#> 2                  SkyModelDataset.dat
#> 3            SkyModelDatasetGround.dat
#>                                                                         path
#> 1 /home/runner/.local/share/R/skymodelr/PragueSkyModelDatasetGroundInfra.dat
#> 2                  /home/runner/.local/share/R/skymodelr/SkyModelDataset.dat
#> 3            /home/runner/.local/share/R/skymodelr/SkyModelDatasetGround.dat
#>         size            modified
#> 1  573546092 2026-06-27 03:23:04
#> 2 2397370196 2026-06-27 03:24:32
#> 3  107632816 2026-06-27 03:22:45
```
