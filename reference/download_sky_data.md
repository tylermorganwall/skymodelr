# Download Prague Sky Model Coefficient Data

This model which allows for sun angles below the horizon, wide spectral
ranges including infrared, polarization, as well rendering at altitude)
need to download a coefficient file. There are three versions of this
dataset, listed here in order of increasing size:

## Usage

``` r
download_sky_data(sea_level = TRUE, wide_spectrum = FALSE)
```

## Arguments

- sea_level:

  Default `TRUE`. Download the sea‑level–only data. Set to `FALSE` to
  download the full‑altitude dataset.

- wide_spectrum:

  Default `FALSE`. If `TRUE`, downloads the wide‑spectrum (55‑channel,
  polarised) version. Valid only when `sea_level = TRUE`.

## Value

Invisibly, the full path to the data file.

## Details

|                                           |                                        |       |
|-------------------------------------------|----------------------------------------|-------|
| Argument combination                      | File                                   | Size  |
| `sea_level = TRUE, wide_spectrum = FALSE` | `SkyModelDatasetGround.dat`            | 107MB |
| `sea_level = TRUE, wide_spectrum = TRUE`  | `PragueSkyModelDatasetGroundInfra.dat` | 574MB |
| `sea_level = FALSE`                       | `SkyModelDataset.dat`                  | 2.4GB |

## Examples

``` r
if (FALSE) { # \dontrun{
  # Standard (11‑channel, sea‑level) coefficients
  download_sky_data()

  # Wide‑spectrum sea‑level coefficients
  download_sky_data(wide_spectrum = TRUE)

  # Full altitude‑range coefficients
  download_sky_data(sea_level = FALSE)
} # }
```
