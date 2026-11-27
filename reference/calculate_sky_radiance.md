# Sample Prague sky radiance over directions and wavelengths

Evaluate the Prague spectral sky model at arbitrary spherical directions
and wavelengths, including grids of atmospheric conditions. Use
`render_mode = "all"` for atmosphere + solar disk, `"atmosphere"` for
atmospheric radiance only, or `"sun"` for the solar disk only.

## Usage

``` r
calculate_sky_radiance(
  phi,
  theta,
  lambda_nm,
  altitude = 0,
  elevation = 10,
  visibility = 50,
  albedo = 0.5,
  azimuth = 90,
  number_cores = 1,
  wide_spectrum = FALSE,
  render_mode = "all",
  backend = "cpu"
)
```

## Arguments

- phi:

  Azimuthal angle of the sample, degrees. Vectorized. Range 0 to 360.

- theta:

  Vertical angle of the sample, degrees. Vectorized. Range -90 to 90.

- lambda_nm:

  Wavelength in nanometers. Vectorized. Must lie within the loaded
  Prague dataset range.

- altitude:

  Default `0`. Vectorized altitude of the viewer in meters. Range 0 to
  15000.

- elevation:

  Default `10`. Vectorized solar elevation angle above/below the horizon
  (degrees). Range -4.2 to 90.

- visibility:

  Default `50`. Vectorized meteorological range in kilometers. Range 20
  to 131.8.

- albedo:

  Default `0.5`. Vectorized ground albedo. Range 0 to 1.

- azimuth:

  Default `90`. Vectorized solar azimuth (degrees). Defaults South.

- number_cores:

  Default `1`. Number of CPU threads. Ignored with `backend = "gpu"`.

- wide_spectrum:

  Default `FALSE`. Whether to use the wide-spectrum (55-channel,
  polarised) coefficients.

- render_mode:

  Default `"all"`. One of `"all"`, `"atmosphere"`, or `"sun"`. Use
  `"all"` for atmosphere + solar disk, `"atmosphere"` for atmospheric
  radiance only, or `"sun"` for the solar disk only.

- backend:

  Default `"cpu"`. Use `"cpu"` for CPU evaluation or `"gpu"` for batch
  evaluation on a macOS Metal device.

## Value

Numeric vector of spectral radiance values in W m^-2 sr^-1 nm^-1, in
input row order.

## Details

Inputs are paired element by element, with shorter vectors recycled by
[`data.frame()`](https://rdrr.io/r/base/data.frame.html). Use
[`expand.grid()`](https://rdrr.io/r/base/expand.grid.html) to form all
combinations, then pass its columns as the corresponding arguments. The
function does not create a Cartesian product automatically.

The GPU backend evaluates geometry, atmospheric radiance, and solar
radiance using the same Prague coefficients as the CPU backend. It uses
single precision, so results can differ slightly, especially near the
solar limb. Coefficients and compiled shaders are cached across calls;
changing datasets replaces the cached GPU coefficients. Large inputs are
processed in bounded batches. The initial upload and shader compilation
can outweigh the benefit for small inputs. An unavailable Metal device
or a GPU failure raises an error; use `backend = "cpu"` to evaluate on
the CPU.

## Examples

``` r
# Compare the same direction and wavelength at two observer altitudes.
lambda_vals = calculate_sky_radiance(
  phi = c(90, 90),
  theta = c(45, 45),
  lambda_nm = 550,
  altitude = c(0, 10000),
  elevation = 20,
  visibility = 80,
  albedo = 0.1
)
#> Error: The Prague sky model coefficient file 'SkyModelDataset.dat' is not installed. Run download_sky_data(sea_level = FALSE, wide_spectrum = FALSE) explicitly to download it.
cbind(altitude = c(0, 10000), radiance = lambda_vals)
#> Error: object 'lambda_vals' not found

# Form a direction-by-wavelength grid and preserve its row ordering.
# On a Mac with Metal, add backend = "gpu" to evaluate the grid on the GPU.
grid = expand.grid(
  phi = seq(0, 350, by = 10),
  theta = seq(0, 90, by = 10),
  lambda_nm = seq(380, 740, by = 20)
)
grid$radiance = calculate_sky_radiance(
  phi = grid$phi,
  theta = grid$theta,
  lambda_nm = grid$lambda_nm,
  render_mode = "atmosphere"
)
#> Error: The Prague sky model coefficient file 'SkyModelDatasetGround.dat' is not installed. Run download_sky_data(sea_level = TRUE, wide_spectrum = FALSE) explicitly to download it.
```
