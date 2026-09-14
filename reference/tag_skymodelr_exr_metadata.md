# Tag skymodelr EXR metadata

Tag skymodelr EXR metadata

## Usage

``` r
tag_skymodelr_exr_metadata(
  sky,
  adopted_white_xy = .skymodelr_d60_xy,
  adopted_white_name = .skymodelr_exr_adopted_neutral_label,
  rgb_colorspace = rayimage::CS_SRGB,
  model_name = NULL,
  prague_rgb_correction = NULL
)
```

## Arguments

- sky:

  Sky image array.

- adopted_white_xy:

  Default `.skymodelr_d60_xy`. Adopted white xy.

- adopted_white_name:

  Default `.skymodelr_exr_adopted_neutral_label`. Adopted white name.

- rgb_colorspace:

  Default
  [`rayimage::CS_SRGB`](https://www.rayimage.dev/reference/colorspace_descriptors.html).
  RGB colorspace.

- model_name:

  Default `NULL`. Sky model name.

- prague_rgb_correction:

  Default `NULL`. Prague RGB correction label.
