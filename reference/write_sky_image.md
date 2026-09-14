# Write a skymodelr image with EXR color metadata

Write a skymodelr image with EXR color metadata

## Usage

``` r
write_sky_image(image, filename, ...)
```

## Arguments

- image:

  Image array or rayimage image.

- filename:

  Destination image path.

- ...:

  Additional arguments passed to
  [`rayimage::ray_write_image()`](https://www.rayimage.dev/reference/ray_write_image.html).
