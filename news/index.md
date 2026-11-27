# Changelog

## skymodelr 0.7.2

### Bugfixes

- Celestial positions use UT1 for Swiss Ephemeris calculations instead
  of treating Terrestrial Time as UT1. This corrects Sun, Moon, and
  planet positions, including visible alignment errors during eclipses.

- Overall: Keeps the macOS Metal source out of automatic source
  discovery so Windows and Linux builds use the CPU stub without
  requiring an Objective-C runtime.

## skymodelr 0.7.1

### New features

- [`calculate_sky_radiance()`](https://tylermorganwall.github.io/skymodelr/reference/calculate_sky_radiance.md)
  Supports `backend = "gpu"` for batched Metal evaluation of direction,
  wavelength, and atmospheric-condition grids on macOS, with reusable
  GPU coefficients. The default remains CPU evaluation.

## skymodelr 0.7.0

### New features

- Overall: Provides a versioned native Metal interface for GPU spectral
  sky, sun, and transmittance queries, with reusable shader functions
  and independently owned coefficient snapshots.

### Bugfixes

- Overall: Preserves small angular distances in Prague transmittance
  queries, removing rounding jumps on sub-metre paths. Avoids reading
  beyond the final piecewise-linear coefficient knot.

## skymodelr (development version)

- Add Hosek support to `generate_sun_disk(hosek = TRUE)`, including
  turbidity and optional manual elevation/azimuth. The independently
  sampled solar disk remains finite at zenith and requires no Prague
  dataset.

- Return nighttime Sun positions from
  [`get_prague_sky_metadata()`](https://tylermorganwall.github.io/skymodelr/reference/get_prague_sky_metadata.md)
  below -4.2 degrees so native renderers can suppress solar radiance
  while retaining atmospheric transmission for celestial lights.

- Keep Moon disk radiance and color continuous as its center crosses the
  horizon. Disk generation retains the finite horizon attenuation/tint
  and leaves partial visibility to the renderer, so the upper limb
  remains visible until it sets. The legacy
  [`generate_moon_latlong()`](https://tylermorganwall.github.io/skymodelr/reference/generate_moon_latlong.md)
  path now does the same: attenuation defaults to no center cutoff, and
  normalization uses the complete disk before clipping so the remaining
  limb does not become artificially bright. Disk visibility is
  independent of the optional atmosphere model’s domain.

- Add
  [`generate_sun_disk()`](https://tylermorganwall.github.io/skymodelr/reference/generate_sun_disk.md)
  and
  [`generate_moon_disk()`](https://tylermorganwall.github.io/skymodelr/reference/generate_sun_disk.md)
  for detailed celestial light textures without an intermediate sky map.
  Both return linear RGB radiance, true-north azimuth, elevation,
  apparent diameter, and an explicit rectilinear projection. Lunar edge
  coverage is applied once; no alpha channel is returned. The Moon
  reuses existing phase, surface texture, earthshine, photometry, and
  atmospheric extinction routines.

## skymodelr 0.3.0

- Added EXR metadata tagging for generated sky maps. EXR output is
  tagged with sRGB/Rec.709 chromaticities and an adopted neutral white,
  defaulting to D60. This records the intended daylight neutral without
  converting pixel values. This complements the Prague RGB tint
  correction by recording the intended D60 adopted neutral in EXR
  metadata.
- Added `prague_rgb_correction`, a Prague-only linear RGB output
  correction that reduces the small magenta / negative-green tint
  observed in generated Prague RGB sky maps. The correction is enabled
  by default for Prague RGB output and can be disabled with
  `prague_rgb_correction = FALSE` to recover raw Prague model RGB
  values.
- Prepared package metadata and source packaging for CRAN release.
- Switched to single GPL-3 package license metadata.
- Added helpers to list and clear cached sky data downloads.
- Improved bundled data and third-party provenance notices.
- Added minimal CRAN-safe tests.
