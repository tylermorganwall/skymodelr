# skymodelr (development version)

* Keep Moon disk radiance and color continuous as its center crosses the
  horizon. Disk generation retains the finite horizon attenuation/tint and
  leaves partial visibility to the renderer, so the upper limb remains visible
  until it sets. The legacy `generate_moon_latlong()` path now does the same:
  attenuation defaults to no center cutoff, and normalization uses the complete
  disk before clipping so the remaining limb does not become artificially bright.
  Disk visibility is independent of the optional atmosphere model's domain.

* Add `generate_sun_disk()` and `generate_moon_disk()` for detailed celestial
  light textures without an intermediate sky map. Both return linear RGB
  radiance, true-north azimuth, elevation, apparent diameter, and an explicit
  rectilinear projection. Lunar edge coverage is applied once; no alpha channel
  is returned. The Moon reuses existing phase, surface texture, earthshine,
  photometry, and atmospheric extinction routines.

# skymodelr 0.3.0

* Added EXR metadata tagging for generated sky maps. EXR output is tagged with
  sRGB/Rec.709 chromaticities and an adopted neutral white, defaulting to D60.
  This records the intended daylight neutral without converting pixel values.
  This complements the Prague RGB tint correction by recording the intended
  D60 adopted neutral in EXR metadata.
* Added `prague_rgb_correction`, a Prague-only linear RGB output correction
  that reduces the small magenta / negative-green tint observed in generated
  Prague RGB sky maps. The correction is enabled by default for Prague RGB
  output and can be disabled with `prague_rgb_correction = FALSE` to recover
  raw Prague model RGB values.
* Prepared package metadata and source packaging for CRAN release.
* Switched to single GPL-3 package license metadata.
* Added helpers to list and clear cached sky data downloads.
* Improved bundled data and third-party provenance notices.
* Added minimal CRAN-safe tests.
