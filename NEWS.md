# skymodelr 0.3.0

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
