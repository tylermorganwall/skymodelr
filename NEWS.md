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
