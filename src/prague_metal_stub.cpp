#include "prague_metal_backend.h"
#include <stdexcept>

bool prague_metal_available() { return false; }

void prague_metal_radiance(const PragueSkyModel &, const std::string &,
                           const PragueRadianceColumns &, size_t, bool, bool,
                           double *, void (*)()) {
  throw std::runtime_error(
      "backend = \"gpu\" requires macOS with an available Metal device. "
      "Use backend = \"cpu\" on this system.");
}
