#ifndef SKYMODELR_PRAGUE_METAL_BACKEND_H
#define SKYMODELR_PRAGUE_METAL_BACKEND_H

#include "PragueSkyModel/PragueSkyModel.h"
#include <cstddef>
#include <string>

// Borrowed, equally sized columns, validated by the R entry point. Angles are
// in degrees, wavelengths in nm, altitude in metres, visibility in kilometres.
struct PragueRadianceColumns {
  const double *phi, *theta, *wavelength, *elevation;
  const double *albedo, *altitude, *visibility, *azimuth;
};

bool prague_metal_available();
void prague_metal_radiance(const PragueSkyModel &model,
                           const std::string &dataset,
                           const PragueRadianceColumns &columns,
                           size_t count, bool atmosphere, bool sun,
                           double *output, void (*check_interrupt)());

#endif
