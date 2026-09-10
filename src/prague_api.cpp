// Registered, versioned native interface. GPL-3.
#include "../inst/include/skymodelr/prague_api.h"
#include "PragueSkyModel/PragueSkyModel.h"
#include <R_ext/Rdynload.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <memory>
#include <stdexcept>

struct skymodelr_prague_handle { PragueSkyModel model; };
namespace {
void error_message(char *output, size_t capacity, const char *message) noexcept {
    if (!output || !capacity) return;
    size_t n = std::min(capacity - 1, std::strlen(message));
    std::memcpy(output, message, n); output[n] = '\0';
}
template <typename F> int guarded(F fn, char *error, size_t capacity) noexcept {
    try { fn(); return 1; }
    catch (const std::exception &e) { error_message(error, capacity, e.what()); }
    catch (...) { error_message(error, capacity, "Unknown error in skymodelr Prague model."); }
    return 0;
}
const PragueSkyModel &model(const skymodelr_prague_handle *handle) {
    if (!handle) throw std::invalid_argument("Prague model is not initialized.");
    return handle->model;
}
void validate_finite(double value) {
    if (!std::isfinite(value)) throw std::invalid_argument("Nonfinite Prague query parameter.");
}
void atmosphere(double visibility, double albedo) {
    validate_finite(visibility); validate_finite(albedo);
    if (visibility < 0 || albedo < 0 || albedo > 1)
        throw std::invalid_argument("Invalid Prague visibility or ground albedo.");
}
PragueSkyModel::Parameters parameters(const skymodelr_prague_parameters *p) {
    if (!p) throw std::invalid_argument("Missing Prague query parameters.");
    for (double angle : {p->theta, p->gamma, p->shadow, p->zero}) {
        validate_finite(angle);
        if (angle < 0 || angle > 3.14159265358979323846)
            throw std::invalid_argument("Prague direction angles must lie between zero and pi.");
    }
    validate_finite(p->elevation); validate_finite(p->altitude); atmosphere(p->visibility, p->albedo);
    if (p->altitude < 0) throw std::invalid_argument("Prague altitude must be nonnegative.");
    return {p->theta, p->gamma, p->shadow, p->zero, p->elevation,
            p->altitude, p->visibility, p->albedo};
}
}

extern "C" {
static skymodelr_prague_handle *create_model(const char *filename, double visibility,
                                                          int cache_spectra, int transmission_table,
                                                          double max_mib,
                                                          char *error, size_t capacity) noexcept {
    std::unique_ptr<skymodelr_prague_handle> result;
    if (!guarded([&] {
        if (!filename || !*filename) throw std::invalid_argument("Missing Prague coefficient filename.");
        validate_finite(visibility);
        if (visibility < 0) throw std::invalid_argument("Prague visibility must be nonnegative.");
        if ((cache_spectra != 0 && cache_spectra != 1) ||
            (transmission_table != 0 && transmission_table != 1))
            throw std::invalid_argument("Prague cache options must be zero or one.");
        if (std::isnan(max_mib) || max_mib < 0)
            throw std::invalid_argument("Prague transmission table limit must be nonnegative MiB or infinity.");
        result.reset(new skymodelr_prague_handle);
        result->model.initialize(filename, visibility, cache_spectra != 0, transmission_table != 0, max_mib);
    }, error, capacity)) return nullptr;
    return result.release();
}
static void destroy_model(skymodelr_prague_handle *handle) noexcept { delete handle; }
static size_t memory_usage(const skymodelr_prague_handle *handle) noexcept {
    return handle ? handle->model.memoryUsage() : 0;
}
static int available(const skymodelr_prague_handle *handle, skymodelr_prague_available *out,
                      char *error, size_t capacity) noexcept {
    return guarded([&] {
        if (!out) throw std::invalid_argument("Missing Prague metadata output.");
        auto a = model(handle).getAvailableData();
        *out = {a.albedoMin, a.albedoMax, a.altitudeMin, a.altitudeMax,
                a.elevationMin, a.elevationMax, a.visibilityMin, a.visibilityMax,
                a.polarisation ? 1 : 0, a.channels, a.channelStart, a.channelWidth};
    }, error, capacity);
}
static int compute_parameters(const skymodelr_prague_handle *handle,
                               const double *position, const double *direction,
                               double elevation, double azimuth, double visibility, double albedo,
                               skymodelr_prague_parameters *out, char *error, size_t capacity) noexcept {
    return guarded([&] {
        if (!position || !direction || !out) throw std::invalid_argument("Missing Prague position, direction, or output.");
        validate_finite(elevation); validate_finite(azimuth); atmosphere(visibility, albedo);
        for (int i = 0; i < 3; ++i) { validate_finite(position[i]); validate_finite(direction[i]); }
        double magnitude = std::hypot(std::hypot(direction[0], direction[1]), direction[2]);
        if (!(magnitude > 0) || !std::isfinite(magnitude))
            throw std::invalid_argument("Prague view direction must have finite nonzero length.");
        double radius = std::hypot(std::hypot(position[0], position[1]), position[2] + 6378000.0);
        if (!std::isfinite(radius) || radius < 6378000.0 - 1e-6)
            throw std::invalid_argument("Prague viewpoint lies inside the Earth.");
        auto p = model(handle).computeParameters({position[0], position[1], position[2]},
            {direction[0], direction[1], direction[2]}, elevation, azimuth, visibility, albedo);
        *out = {p.theta, p.gamma, p.shadow, p.zero, p.elevation, p.altitude, p.visibility, p.albedo};
    }, error, capacity);
}
static int spectrum(const skymodelr_prague_handle *handle, int quantity,
                     const skymodelr_prague_parameters *input, const double *wavelengths,
                     size_t count, double distance, int attenuate_sun, double *out,
                     char *error, size_t capacity) noexcept {
    return guarded([&] {
        const auto &m = model(handle);
        auto p = parameters(input);
        if (count && (!wavelengths || !out)) throw std::invalid_argument("Missing Prague spectral input or output.");
        for (size_t i = 0; i < count; ++i) validate_finite(wavelengths[i]);
        switch (quantity) {
        case SKYMODELR_PRAGUE_SKY: m.skyRadianceSpectrum(p, wavelengths, count, out); break;
        case SKYMODELR_PRAGUE_TRANSMITTANCE:
            if (std::isnan(distance) || distance < 0) throw std::invalid_argument("Invalid Prague transmission distance.");
            m.transmittanceSpectrum(p, wavelengths, count,
                std::min(distance, std::numeric_limits<double>::max()), out);
            break;
        case SKYMODELR_PRAGUE_SUN:
            if (attenuate_sun != 0 && attenuate_sun != 1) throw std::invalid_argument("attenuate_sun must be zero or one.");
            for (size_t i = 0; i < count; ++i) out[i] = m.sunRadiance(p, wavelengths[i], attenuate_sun != 0);
            break;
        case SKYMODELR_PRAGUE_POLARISATION:
            for (size_t i = 0; i < count; ++i) out[i] = m.polarisation(p, wavelengths[i]);
            break;
        default: throw std::invalid_argument("Unknown Prague spectral quantity.");
        }
    }, error, capacity);
}
static const skymodelr_prague_api *get_api() noexcept {
    static const skymodelr_prague_api api = {SKYMODELR_PRAGUE_ABI_VERSION, sizeof(skymodelr_prague_api),
        create_model, destroy_model, available, memory_usage, compute_parameters, spectrum};
    return &api;
}
}

// [[Rcpp::init]]
void skymodelr_register_prague_api(DllInfo *dll) {
    (void)dll;
    R_RegisterCCallable("skymodelr", "prague_get_api", reinterpret_cast<DL_FUNC>(&get_api));
}
