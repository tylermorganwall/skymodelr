// Public C++ convenience wrapper for skymodelr's registered C API. GPL-3.
// Contains no Prague implementation or coefficients; all evaluation executes
// in skymodelr. Construct on R's main thread, then share const queries freely.
#ifndef SKYMODELR_PRAGUE_HPP
#define SKYMODELR_PRAGUE_HPP
#include "prague_api.h"
#include <R_ext/Rdynload.h>
#include <limits>
#include <stdexcept>
#include <string>

namespace skymodelr {
class PragueSkyModel {
public:
    struct Vector3 {
        double x = 0, y = 0, z = 0;
        Vector3() = default;
        Vector3(double x, double y, double z) : x(x), y(y), z(z) {}
        Vector3 operator+(Vector3 v) const { return {x + v.x, y + v.y, z + v.z}; }
        Vector3 operator-(Vector3 v) const { return {x - v.x, y - v.y, z - v.z}; }
        Vector3 operator*(double s) const { return {x * s, y * s, z * s}; }
        Vector3 operator/(double s) const { return {x / s, y / s, z / s}; }
    };
    using Parameters = skymodelr_prague_parameters_v1;
    using AvailableData = skymodelr_prague_available_v1;

    PragueSkyModel() : api(resolve()) {}
    PragueSkyModel(const std::string &filename, double visibility = 0) : PragueSkyModel() {
        initialize(filename, visibility);
    }
    ~PragueSkyModel() { if (handle) api->destroy(handle); }
    PragueSkyModel(const PragueSkyModel &) = delete;
    PragueSkyModel &operator=(const PragueSkyModel &) = delete;

    // Initialization/reinitialization requires exclusive access. On failure,
    // a previously initialized model remains usable.
    void initialize(const std::string &filename, double visibility = 0) {
        char error[512];
        auto *next = api->create(filename.c_str(), visibility, error, sizeof(error));
        if (!next) throw std::runtime_error(error);
        if (handle) api->destroy(handle);
        handle = next;
    }
    AvailableData getAvailableData() const {
        AvailableData result; char error[512];
        if (!api->available(handle, &result, error, sizeof(error))) throw std::runtime_error(error);
        return result;
    }
    size_t memoryUsage() const { return sizeof(*this) + api->memory_usage(handle); }
    Parameters computeParameters(Vector3 p, Vector3 w, double elevation, double azimuth,
                                  double visibility, double albedo) const {
        const double position[3] = {p.x, p.y, p.z}, direction[3] = {w.x, w.y, w.z};
        Parameters result; char error[512];
        if (!api->parameters(handle, position, direction, elevation, azimuth, visibility,
                              albedo, &result, error, sizeof(error))) throw std::runtime_error(error);
        return result;
    }
    void skyRadianceSpectrum(const Parameters &p, const double *w, size_t n, double *out) const {
        spectrum(SKYMODELR_PRAGUE_SKY, p, w, n, 0, true, out);
    }
    void transmittanceSpectrum(const Parameters &p, const double *w, size_t n,
                                double distance, double *out) const {
        spectrum(SKYMODELR_PRAGUE_TRANSMITTANCE, p, w, n, distance, true, out);
    }
    double skyRadiance(const Parameters &p, double wavelength) const {
        double result; skyRadianceSpectrum(p, &wavelength, 1, &result); return result;
    }
    double sunRadiance(const Parameters &p, double wavelength, bool attenuate = true) const {
        double result;
        spectrum(SKYMODELR_PRAGUE_SUN, p, &wavelength, 1, 0, attenuate, &result);
        return result;
    }
    double transmittance(const Parameters &p, double wavelength, double distance) const {
        double result; transmittanceSpectrum(p, &wavelength, 1, distance, &result); return result;
    }
    double polarisation(const Parameters &p, double wavelength) const {
        double result; spectrum(SKYMODELR_PRAGUE_POLARISATION, p, &wavelength, 1, 0, true, &result);
        return result;
    }
private:
    const skymodelr_prague_api_v1 *api;
    skymodelr_prague_handle_v1 *handle = nullptr;
    static const skymodelr_prague_api_v1 *resolve() {
        // No lazy symbol resolution in worker queries or destructors.
        auto get = reinterpret_cast<skymodelr_prague_get_api_v1_fn>(
            R_GetCCallable("skymodelr", "prague_get_api_v1"));
        const auto *result = get();
        if (!result || result->abi_version != 1 || result->struct_size < sizeof(*result))
            throw std::runtime_error("Incompatible skymodelr Prague native API; update skymodelr.");
        return result;
    }
    void spectrum(int kind, const Parameters &p, const double *w, size_t n,
                    double distance, bool attenuate, double *out) const {
        char error[512];
        if (!api->spectrum(handle, kind, &p, w, n, distance, attenuate ? 1 : 0,
                            out, error, sizeof(error))) throw std::runtime_error(error);
    }
};
}
#endif
