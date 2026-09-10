/* skymodelr native Prague API. GPL-3. */
#ifndef SKYMODELR_PRAGUE_API_H
#define SKYMODELR_PRAGUE_API_H
#include <stddef.h>
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

/* Only this ABI is provided. Bump the version for incompatible layout or
 * signature changes; consumers check it and struct_size before using the table.
 * Resolve prague_get_api through R_GetCCallable("skymodelr", ...) on R's
 * main thread after loading the namespace. All table functions are native,
 * never call R, and contain exceptions within the provider DLL.
 * Handles own immutable coefficients; queries may share a handle across
 * threads. Do not destroy a handle until every query using it has finished.
 * Keep the provider namespace loaded for the lifetime of the table/handles.
 * Status functions return 1 on success, 0 on failure. Errors are copied into
 * the caller's optional buffer, with a terminator when error_capacity > 0.
 */
#define SKYMODELR_PRAGUE_ABI_VERSION 3
typedef struct skymodelr_prague_handle skymodelr_prague_handle;
typedef struct {
    double theta, gamma, shadow, zero, elevation, altitude, visibility, albedo;
} skymodelr_prague_parameters;
/* Dataset bounds retain upstream units: solar elevation is in degrees,
 * altitude in meters, visibility in kilometers, spectral channels in nm.
 * Query parameters above/below use radians for all angles. */
typedef struct {
    double albedoMin, albedoMax, altitudeMin, altitudeMax;
    double elevationMin, elevationMax, visibilityMin, visibilityMax;
    int32_t polarisation, channels;
    double channelStart, channelWidth;
} skymodelr_prague_available;
enum {
    SKYMODELR_PRAGUE_SKY = 0, SKYMODELR_PRAGUE_SUN = 1,
    SKYMODELR_PRAGUE_TRANSMITTANCE = 2, SKYMODELR_PRAGUE_POLARISATION = 3
};
typedef struct {
    uint32_t abi_version;
    size_t struct_size;
    /* Options are fixed at creation. cache_spectra and transmission_table must
     * each be 0 or 1. max_mib is nonnegative: zero skips table expansion and
     * +Inf removes the memory cap. Both accelerations are exact. A table that
     * exceeds the cap or cannot be allocated uses the compressed evaluator. */
    skymodelr_prague_handle *(*create)(const char *filename, double visibility,
                                      int cache_spectra, int transmission_table,
                                      double max_mib, char *error, size_t error_capacity);
    void (*destroy)(skymodelr_prague_handle *);
    int (*available)(const skymodelr_prague_handle *, skymodelr_prague_available *,
                      char *error, size_t error_capacity);
    size_t (*memory_usage)(const skymodelr_prague_handle *);
    /* Model coordinates: +Z is up, X/Y are tangent-plane meters, Earth center
     * is (0, 0, -6378000). Angles are radians, visibility kilometers. */
    int (*parameters)(const skymodelr_prague_handle *, const double position[3],
                       const double direction[3], double solar_elevation,
                       double solar_azimuth, double visibility, double albedo,
                       skymodelr_prague_parameters *, char *error, size_t error_capacity);
    /* Wavelengths are nm; output has count doubles. Distance is meters;
     * use DBL_MAX for infinity. attenuate_sun is 0 or 1 (Sun only).
     * Sky and transmittance evaluate a spectral batch with shared lookups. */
    int (*spectrum)(const skymodelr_prague_handle *, int quantity,
                     const skymodelr_prague_parameters *, const double *wavelengths,
                     size_t count, double distance, int attenuate_sun, double *output,
                     char *error, size_t error_capacity);
} skymodelr_prague_api;
typedef const skymodelr_prague_api *(*skymodelr_prague_get_api_fn)(void);

#ifdef __cplusplus
}
#endif
#endif
