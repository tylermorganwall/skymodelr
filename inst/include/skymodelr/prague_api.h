/* skymodelr native Prague API, version 1. GPL-3. */
#ifndef SKYMODELR_PRAGUE_API_H
#define SKYMODELR_PRAGUE_API_H
#include <stddef.h>
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

/* ABI v1 is frozen. New incompatible interfaces get a new callable name.
 * Resolve prague_get_api_v1 through R_GetCCallable("skymodelr", ...) on R's
 * main thread after loading the namespace. All table functions are native,
 * never call R, and contain exceptions within the provider DLL.
 * Handles own immutable coefficients; queries may share a handle across
 * threads. Do not destroy a handle until every query using it has finished.
 * Keep the provider namespace loaded for the lifetime of the table/handles.
 * Status functions return 1 on success, 0 on failure. Errors are copied into
 * the caller's optional buffer, with a terminator when error_capacity > 0.
 */
typedef struct skymodelr_prague_handle_v1 skymodelr_prague_handle_v1;
typedef struct {
    double theta, gamma, shadow, zero, elevation, altitude, visibility, albedo;
} skymodelr_prague_parameters_v1;
/* Dataset bounds retain upstream units: solar elevation is in degrees,
 * altitude in meters, visibility in kilometers, spectral channels in nm.
 * Query parameters above/below use radians for all angles. */
typedef struct {
    double albedoMin, albedoMax, altitudeMin, altitudeMax;
    double elevationMin, elevationMax, visibilityMin, visibilityMax;
    int32_t polarisation, channels;
    double channelStart, channelWidth;
} skymodelr_prague_available_v1;
enum {
    SKYMODELR_PRAGUE_SKY = 0, SKYMODELR_PRAGUE_SUN = 1,
    SKYMODELR_PRAGUE_TRANSMITTANCE = 2, SKYMODELR_PRAGUE_POLARISATION = 3
};
typedef struct {
    uint32_t abi_version;
    size_t struct_size;
    skymodelr_prague_handle_v1 *(*create)(const char *filename, double visibility,
                                         char *error, size_t error_capacity);
    void (*destroy)(skymodelr_prague_handle_v1 *);
    int (*available)(const skymodelr_prague_handle_v1 *, skymodelr_prague_available_v1 *,
                      char *error, size_t error_capacity);
    size_t (*memory_usage)(const skymodelr_prague_handle_v1 *);
    /* Model coordinates: +Z is up, X/Y are tangent-plane meters, Earth center
     * is (0, 0, -6378000). Angles are radians, visibility kilometers. */
    int (*parameters)(const skymodelr_prague_handle_v1 *, const double position[3],
                       const double direction[3], double solar_elevation,
                       double solar_azimuth, double visibility, double albedo,
                       skymodelr_prague_parameters_v1 *, char *error, size_t error_capacity);
    /* Wavelengths are nm; output has count doubles. Distance is meters;
     * use DBL_MAX for infinity. attenuate_sun is 0 or 1 (Sun only).
     * Sky and transmittance evaluate a spectral batch with shared lookups. */
    int (*spectrum)(const skymodelr_prague_handle_v1 *, int quantity,
                     const skymodelr_prague_parameters_v1 *, const double *wavelengths,
                     size_t count, double distance, int attenuate_sun, double *output,
                     char *error, size_t error_capacity);
} skymodelr_prague_api_v1;
typedef const skymodelr_prague_api_v1 *(*skymodelr_prague_get_api_v1_fn)(void);
#ifdef __cplusplus
}
#endif
#endif
