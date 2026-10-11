/* Prague Metal data ABI. GPL-3; coefficient/model provenance remains upstream.
 * No Objective-C or Metal framework types cross the package boundary. */
#ifndef SKYMODELR_PRAGUE_METAL_H
#define SKYMODELR_PRAGUE_METAL_H
#include "prague_api.h"
#ifdef __cplusplus
extern "C" {
#endif
#define SKYMODELR_PRAGUE_METAL_ABI_VERSION 1
/* All offsets index float elements in one immutable, packed buffer. */
typedef struct {
    uint32_t offset, count;
} skymodelr_prague_gpu_axis;
typedef struct {
    skymodelr_prague_gpu_axis visibility, albedo, altitude, elevation;
    skymodelr_prague_gpu_axis gamma, alpha, emphasis, trans_altitude, trans_visibility, sun;
    uint32_t radiance, trans_u, trans_v, rank, sun_offset, sun_stride;
    uint32_t zenith_offset, zenith_stride, emphasis_offset, config_stride;
    uint32_t trans_rank, a_dim, d_dim, channels;
    float channel_start, channel_width;
} skymodelr_prague_gpu_layout;
typedef struct skymodelr_prague_gpu_handle skymodelr_prague_gpu_handle;
typedef struct {
    float theta, gamma, shadow, zero, elevation, altitude, visibility, albedo;
} skymodelr_prague_gpu_parameters;
typedef struct {
    skymodelr_prague_gpu_parameters parameters;
    float wavelength, distance;
    uint32_t quantity, attenuate_sun;
} skymodelr_prague_gpu_query;
/* Three float4 records: position.xyz (metres), direction.xyz, then
 * elevation/azimuth (radians), visibility (km), albedo. Unused w lanes are zero. */
typedef struct { float position[4], direction[4], settings[4]; } skymodelr_prague_gpu_geometry;
typedef struct {
    uint32_t abi_version;
    size_t struct_size;
    /* Snapshot owns packed data independently of the CPU model. Its pointers
     * remain valid until destroy; create is an explicit allocation boundary. */
    skymodelr_prague_gpu_handle *(*create)(const skymodelr_prague_handle *, char *, size_t);
    void (*destroy)(skymodelr_prague_gpu_handle *);
    const skymodelr_prague_gpu_layout *(*layout)(const skymodelr_prague_gpu_handle *);
    const float *(*data)(const skymodelr_prague_gpu_handle *, size_t *count);
    /* Null-terminated MSL with device-callable functions and a standalone
     * prague_spectrum_batch kernel. Lifetime equals provider DLL lifetime. */
    const char *(*source)(void);
} skymodelr_prague_gpu_api;
typedef const skymodelr_prague_gpu_api *(*skymodelr_prague_gpu_get_api_fn)(void);
#ifdef __cplusplus
}
#endif
#endif
