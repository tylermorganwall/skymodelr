# Native Prague evaluation in Metal

R callers can use `calculate_sky_radiance(..., backend = "gpu")` on macOS.
Pass columns from `expand.grid()` to evaluate all requested combinations; each
input row produces one spectral radiance value in the same order. This path
calls the device functions below for geometry, sky, and attenuated sunlight,
caches the compiled pipeline and one dataset's GPU coefficients, and processes
large inputs in bounded batches. All three `render_mode` values are supported.
The default `backend = "cpu"` retains CPU evaluation. Explicit GPU requests
fail if Metal is unavailable; results use single precision and may differ
slightly from CPU values, particularly at the solar limb.

A renderer can evaluate the Prague spectral model inside its GPU work queues,
without a CPU callback for every sky or transmittance query. This API exports
compressed model coefficients and reusable Metal Shading Language functions;
it does not bake the sky into an RGB image or angular lookup table. The CPU
interface and non-Apple builds do not require a Metal SDK.

On R's main thread, resolve `prague_gpu_get_api` with `R_GetCCallable` after
loading skymodelr. Check `abi_version` and `struct_size` from
`skymodelr/prague_metal.h`. `create` takes an initialized CPU model handle and
returns an independently owned packed snapshot. `layout` describes offsets in
the float array returned by `data`. The C++ wrapper's `metalSnapshot()` caches
one immutable snapshot per model; reinitialization invalidates that cache.
Keep skymodelr loaded until every snapshot is destroyed. Never destroy a
snapshot while a consumer uses its borrowed pointers.

Copy the layout and coefficients into consumer-owned Metal buffers once, then
reuse them. The `source` function supplies null-terminated MSL to prepend to a
consumer's shader source after `metal_stdlib` and `using namespace metal`.
The same source is installed as the C++ string include
`skymodelr/prague_metal.inc`. Compile with safe math (fast math disabled).

Device functions in `skymodelr_gpu` are `parameters`, `sky`, `sun`, and
`transmission`. Parameters and wavelengths use the same units as the CPU API:
metres in the local Z-up Earth frame, radians, visibility in kilometres, and
wavelengths in nanometres. Queries must satisfy the CPU API's finite-input and
above-Earth preconditions. The GPU does not throw exceptions or clamp arbitrary
invalid inputs. Polarisation is not currently exported.

Two standalone kernels support bulk consumers:

- `prague_parameters_batch`: buffers 0/1/2 are geometry input, parameter output,
  and a uint32 query count. Public C record types specify their layouts.
- `prague_spectrum_batch`: buffers 0/1/2/3/4 are layout, coefficients, query
  records, float output, and uint32 query count. Quantities are sky=0, sun=1,
  transmission=2; other values are not supported. Infinity is accepted for
  transmittance distance. Solar queries optionally include attenuation.

The evaluator retains spectral coefficients and interpolation order but uses
float arithmetic; it is not bit-identical to the CPU's double implementation.
Earth-scale geometry uses rationalized heights/roots and atan2 angular distances
to retain small atmospheric distances. GPU ray transport, spectral-to-RGB
quadrature, visibility, and finite-segment source accounting belong to the
consumer, not this model API. Standalone batch readback is useful for bulk tools;
renderers should call the functions directly within GPU queues to avoid handoffs.

For repeated wavelengths at one point, `prepare_sky` / `sky_prepared` and
`prepare_transmission` / `transmission_prepared` share geometry and interpolation
coordinates. Consumers may distribute wavelengths across SIMD lanes and reduce
their own RGB quadrature weights; skymodelr imposes no RGB color space or fixed
wavelength count. Weighted corner accumulation can change float rounding order.
