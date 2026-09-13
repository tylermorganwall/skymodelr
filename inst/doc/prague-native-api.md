# Native Prague atmosphere interface

Skymodelr 0.6.2 supplies one native interface for packages that
query the Prague model from C or C++. The implementation is compiled once, in
skymodelr. Consumers include the installed headers and obtain a function table
through R's registered C-callable mechanism. 

This follows [R's documented native package linking mechanism](https://stat.ethz.ch/CRAN/doc/manuals/r-release/R-exts.html#Linking-to-native-routines-in-other-packages).
The implementation remains in `src/PragueSkyModel/`.

## Package setup

A consuming package needs both a runtime import and header dependency in its
`DESCRIPTION`:

```
Imports: skymodelr (>= 0.6.3)
LinkingTo: skymodelr (>= 0.6.3)
```

Keep any other dependencies already present in those fields. Import a public R
function that the consumer uses, for example
`importFrom(skymodelr,get_prague_sky_metadata)` in `NAMESPACE`. This loads the
provider namespace and DLL before native queries begin. `LinkingTo` alone only
makes headers available at compile time. Skymodelr requires R >= 4.3.

Install the updated skymodelr package before building the consumer. The Prague
coefficients are separate data, installed by `download_sky_data()`. Consumers
can use `get_prague_sky_metadata()` to resolve the coefficient file and obtain
scene metadata; native initialization itself never downloads data.

## C++ example

The installed `skymodelr/prague.hpp` header provides a small owning wrapper. It
contains no model implementation or coefficients. This Rcpp example can be
compiled with `Rcpp::sourceCpp()` after `library(skymodelr)`:

```cpp
// [[Rcpp::depends(skymodelr)]]
#include <Rcpp.h>
#include <skymodelr/prague.hpp>

// [[Rcpp::export]]
Rcpp::NumericVector zenith_spectrum(std::string filename) {
    // Resolve R's callable on the main thread and load coefficients here.
    skymodelr::PragueSkyModel model(filename, 50);
    auto parameters = model.computeParameters(
        {0, 0, 0}, {0, 0, 1}, 0.3, 0.7, 50, 0.5);
    const double wavelengths[] = {440, 550, 680};
    double radiance[3];
    model.skyRadianceSpectrum(parameters, wavelengths, 3, radiance);
    return Rcpp::NumericVector(radiance, radiance + 3);
}
```

The native coordinate system has +Z up, with horizontal X/Y coordinates in
meters and the Earth's center at `(0, 0, -6378000)`. Angles are radians,
visibility is kilometers, altitude and travel distance are meters, wavelengths
are nanometers, and ground albedo is between zero and one. These angular units
apply to queries; `getAvailableData()` retains the upstream coefficient format
and reports solar-elevation bounds in **degrees**. `computeParameters()`
accounts for the spherical Earth; callers map their own scene coordinates into
this system. `getAvailableData()` reports the loaded coefficient domain.
Queries return the model's spectral radiance, dimensionless transmittance, or
polarisation; RGB conversion and application-specific calibration remain the
consumer's responsibility.

The wrapper exposes scalar sky, Sun, transmission and polarisation queries,
and batched sky and transmission queries. Batched evaluation shares lookups
across wavelengths. `sunRadiance(..., false)` returns intrinsic solar radiance
without atmospheric attenuation. `memoryUsage()` reports model storage.
Initialization can load one visibility subset or all visibilities (zero).
A failed reinitialization leaves an existing model usable.

The C++ wrapper enables two exact accelerations. `initialize(filename, visibility, cache_spectra, transmission_table)`
and the corresponding constructor accept two optional booleans, both defaulting
to `true`. The spectral cache stores 16 exact queries per thread for batches of
up to 16 wavelengths. Model generations invalidate entries on reinitialization.
Other queries use the ordinary evaluator.

The transmission table precomputes the compressed rank reconstruction at the
dataset's native knots in double precision. It preserves the native interpolation,
clamps, and final square; it introduces no new model approximation. For a single
visibility it expands only the needed one or two slices. At 50 km visibility on
the full visible dataset this adds about 128 MiB per model, shared by threads.
Queries outside those slices use the compressed evaluator. Expansion is skipped
if it would exceed the memory budget or allocation fails. Pass `false` for the
table option to avoid this memory cost. Options are fixed before workers start.

A fifth constructor/initialization argument,
`transmission_table_max_mib`, controls the extra table memory budget per model.
It defaults to 512 MiB, accepts nonnegative fractional values, and accepts zero
to skip expansion or positive infinity to remove the cap. The entire required
table must fit; the fallback retains exact compressed evaluation. For example,
`PragueSkyModel model(filename, 50, true, true, 1024)` permits up to 1 GiB. A
55-channel extended-spectrum table needs about 641.5 MiB for two visibility
slices, so that table fits with this limit but not with the default.

## Threading, lifetime, and errors

Construct wrappers on R's main thread, after loading skymodelr: construction
calls `R_GetCCallable()`. Complete coefficient initialization before starting
workers. All subsequent const queries use native function pointers, call no R
API, and can share the same immutable model across threads. Initialization and
destruction require exclusive access. Wait for workers before destroying their
model, and keep skymodelr's namespace and DLL loaded while any model or function
table remains in use.

Allocation and destruction happen inside skymodelr's DLL. C++ exceptions are
caught there and converted to status codes with caller-owned error text. The
C++ wrapper converts failures into local `std::runtime_error` exceptions;
consumers must handle worker exceptions through their own thread/error system.
No C++ exceptions, STL containers, or ownership of provider allocations cross
the native boundary.

## C ABI and maintenance

`skymodelr/prague_api.h` is a C-compatible header describing the opaque handle,
plain parameter records and immutable function table. On the R main thread,
resolve `R_GetCCallable("skymodelr", "prague_get_api")` and cast the result to
`skymodelr_prague_get_api_fn`. Call the getter and check
`abi_version == SKYMODELR_PRAGUE_ABI_VERSION` (currently 3) and
`struct_size >= sizeof(skymodelr_prague_api)` before using any function pointer.

Create a handle with
`create(filename, visibility, cache_spectra, transmission_table, max_mib, error, error_capacity)`,
issue `parameters()` and `spectrum()` queries,
and release it with `destroy()`. `spectrum()` accepts an array of wavelengths
and a matching output buffer for all four quantities. Status-returning functions
return one on success and zero on error; `create()` returns null on error.
An optional error buffer receives a terminated message on failure when its
capacity is positive. Output buffers should only be read on success. Use
`DBL_MAX` for an infinite transmission distance (`+Inf` is also accepted).

There is one flat function table and one creation function. Both cache flags
must be zero or one, and the table budget is in MiB. Creation options are fixed
before threads can query the model. To compare against the original evaluator,
create another handle through the same API with both cache flags set to zero.

This development interface replaces the earlier experimental tables and getter
names. Cleanly rebuild skymodelr first and then all consumers, and restart R to
unload older DLLs. No compatibility tables or aliases are retained. The ABI
identifier and structure-size checks detect mismatched SDKs and providers;
incompatible future changes must increment the ABI identifier. The registration
hook uses `Rcpp::init` so `Rcpp::compileAttributes()` preserves it. There is no
need to export implementation symbols from the DLL.

The installed-header integration test compiles a separate consumer, exercises
parallel queries and error containment, and compares accelerated and unaccelerated
handles through this one interface. Memory-budget tests cover allocation thresholds,
zero, infinity, and the extended-spectrum dataset.
