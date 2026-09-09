# Native Prague atmosphere interface

Skymodelr 0.6.1 supplies a versioned native interface for packages that
query the Prague model from C or C++. The implementation is compiled once, in
skymodelr. Consumers include the installed headers and obtain a function table
through R's registered C-callable mechanism. 

This follows [R's documented native package linking mechanism](https://stat.ethz.ch/CRAN/doc/manuals/r-release/R-exts.html#Linking-to-native-routines-in-other-packages).
The implementation remains in `src/PragueSkyModel/`.

## Package setup

A consuming package needs both a runtime import and header dependency in its
`DESCRIPTION`:

```
Imports: skymodelr (>= 0.6.1)
LinkingTo: skymodelr (>= 0.6.1)
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
resolve `R_GetCCallable("skymodelr", "prague_get_api_v1")` and cast the result to
`skymodelr_prague_get_api_v1_fn`. Call the getter and check `abi_version == 1`
and `struct_size >= sizeof(skymodelr_prague_api_v1)` before use.

Create a handle with `create()`, issue `parameters()` and `spectrum()` queries,
and release it with `destroy()`. `spectrum()` accepts an array of wavelengths
and a matching output buffer for all four quantities. Status-returning functions
return one on success and zero on error; `create()` returns null on error.
An optional error buffer receives a terminated message on failure when its
capacity is positive. Output buffers should only be read on success. Use
`DBL_MAX` for an infinite transmission distance (`+Inf` is also accepted).

The version-one table, records, function signatures and callable name are
frozen. An incompatible change must publish a new versioned callable while
retaining version one for existing consumers. The registration hook uses
`Rcpp::init` so `Rcpp::compileAttributes()` preserves it. There is no need to
export implementation symbols from the DLL.

The installed-header integration test compiles a separate consumer, exercises
parallel queries and checks error containment. 
