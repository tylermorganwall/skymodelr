#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include "prague_metal_backend.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {

// One thread evaluates one row, including geometry and both radiance terms.
// Keep the host and shader records scalar-only to avoid float3 padding rules.
struct RadianceRow {
  float phi, theta, wavelength, elevation, albedo, altitude, visibility, azimuth;
};
static_assert(sizeof(RadianceRow) == 8 * sizeof(float), "Metal row layout");

const char *radiance_source =
    "#include <metal_stdlib>\nusing namespace metal;\n"
#include "../inst/include/skymodelr/prague_metal.inc"
    R"MSL(
struct RadianceRow {
  float phi, theta, wavelength, elevation, albedo, altitude, visibility, azimuth;
};
kernel void skymodelr_radiance_batch(
    device const skymodelr_gpu::Layout &model [[buffer(0)]],
    device const float *data [[buffer(1)]],
    device const RadianceRow *rows [[buffer(2)]],
    device float *output [[buffer(3)]],
    constant uint &count [[buffer(4)]],
    constant uint &mode [[buffer(5)]], uint i [[thread_position_in_grid]]) {
  if (i >= count) return;
  auto q = rows[i];
  float radians = skymodelr_gpu::pi / 180.f;
  float phi = q.phi * radians, theta = q.theta * radians;
  // Match the CPU's protection against the exact solar zenith singularity.
  float elevation = min(q.elevation, 89.9999f) * radians;
  float3 direction(cos(phi) * cos(theta), sin(phi) * cos(theta), sin(theta));
  // cos(float(pi/2)) leaves a horizontal component at the poles. The
  // transmission fit's fourth-root distance coordinate magnifies that error.
  if (abs(q.theta) == 90.f) direction = float3(0, 0, sign(q.theta));
  auto p = skymodelr_gpu::parameters(float3(0, 0, q.altitude), direction,
      elevation, q.azimuth * radians, q.visibility, q.albedo);
  float result = 0;
  if (mode & 1) result += skymodelr_gpu::sky(p, q.wavelength, model, data);
  if (mode & 2) result += skymodelr_gpu::sun(p, q.wavelength, true, model, data);
  output[i] = isfinite(result) ? result : 0.f;
}
)MSL";

// Called only from R's main thread. ARC owns all Metal resources, including
// partially constructed state after an error. Cache one dataset at a time.
struct MetalRadianceState {
  id<MTLDevice> device;
  id<MTLCommandQueue> queue;
  id<MTLComputePipelineState> pipeline;
  id<MTLBuffer> layout_buffer, coefficients;
  skymodelr_prague_gpu_layout layout;
  std::string dataset;
};

MetalRadianceState state;

std::runtime_error metal_error(const char *operation, NSError *error) {
  const char *detail = error ? [[error localizedDescription] UTF8String] : nullptr;
  return std::runtime_error(std::string(operation) +
                            (detail ? std::string(": ") + detail : ""));
}

void prepare_metal(const PragueSkyModel &model, const std::string &dataset) {
  if (!state.device) state.device = MTLCreateSystemDefaultDevice();
  if (!state.device) {
    throw std::runtime_error(
        "backend = \"gpu\" requires macOS with an available Metal device. "
        "Use backend = \"cpu\" on this system.");
  }
  if (!state.queue) {
    state.queue = [state.device newCommandQueue];
    if (!state.queue) throw metal_error("Could not create Metal command queue", nil);
  }
  if (!state.pipeline) {
    MTLCompileOptions *options = [MTLCompileOptions new];
#if __MAC_OS_X_VERSION_MAX_ALLOWED >= 150000
    if (@available(macOS 15.0, *)) {
      options.mathMode = MTLMathModeSafe;
      options.mathFloatingPointFunctions = MTLMathFloatingPointFunctionsPrecise;
    } else
#endif
    {
      options.fastMathEnabled = NO;
    }
    NSError *error = nil;
    id<MTLLibrary> library = [state.device
        newLibraryWithSource:[NSString stringWithUTF8String:radiance_source]
        options:options error:&error];
    if (!library) throw metal_error("Could not compile Prague Metal shaders", error);
    id<MTLFunction> function = [library newFunctionWithName:@"skymodelr_radiance_batch"];
    state.pipeline = [state.device newComputePipelineStateWithFunction:function error:&error];
    if (!state.pipeline) throw metal_error("Could not create Prague Metal pipeline", error);
  }
  if (state.dataset == dataset) return;

  // Release the old GPU copy before exporting a potentially large new dataset.
  // Publish the cache key last, so an allocation failure remains retryable.
  state.dataset.clear();
  state.coefficients = nil;
  state.layout_buffer = nil;
  std::vector<float> packed;
  model.exportMetal(state.layout, packed);
  const size_t bytes = packed.size() * sizeof(float);
  if (bytes > state.device.maxBufferLength) {
    throw std::runtime_error("Prague coefficients exceed this Metal device's buffer limit. "
                             "Use backend = \"cpu\" for this dataset.");
  }
  state.coefficients = [state.device newBufferWithBytes:packed.data()
      length:bytes options:MTLResourceStorageModeShared];
  state.layout_buffer = [state.device newBufferWithBytes:&state.layout
      length:sizeof(state.layout) options:MTLResourceStorageModeShared];
  if (!state.coefficients || !state.layout_buffer) {
    state.coefficients = nil;
    state.layout_buffer = nil;
    throw metal_error("Could not allocate Prague Metal coefficients", nil);
  }
  state.dataset = dataset;
}

} // namespace

bool prague_metal_available() {
  @autoreleasepool {
    return MTLCreateSystemDefaultDevice() != nil;
  }
}

void prague_metal_radiance(const PragueSkyModel &model,
                           const std::string &dataset,
                           const PragueRadianceColumns &columns,
                           size_t count, bool atmosphere, bool sun,
                           double *output, void (*check_interrupt)()) {
  if (count == 0) return;
  @autoreleasepool {
    prepare_metal(model, dataset);

    // Bound temporary memory and GPU command duration independently of grid
    // size. Reuse these buffers and allow R interrupts between completed batches.
    const size_t capacity = std::min(count, size_t(65536));
    id<MTLBuffer> input = [state.device newBufferWithLength:capacity * sizeof(RadianceRow)
        options:MTLResourceStorageModeShared];
    id<MTLBuffer> result = [state.device newBufferWithLength:capacity * sizeof(float)
        options:MTLResourceStorageModeShared];
    if (!input || !result) throw metal_error("Could not allocate Metal batch buffers", nil);
    auto *rows = static_cast<RadianceRow *>([input contents]);
    const uint32_t mode = (atmosphere ? 1u : 0u) | (sun ? 2u : 0u);
    const NSUInteger width = state.pipeline.threadExecutionWidth;

    for (size_t start = 0; start < count; start += capacity) {
      check_interrupt();
      const uint32_t batch_count = std::min(capacity, count - start);
      for (size_t j = 0; j < batch_count; ++j) {
        const size_t i = start + j;
        // Preserve the CPU spectral bin when a double just below a channel
        // boundary would otherwise round up in single precision.
        const double wavelength = columns.wavelength[i];
        const auto &m = state.layout;
        const unsigned channel = (wavelength - m.channel_start) / m.channel_width;
        const float upper = m.channel_start + (channel + 1) * m.channel_width;
        const float gpu_wavelength = std::min(float(wavelength), std::nextafter(upper, -INFINITY));
        rows[j] = {float(columns.phi[i]), float(columns.theta[i]), gpu_wavelength,
                   float(columns.elevation[i]), float(columns.albedo[i]),
                   float(columns.altitude[i]), float(columns.visibility[i]),
                   float(columns.azimuth[i])};
      }
      @autoreleasepool {
        id<MTLCommandBuffer> command = [state.queue commandBuffer];
        id<MTLComputeCommandEncoder> encoder = [command computeCommandEncoder];
        if (!command || !encoder) throw metal_error("Could not encode Metal batch", nil);
        [encoder setComputePipelineState:state.pipeline];
        [encoder setBuffer:state.layout_buffer offset:0 atIndex:0];
        [encoder setBuffer:state.coefficients offset:0 atIndex:1];
        [encoder setBuffer:input offset:0 atIndex:2];
        [encoder setBuffer:result offset:0 atIndex:3];
        [encoder setBytes:&batch_count length:sizeof(batch_count) atIndex:4];
        [encoder setBytes:&mode length:sizeof(mode) atIndex:5];
        [encoder dispatchThreadgroups:MTLSizeMake((batch_count + width - 1) / width, 1, 1)
            threadsPerThreadgroup:MTLSizeMake(width, 1, 1)];
        [encoder endEncoding];
        [command commit];
        [command waitUntilCompleted];
        if (command.status != MTLCommandBufferStatusCompleted) {
          throw metal_error("Prague Metal batch failed", command.error);
        }
        const auto *values = static_cast<const float *>([result contents]);
        std::copy(values, values + batch_count, output + start);
      }
    }
  }
}
