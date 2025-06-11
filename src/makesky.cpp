#include <Rcpp.h>
#include <RcppThread.h>

extern "C" {
#include "ArHosekSkyModel.h"
}

// OpenEXR
#include <OpenEXR/ImfRgbaFile.h>
#include <OpenEXR/ImfRgba.h>
#include <OpenEXR/ImfChannelList.h>

#include <Imath/ImathVec.h>
#include <Imath/ImathFun.h>
using namespace Imf;
using namespace Imath;

// ───────────────────────────── Helper maths ──────────────────────────────
static inline float clamp_float(float x, float a, float b) {
  return (x < a) ? a : (x > b) ? b : x;
}

static inline Vec3<float> latLongSquareToSphere(float u, float v) {
  float theta = static_cast<float>(M_PI) * v;       // 0…π
  float phi   = static_cast<float>(2.0 * M_PI) * u; // 0…2π
  return { std::sin(theta) * std::cos(phi),
           std::sin(theta) * std::sin(phi),
           std::cos(theta) };
}

static inline Vec3<float> equalAreaSquareToSphere(float u, float v) {
  float a = 2.f * u - 1.f;
  float b = 2.f * v - 1.f;
  float r, phi;
  if (std::fabs(a) > std::fabs(b)) {
    r   = a;
    phi = static_cast<float>(M_PI_4) * (b / r);
  } else {
    r   = b;
    phi = static_cast<float>(M_PI_2) - static_cast<float>(M_PI_4) * (a / r);
  }
  float theta = std::acos(r);
  float sinT  = std::sin(theta);
  return { std::cos(phi) * sinT,
           std::sin(phi) * sinT,
           std::cos(theta) };
}

// ───────────────────────────── C++ entrypoint ────────────────────────────
// [[Rcpp::export]]
void makesky_rcpp(std::string  outfile,
                  double       albedo    = 0.5,
                  double       turbidity = 3.0,
                  double       elevation = 10.0,
                  unsigned int resolution = 2048,
                  unsigned int numbercores = 1,
                  bool         square_projection = false) {

  /* ── user‑input guardrails ─────────────────────────────────────── */
  if (albedo    < 0.0 || albedo    > 1.0)  Rcpp::stop("albedo in [0,1]");
  if (turbidity < 1.7 || turbidity > 10.0) Rcpp::stop("turbidity in [1.7,10]");
  if (elevation < 0.0 || elevation > 90.0) Rcpp::stop("elevation in [0,90]");
  if (resolution == 0)                     Rcpp::stop("resolution ≥ 1");

  float elevRad = static_cast<float>(elevation) * static_cast<float>(M_PI) / 180.0f;


  constexpr int num_channels = 9;
  // Three wavelengths around red, three around green, and three around blue.
  constexpr double lambda[num_channels] = {630, 680, 710, 500, 530, 560, 460, 480, 490};

  ArHosekSkyModelState *skymodel_state[num_channels];
  for (int i = 0; i < num_channels; ++i) {
    skymodel_state[i] =
      arhosekskymodelstate_alloc_init(elevRad, turbidity, albedo);
    if (!skymodel_state[i]) Rcpp::stop("Hosek–Wilkie initialisation failed");

  }

  Imath::Vec3<float> sunDir(0.f, std::sin(elevRad), std::cos(elevRad));
  int nTheta = resolution, nPhi = 2 * nTheta;
  std::vector<float> img(3 * nTheta * nPhi, 0.f);

  RcppThread::parallelFor(0u, resolution,                                    // NEW begin/end
                          [&](size_t t) {
                          float theta = (float(t) + 0.5f) / static_cast<float>(nTheta) * static_cast<float>(M_PI);
                          if (theta > static_cast<float>(M_PI) / 2.f) return;
                            for (int p = 0; p < nPhi; ++p) {
                              float phi = (float(p) + 0.5f) / nPhi * 2.f * static_cast<float>(M_PI);

                              // Vector corresponding to the direction for this pixel.
                              Imath::Vec3<float> v(std::cos(phi) * std::sin(theta),
                                                   std::cos(theta),
                                                   std::sin(phi) * std::sin(theta));
                              // Compute the angle between the pixel's direction and the sun
                              // direction.
                              //
                              float gamma = std::acos(clamp_float(v.dot(sunDir), -1.f, 1.f));

                              for (int c = 0; c < num_channels; ++c) {
                                float val = arhosekskymodel_solar_radiance(
                                  skymodel_state[c], theta, gamma, lambda[c]);
                                // For each of red, green, and blue, average the three
                                // values for the three wavelengths for the color.
                                // TODO: do a better spectral->RGB conversion.
                                img[3 * (t * nPhi + p) + c / 3] += val / 3.f;
                              }
                            }
                          },
                          numbercores);

  /* ── pack + write EXR ──────────────────────────────────────────── */
  const int width  = nPhi;          // 2×resolution
  const int height = nTheta;        // resolution
  const int nPix   = width * height;

  std::vector<Imf::Rgba> px(nPix);
  for (int p = 0; p < nPix; ++p) {
      px[p] = Imf::Rgba(img[3*p], img[3*p+1], img[3*p+2]);
  }

  try {
    Imf::RgbaOutputFile out(outfile.c_str(),
                            width, height,
                            Imf::WRITE_RGB);

    out.setFrameBuffer(px.data(), 1, width);
    out.writePixels(height);
  } catch (const std::exception &e) {
    for(int i = 0; i < num_channels; i++) {
      arhosekskymodelstate_free(skymodel_state[i]);
    }
    Rcpp::stop("OpenEXR write failed: %s", e.what());
  }

  for(int i = 0; i < num_channels; i++) {
    arhosekskymodelstate_free(skymodel_state[i]);
  }
}

