#include <Rcpp.h>
#include <RcppThread.h>

extern "C" {
#include "ArHosekSkyModel.h"
}

#include "PragueSkyModel.h"

// OpenEXR
#include <OpenEXR/ImfRgbaFile.h>
#include <OpenEXR/ImfRgba.h>
#include <OpenEXR/ImfChannelList.h>

#include <Imath/ImathVec.h>
#include <Imath/ImathFun.h>
#include <Imath/ImathBox.h>

#include <cmath>
#include <vector>
#include <string>

using namespace Imf;
using namespace Imath;

static inline float clamp_float(float x, float a, float b)
{ return (x < a) ? a : (x > b) ? b : x; }

// Forward conversion to Prague’s tiny Vector3 wrapper
static inline PragueSkyModel::Vector3 toPrague(const Vec3<float>& v)
  { return { double(v.x), double(v.y), double(v.z) }; }

// Static instance so the .dat stays in RAM across calls
static PragueSkyModel  prague_model;
static std::string     prague_loaded_file;

// [[Rcpp::export]]
void makesky_rcpp(std::string  outfile,
                  double       albedo            = 0.5,
                  double       turbidity         = 3.0,
                  double       elevation         = 10.0,
				          double       azimuth_deg       = 90,
                  unsigned int resolution        = 2048,
                  unsigned int numbercores       = 1,
                  bool         square_projection = false,
                  std::string  model             = "hosek",  // hosek | prague
                  std::string  prg_dataset       = "",
				          double       altitude          = 0.0,
                  double       visibility        = 50.0)     // Prague only (km)
{
  if (albedo < 0.0 || albedo > 1.0) {
    Rcpp::stop("albedo must be in [0,1]");
  }
  if (resolution == 0) {
    Rcpp::stop("resolution must be ≥ 1");
  }

  float elev_rad = static_cast<float>(elevation * M_PI / 180.0);
  float az_rad = static_cast<float>(azimuth_deg * M_PI / 180.0);

  // Prague model initialisation
  if (model == "prague") {
    if (prg_dataset.empty()) {
      Rcpp::stop("prg_dataset must be supplied when model=\"prague\"");
    }

    if (prague_loaded_file != prg_dataset) {
      prague_model.initialize(prg_dataset);
      prague_loaded_file = prg_dataset;
    }
    auto avail = prague_model.getAvailableData();
    if (albedo    < avail.albedoMin    || albedo    > avail.albedoMax ||
        visibility < avail.visibilityMin|| visibility > avail.visibilityMax) {
      Rcpp::stop("albedo or visibility outside dataset range");
    }
  } else if (model != "hosek") {
    Rcpp::stop("unknown model \"%s\"", model.c_str());
  }

  // Spectral sampling and Hosek state
  constexpr int    N_LAMBDA = 9;
  constexpr double lambda_nm[N_LAMBDA] =
    { 630,680,710, 500,530,560, 460,480,490 };

  ArHosekSkyModelState* hosek[N_LAMBDA] = { nullptr };
  if (model == "hosek") {
    if (elevation < 0.0 || elevation > 90.0) {
      Rcpp::stop("elevation must be in [0,90]");
    }
    if (turbidity < 1.7 || turbidity > 10.0) {
      Rcpp::stop("turbidity must be in [1.7,10]");
    }
    for (int i = 0; i < N_LAMBDA; ++i) {
      hosek[i] = arhosekskymodelstate_alloc_init(
        elev_rad, turbidity, albedo);
      if (!hosek[i]) {
        Rcpp::stop("Hosek–Wilkie initialisation failed");
      }
    }
  } else {
    if (elevation < -4.2|| elevation > 90.0) {
      Rcpp::stop("elevation must be in [-4.2, 90]");
    }
  }

  // Initialize image buffer
  const int nTheta = resolution;          // rows  (latitude)
  const int nPhi   = 2 * resolution;      // cols  (longitude)
  std::vector<float> img(3 * nTheta * nPhi, 0.0f);

  //   Vec3<float> sun_dir(0.0f, std::sin(elev_rad), std::cos(elev_rad));
  Vec3<float> sun_dir(std::cos(az_rad) * std::cos(elev_rad),  // +X (South)
                      std::sin(elev_rad),                     // +Y (Up)
                      std::sin(az_rad) * std::cos(elev_rad)); // +Z (West)
  // Render loop (parallel over rows)
  RcppThread::parallelFor(
    0u, size_t(nTheta),
    [&](size_t t) {
      float theta = (static_cast<float>(t) + 0.5f) /
        static_cast<float>(nTheta) * static_cast<float>(M_PI);

      if (theta > static_cast<float>(M_PI_2)) {
        return;
      }

      for (int p = 0; p < nPhi; ++p) {
        float phi = (static_cast<float>(p) + 0.5f) / nPhi *
          2.0f * static_cast<float>(M_PI);

        Vec3<float> v(std::cos(phi) * std::sin(theta),
                      std::cos(theta),
                      std::sin(phi) * std::sin(theta));

        if (model == "hosek") {
          float gamma = std::acos(
            clamp_float(v.dot(sun_dir), -1.0f, 1.0f));
          for (int c = 0; c < N_LAMBDA; ++c) {
            float L = arhosekskymodel_solar_radiance(
              hosek[c], theta, gamma, lambda_nm[c]);
            img[3 * (t * nPhi + p) + c / 3] += L / 3.0f;
          }
        } else {
          Vec3<float> v_zup(std::cos(phi) * std::sin(theta),
                            std::sin(phi) * std::sin(theta),
                            std::cos(theta));
          auto P = prague_model.computeParameters(
            /*viewPoint*/   { 0.0, 0.0, altitude },
            /*viewDir  */   toPrague(v_zup),
            /*sunElev  */   elev_rad,
            /*sunAzim  */   az_rad,
            /*visibility*/  visibility,
            /*albedo   */   albedo);
            for (int c = 0; c < N_LAMBDA; ++c) {
              double L = prague_model.skyRadiance(P, lambda_nm[c]) +
                prague_model.sunRadiance(P, lambda_nm[c]);
              img[3 * (t * nPhi + p) + c / 3] += static_cast<float>(L / 3.0);
            }
        }
      }
    },
    numbercores);

  // Pack and write EXR
  const int width  = nPhi;
  const int height = nTheta;
  const int nPix   = width * height;

  std::vector<Rgba> px(nPix);
  for (int p = 0; p < nPix; ++p) {
    px[p] = Rgba(img[3*p], img[3*p+1], img[3*p+2]);
  }

  try {
    RgbaOutputFile out(outfile.c_str(), width, height, WRITE_RGB);
    out.setFrameBuffer(px.data(), 1, width);
    out.writePixels(height);
  } catch (const std::exception& e) {
    if (model == "hosek")
      for (int i = 0; i < N_LAMBDA; ++i)
        arhosekskymodelstate_free(hosek[i]);
    Rcpp::stop("OpenEXR write failed: %s", e.what());
  }


  if (model == "hosek") {
    for (int i = 0; i < N_LAMBDA; ++i) {
      arhosekskymodelstate_free(hosek[i]);
    }
  }
}

// [[Rcpp::export]]
Rcpp::NumericMatrix calculate_raw_prague(Rcpp::NumericVector phi,
                           Rcpp::NumericVector theta,
                           Rcpp::NumericVector elevation,
                           Rcpp::NumericVector albedo,
                           Rcpp::NumericVector altitude,
                           Rcpp::NumericVector visibility,
                           double       azimuth_deg       = 90,
                            unsigned int numbercores       = 1,
                            std::string  prg_dataset       = "")     // Prague only (km)
{
  float az_rad = static_cast<float>(azimuth_deg * M_PI / 180.0);

  // Prague model initialisation
  if (prg_dataset.empty()) {
    Rcpp::stop("prg_dataset must be supplied when model=\"prague\"");
  }

  if (prague_loaded_file != prg_dataset) {
    prague_model.initialize(prg_dataset);
    prague_loaded_file = prg_dataset;
  }
  // auto avail = prague_model.getAvailableData();
  // if (albedo    < avail.albedoMin    || albedo    > avail.albedoMax ||
  //     visibility < avail.visibilityMin|| visibility > avail.visibilityMax) {
  //   Rcpp::stop("albedo or visibility outside dataset range");
  // }

  // Spectral sampling and Hosek state
  constexpr int    N_LAMBDA = 9;
  constexpr double lambda_nm[N_LAMBDA] =
    { 630,680,710, 500,530,560, 460,480,490 };

  // Initialize image buffer
  std::vector<float> img(3 * phi.size(), 0.0f);

  // Render loop (parallel over rows)
  RcppThread::parallelFor(
    0u, size_t(phi.size()),
    [&](size_t t) {
      float theta_rad = theta[t] * static_cast<float>(M_PI)/180;
      float phi_rad = phi[t] * static_cast<float>(M_PI)/180;
      float elev_rad = elevation[t] * static_cast<float>(M_PI)/180;

      if (theta_rad > static_cast<float>(M_PI_2)) {
        return;
      }

      Vec3<float> sun_dir(std::cos(az_rad) * std::cos(elev_rad),  // +X (South)
                          std::sin(elev_rad),                     // +Y (Up)
                          std::sin(az_rad) * std::cos(elev_rad)); // +Z (West)

      Vec3<float> v(std::cos(phi_rad) * std::sin(theta_rad),
                    std::cos(theta_rad),
                    std::sin(phi_rad) * std::sin(theta_rad));

      Vec3<float> v_zup(std::cos(phi_rad) * std::sin(theta_rad),
                        std::sin(phi_rad) * std::sin(theta_rad),
                        std::cos(theta_rad));
      auto P = prague_model.computeParameters(
        /*viewPoint*/   { 0.0, 0.0, altitude[t] },
        /*viewDir  */   toPrague(v_zup),
        /*sunElev  */   elev_rad,
        /*sunAzim  */   az_rad,
        /*visibility*/  visibility[t],
        /*albedo   */   albedo[t]);
        for (int c = 0; c < N_LAMBDA; ++c) {
          double L = prague_model.skyRadiance(P, lambda_nm[c]) +
            prague_model.sunRadiance(P, lambda_nm[c]);
          img[3 * t + c / 3] += static_cast<float>(L / 3.0);
        }
    },
    numbercores);

  Rcpp::NumericMatrix rgb(phi.size(),3);
  for (int p = 0; p < phi.size(); p++) {
    rgb(p,0) = img[3*p];
    rgb(p,1) = img[3*p + 1];
    rgb(p,2) = img[3*p + 2];
  }
  return(rgb);
}
