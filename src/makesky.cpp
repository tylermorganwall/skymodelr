#include <Rcpp.h>
#include <RcppThread.h>

extern "C" {
#include "ArHosekSkyModel.h"
}

#include "PragueSkyModel.h"

#include <Imath/ImathVec.h>
#include <Imath/ImathFun.h>
#include <Imath/ImathBox.h>
#include "CIE1931Data.h"

#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <iterator>

using namespace Imath;

static inline float clamp_float(float x, float a, float b)
{ return (x < a) ? a : (x > b) ? b : x; }

// Forward conversion to Prague’s tiny Vector3 wrapper
static inline PragueSkyModel::Vector3 toPrague(const Vec3<float>& v)
  { return { double(v.x), double(v.y), double(v.z) }; }

// Static instance so the .dat stays in RAM across calls
static PragueSkyModel  prague_model;
static std::string     prague_loaded_file;

static inline void xyz_to_srgb(double X, double Y, double Z,
                               double& R, double& G, double& B)
{
  // IEC 61966-2-1 sRGB (D65), 1931 2°
  const double M[9] = {
     3.2404542, -1.5371385, -0.4985314,
    -0.9692660,  1.8760108,  0.0415560,
     0.0556434, -0.2040259,  1.0572252
  };
  R = M[0]*X + M[1]*Y + M[2]*Z;
  G = M[3]*X + M[4]*Y + M[5]*Z;
  B = M[6]*X + M[7]*Y + M[8]*Z;
}

static inline double clamp_floor(double v, double floor_v = -1e-10) {
  return (v < floor_v) ? floor_v : v;
}

struct CIETristimulus {
  double x;
  double y;
  double z;
};



static const std::vector<double>& default_lambda_low()
{
  static const double lambda_low_arr[] = {
    380, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700
  };
  static const std::vector<double> lambda_low(
    std::begin(lambda_low_arr), std::end(lambda_low_arr));
  return lambda_low;
}

static CIETristimulus interpolate_cie_xyz(double lambda_nm_value)
{
  constexpr int n = sizeof(cie_lambda_nm) / sizeof(double);
  if (lambda_nm_value < cie_lambda_nm[0] ||
      lambda_nm_value > cie_lambda_nm[n - 1]) {
    Rcpp::stop("lambda_nm entries must be within [%f, %f]",
               cie_lambda_nm[0], cie_lambda_nm[n - 1]);
  }
  auto upper = std::lower_bound(std::begin(cie_lambda_nm),
                                std::end(cie_lambda_nm),
                                lambda_nm_value);
  if (upper != std::end(cie_lambda_nm) && *upper == lambda_nm_value) {
    const int idx = int(upper - std::begin(cie_lambda_nm));
    return { cie_xbar[idx], cie_ybar[idx], cie_zbar[idx] };
  }

  const int idx_hi = int(upper - std::begin(cie_lambda_nm));
  const int idx_lo = idx_hi - 1;
  const double lam_lo = cie_lambda_nm[idx_lo];
  const double lam_hi = cie_lambda_nm[idx_hi];
  const double t = (lambda_nm_value - lam_lo) / (lam_hi - lam_lo);
  const double x = cie_xbar[idx_lo] + t * (cie_xbar[idx_hi] - cie_xbar[idx_lo]);
  const double y = cie_ybar[idx_lo] + t * (cie_ybar[idx_hi] - cie_ybar[idx_lo]);
  const double z = cie_zbar[idx_lo] + t * (cie_zbar[idx_hi] - cie_zbar[idx_lo]);
  return { x, y, z };
}

static std::vector<double> trapezoid_weights(const std::vector<double>& lambda_values)
{
  const size_t n = lambda_values.size();
  std::vector<double> weights(n, 0.0);
  if (n == 0) return weights;
  if (n == 1) {
    weights[0] = 1.0;
    return weights;
  }
  for (size_t i = 0; i < n; ++i) {
    if (i == 0) {
      weights[i] = (lambda_values[1] - lambda_values[0]) / 2.0;
    } else if (i == n - 1) {
      weights[i] = (lambda_values[i] - lambda_values[i - 1]) / 2.0;
    } else {
      weights[i] = (lambda_values[i + 1] - lambda_values[i - 1]) / 2.0;
    }
  }
  return weights;
}

// [[Rcpp::export]]
Rcpp::NumericVector makesky_rcpp(
	double       albedo            = 0.5,
    double       turbidity         = 3.0,
    double       elevation         = 10.0,
	double       azimuth_deg       = 90,
    unsigned int resolution        = 2048,
    unsigned int numbercores       = 1,
    std::string  model             = "hosek",  // hosek | prague
    std::string  prg_dataset       = "",
	double       altitude          = 0.0,
    double       visibility        = 50.0, // Prague only (km)
    bool         render_solar_disk = true,
    Rcpp::Nullable<Rcpp::NumericVector> lambda_nm = R_NilValue
) {
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

  // Spectral sampling setup
  std::vector<double> lambda_values;
  const auto& default_values = default_lambda_low();
  if (lambda_nm.isNull()) {
    lambda_values.assign(default_values.begin(), default_values.end());
  } else {
    Rcpp::NumericVector lambda_vec(lambda_nm.get());
    if (lambda_vec.size() == 0) {
      lambda_values.assign(default_values.begin(), default_values.end());
    } else {
      lambda_values.assign(lambda_vec.begin(), lambda_vec.end());
    }
  }
  if (lambda_values.size() < 2) {
    Rcpp::stop("lambda_nm must contain at least two wavelengths");
  }
  for (size_t i = 1; i < lambda_values.size(); ++i) {
    if (!(lambda_values[i] > lambda_values[i - 1])) {
      Rcpp::stop("lambda_nm must be strictly increasing");
    }
  }
  const auto lambda_weights = trapezoid_weights(lambda_values);
  std::vector<CIETristimulus> cie_samples;
  cie_samples.reserve(lambda_values.size());
  for (double lam : lambda_values) {
    cie_samples.push_back(interpolate_cie_xyz(lam));
  }
  const size_t N_LAMBDA = lambda_values.size();

	ArHosekSkyModelState* hosek = nullptr;
	if (model == "hosek") {
		if (elevation < 0.0 || elevation > 90.0) {
			Rcpp::stop("elevation must be in [0,90]");
		}
		if (turbidity < 1.7 || turbidity > 10.0) {
			Rcpp::stop("turbidity must be in [1.7,10]");
		}
		hosek = arhosekskymodelstate_alloc_init(
			elevation * M_PI / 180.0, turbidity, albedo
		);
		if (!hosek) Rcpp::stop("Hosek–Wilkie initialisation failed");

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

		double X = 0.0, Y = 0.0, Z = 0.0;

        if (model == "hosek") {
          float gamma = std::acos(
            clamp_float(v.dot(sun_dir), -1.0f, 1.0f));
          for (size_t c = 0; c < N_LAMBDA; ++c) {
			const double lam = lambda_values[c];
			float L = render_solar_disk ? 
			  arhosekskymodel_solar_radiance(hosek, theta, gamma, lam) : 
			  arhosekskymodel_radiance      (hosek, theta, gamma, lam);
			// Accumulate XYZ
            double Li = static_cast<double>(L);
            X += Li * cie_samples[c].x * lambda_weights[c];
            Y += Li * cie_samples[c].y * lambda_weights[c];
            Z += Li * cie_samples[c].z * lambda_weights[c];
          }
        } else {
		// Prague uses Z-up vector
          Vec3<float> v_zup(std::cos(phi) * std::sin(theta),
                            std::sin(phi) * std::sin(theta),
                            std::cos(theta));
          auto P = prague_model.computeParameters(
            { 0.0, 0.0, altitude }, 
			toPrague(v_zup),
            elev_rad, 
			az_rad, 
			visibility, 
			albedo
		);

          for (size_t i = 0; i < N_LAMBDA; ++i) {
            double Li = prague_model.skyRadiance(P, lambda_values[i]);
            if (render_solar_disk) Li += prague_model.sunRadiance(P, lambda_values[i]);
            X += Li * cie_samples[i].x * lambda_weights[i];
            Y += Li * cie_samples[i].y * lambda_weights[i];
            Z += Li * cie_samples[i].z * lambda_weights[i];
          }
        }
		double R, G, B;
  		xyz_to_srgb(X, Y, Z, R, G, B);
		// Avoid tiny negative noise; don't hard clip highlights here
        R = std::max(0.0, clamp_floor(R));
        G = std::max(0.0, clamp_floor(G));
        B = std::max(0.0, clamp_floor(B));

        const int idx = 3 * (int(t) * nPhi + p);
        img[idx + 0] = static_cast<float>(R);
        img[idx + 1] = static_cast<float>(G);
        img[idx + 2] = static_cast<float>(B);
      }
    },
    numbercores);

  const int width  = nPhi;
  const int height = nTheta;
  const int nPix   = width * height;

  if (model == "hosek") {
	arhosekskymodelstate_free(hosek);
  }

  Rcpp::NumericVector result(nPix * 3);
  result.attr("dim") = Rcpp::IntegerVector::create(height, width, 3);

  for (int t = 0; t < height; ++t) {
    for (int p = 0; p < width; ++p) {
      const int pix_index = t * width + p;
      for (int c = 0; c < 3; ++c) {
        const int arr_index = t + height * (p + width * c);
        result[arr_index] = static_cast<double>(img[3 * pix_index + c]);
      }
    }
  }

  return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix calculate_raw_prague(Rcpp::NumericVector phi,
                           Rcpp::NumericVector theta,
                           Rcpp::NumericVector elevation,
                           Rcpp::NumericVector albedo,
                           Rcpp::NumericVector altitude,
                           Rcpp::NumericVector visibility,
						   Rcpp::NumericVector azimuth,
                           unsigned int numbercores       = 1,
                            std::string  prg_dataset       = "",
                            bool render_solar_disk = true)
{

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
  	  float az_rad = azimuth[t] * static_cast<float>(M_PI) / 180;

      if (theta_rad > static_cast<float>(M_PI_2)) {
        return;
      }

      Vec3<float> v_zup(std::cos(phi_rad) * std::cos(theta_rad),
                        std::sin(phi_rad) * std::cos(theta_rad),
                        std::sin(theta_rad));
      auto P = prague_model.computeParameters(
        /*viewPoint*/   { 0.0, 0.0, altitude[t] },
        /*viewDir  */   toPrague(v_zup),
        /*sunElev  */   elev_rad,
        /*sunAzim  */   az_rad,
        /*visibility*/  visibility[t],
        /*albedo   */   albedo[t]);
        for (int c = 0; c < N_LAMBDA; ++c) {
          if(render_solar_disk) {
            double L = prague_model.skyRadiance(P, lambda_nm[c]) +
              prague_model.sunRadiance(P, lambda_nm[c]);
            img[3 * t + c / 3] += static_cast<float>(L / 3.0);
          } else {
            double L = prague_model.skyRadiance(P, lambda_nm[c]);
            img[3 * t + c / 3] += static_cast<float>(L / 3.0);
          }
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
