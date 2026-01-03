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

static inline double clamp_double(double x, double a, double b)
{ return (x < a) ? a : (x > b) ? b : x; }

// Forward conversion to Prague’s tiny Vector3 wrapper
static inline PragueSkyModel::Vector3 toPrague(const Vec3<double>& v)
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

struct CIETristimulus {
  double x;
  double y;
  double z;
};



static const std::vector<double>& default_lambda_low()
{
  static const double lambda_low_arr[] = {
    360, 380, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700
  };
  static const std::vector<double> lambda_low(
    std::begin(lambda_low_arr), std::end(lambda_low_arr));
  return lambda_low;
}

static CIETristimulus interpolate_cie_xyz(double lambda_nm_value)
{
  constexpr int n = CIE_N;
  if (lambda_nm_value < cie_lambda_nm[0] ||
      lambda_nm_value > cie_lambda_nm[n - 1]) {
    Rcpp::stop("lambda_nm entries must be within [%f, %f]",
               cie_lambda_nm[0], cie_lambda_nm[n - 1]);
  }
  const auto begin = std::begin(cie_lambda_nm);
  const auto end = std::end(cie_lambda_nm);

  auto upper = std::lower_bound(begin, end, lambda_nm_value);

  if (upper == begin) {
    return { cie_xbar[0], cie_ybar[0], cie_zbar[0] };
  }
  if (upper == end) {
    const int last = n - 1;
    return { cie_xbar[last], cie_ybar[last], cie_zbar[last] };
  }

  const int idx_hi = int(upper - begin);
  const int idx_lo = idx_hi - 1;

  if (*upper == lambda_nm_value) {
    return { cie_xbar[idx_hi], cie_ybar[idx_hi], cie_zbar[idx_hi] };
  }

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
    Rcpp::Nullable<Rcpp::NumericVector> lambda_nm = R_NilValue,
	bool         below_horizon     = true
) {
  if (albedo < 0.0 || albedo > 1.0) {
    Rcpp::stop("albedo must be in [0,1]");
  }
  if (resolution == 0) {
    Rcpp::stop("resolution must be ≥ 1");
  }

  double elev_rad = elevation * M_PI / 180.0;
  double az_rad = azimuth_deg * M_PI / 180.0;

  // Prague model initialisation
  PragueSkyModel::AvailableData avail{};
  if (model == "prague") {
    if (prg_dataset.empty()) {
      Rcpp::stop("prg_dataset must be supplied when model=\"prague\"");
    }

    if (prague_loaded_file != prg_dataset) {
      prague_model.initialize(prg_dataset);
      prague_loaded_file = prg_dataset;
    }
    avail = prague_model.getAvailableData();
    if (albedo    < avail.albedoMin    || albedo    > avail.albedoMax ||
        visibility < avail.visibilityMin|| visibility > avail.visibilityMax) {
      Rcpp::stop("albedo or visibility outside dataset range");
    }
  } else if (model != "hosek") {
    Rcpp::stop("unknown model \"%s\"", model.c_str());
  }

  // Spectral sampling setup
  std::vector<double> lambda_values;
  std::vector<double> lambda_weights;
  std::vector<CIETristimulus> cie_samples;

  if (model == "prague") {
    // Prague sampling is locked to dataset channel centers within CIE bounds (380–740 nm).
    lambda_values.reserve(avail.channels);
    for (int i = 0; i < avail.channels; ++i) {
      const double lam = avail.channelStart + (i + 0.5) * avail.channelWidth;
      if (lam >= 380.0 && lam <= 740.0) {
        lambda_values.push_back(lam);
      }
    }
    if (lambda_values.empty()) {
      Rcpp::stop("No Prague spectral channels fall within [380, 740] nm.");
    }
    // Prague channels are evenly spaced bands.
    lambda_weights.assign(lambda_values.size(), avail.channelWidth);
    cie_samples.reserve(lambda_values.size());
    for (double lam : lambda_values) {
      cie_samples.push_back(interpolate_cie_xyz(lam));
    }

  } else {
    // Hosek path (unchanged behavior)
    const auto &default_values = default_lambda_low();

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

    lambda_weights = trapezoid_weights(lambda_values);

    cie_samples.reserve(lambda_values.size());
    for (double lam : lambda_values) {
      cie_samples.push_back(interpolate_cie_xyz(lam));
    }
  }

  const size_t N_LAMBDA = lambda_values.size();

  ArHosekSkyModelState *hosek = nullptr;
  if (model == "hosek") {
    if (elevation < 0.0 || elevation > 90.0) {
      Rcpp::stop("elevation must be in [0,90]");
    }
    if (turbidity < 1.7 || turbidity > 10.0) {
      Rcpp::stop("turbidity must be in [1.7,10]");
    }
    hosek = arhosekskymodelstate_alloc_init(elevation * M_PI / 180.0, turbidity,
                                            albedo);
    if (!hosek)
      Rcpp::stop("Hosek–Wilkie initialisation failed");

  } else {
    if (elevation < -4.2 || elevation > 90.0) {
      Rcpp::stop("elevation must be in [-4.2, 90]");
    }
	}

  // Initialize image buffer
  const int nTheta = resolution;          // rows  (latitude)
  const int nPhi   = 2 * resolution;      // cols  (longitude)
  std::vector<double> img(3 * nTheta * nPhi, 0.0);
  std::vector<double> img_band(nTheta * nPhi, 0.0);

  //   Vec3<double> sun_dir(0.0, std::sin(elev_rad), std::cos(elev_rad));
  Vec3<double> sun_dir(std::cos(az_rad) * std::cos(elev_rad),  // +X (South)
                       std::sin(elev_rad),                     // +Y (Up)
                       std::sin(az_rad) * std::cos(elev_rad)); // +Z (West)
  // Render loop (parallel over rows)
  RcppThread::parallelFor(
    0u, size_t(nTheta),
    [&](size_t t) {
      double theta = (static_cast<double>(t) + 0.5) /
        static_cast<double>(nTheta) * M_PI;

      if (!below_horizon && theta > M_PI_2) {
        return;
      }

      for (int p = 0; p < nPhi; ++p) {
        double phi = (static_cast<double>(p) + 0.5) / nPhi *
          2.0 * M_PI;

        Vec3<double> v(std::cos(phi) * std::sin(theta),
                       std::cos(theta),
                       std::sin(phi) * std::sin(theta));

		double X = 0.0, Y = 0.0, Z = 0.0;
        double L_band = 0.0;
        const bool view_above_horizon = theta <= M_PI_2;

        if (model == "hosek") {
          double gamma = std::acos(
            clamp_double(v.dot(sun_dir), -1.0, 1.0));
          for (size_t c = 0; c < N_LAMBDA; ++c) {
				const double lam = lambda_values[c];
				double L = (render_solar_disk && view_above_horizon) ?
				  arhosekskymodel_solar_radiance(hosek, theta, gamma, lam) :
				  arhosekskymodel_radiance      (hosek, theta, gamma, lam);
			// Accumulate XYZ
            double Li = static_cast<double>(L);
            const double d_lambda = lambda_weights[c];
            X += Li * cie_samples[c].x * d_lambda;
            Y += Li * cie_samples[c].y * d_lambda;
            Z += Li * cie_samples[c].z * d_lambda;
            L_band += Li * d_lambda;
          }
        } else {
        // Prague uses Z-up vector
          Vec3<double> v_zup(std::cos(phi) * std::sin(theta),
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
            if (render_solar_disk && view_above_horizon) Li += prague_model.sunRadiance(P, lambda_values[i]);
            const double d_lambda = lambda_weights[i];
            X += Li * cie_samples[i].x * d_lambda;
            Y += Li * cie_samples[i].y * d_lambda;
            Z += Li * cie_samples[i].z * d_lambda;
            L_band += Li * d_lambda;
          }
        }
		if (!std::isfinite(X) || !std::isfinite(Y) || !std::isfinite(Z)) {
			X = 0.0;
			Y = 0.0;
			Z = 0.0;
        }
		double R, G, B;
  		xyz_to_srgb(X, Y, Z, R, G, B);

		if (!std::isfinite(R) || !std::isfinite(G) ||
			!std::isfinite(B)) {
			R = 0.0;
			G = 0.0;
			B = 0.0;
		}

        const int idx = 3 * (int(t) * nPhi + p);
        img[idx + 0] = R;
        img[idx + 1] = G;
        img[idx + 2] = B;
        img_band[int(t) * nPhi + p] = L_band;
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
          result[arr_index] = img[3 * pix_index + c];
        }
      }
    }

  Rcpp::NumericVector band(nPix);
  band.attr("dim") = Rcpp::IntegerVector::create(height, width);
  for (int t = 0; t < height; ++t) {
    for (int p = 0; p < width; ++p) {
      const int pix_index = t * width + p;
      band[pix_index] = img_band[pix_index];
    }
  }
  result.attr("L_band") = band;

  return result;
}

// [[Rcpp::export]]
double calculate_sun_brightness_rcpp(
    double      albedo       = 0.5,
    double      turbidity    = 3.0,
    double      elevation    = 10.0,
    double      azimuth_deg  = 90.0,
    std::string model        = "hosek",
    std::string prg_dataset  = "",
    double      altitude     = 0.0,
    double      visibility   = 50.0,
    Rcpp::Nullable<Rcpp::NumericVector> lambda_nm = R_NilValue
) {
  if (albedo < 0.0 || albedo > 1.0) {
    Rcpp::stop("albedo must be in [0,1]");
  }

  const double elev_rad = elevation * M_PI / 180.0;
  const double az_rad = azimuth_deg * M_PI / 180.0;

  // Prague model initialisation
  PragueSkyModel::AvailableData avail{};
  if (model == "prague") {
    if (prg_dataset.empty()) {
      Rcpp::stop("prg_dataset must be supplied when model=\"prague\"");
    }

    if (prague_loaded_file != prg_dataset) {
      prague_model.initialize(prg_dataset);
      prague_loaded_file = prg_dataset;
    }
    avail = prague_model.getAvailableData();
    if (albedo < avail.albedoMin || albedo > avail.albedoMax ||
        visibility < avail.visibilityMin || visibility > avail.visibilityMax) {
      Rcpp::stop("albedo or visibility outside dataset range");
    }
  } else if (model != "hosek") {
    Rcpp::stop("unknown model \"%s\"", model.c_str());
  }

  // Spectral sampling setup
  std::vector<double> lambda_values;
  std::vector<double> lambda_weights;
  std::vector<double> cie_y_samples;

  if (model == "prague") {
    // Prague sampling is locked to dataset channel centers within CIE bounds (380–740 nm).
    lambda_values.reserve(avail.channels);
    for (int i = 0; i < avail.channels; ++i) {
      const double lam = avail.channelStart + (i + 0.5) * avail.channelWidth;
      if (lam >= 380.0 && lam <= 740.0) {
        lambda_values.push_back(lam);
      }
    }
    if (lambda_values.empty()) {
      Rcpp::stop("No Prague spectral channels fall within [380, 740] nm.");
    }
    // Prague channels are evenly spaced bands.
    lambda_weights.assign(lambda_values.size(), avail.channelWidth);
    cie_y_samples.reserve(lambda_values.size());
    for (double lam : lambda_values) {
      cie_y_samples.push_back(interpolate_cie_xyz(lam).y);
    }
  } else {
    // Hosek path (unchanged behavior)
    const auto &default_values = default_lambda_low();

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

    lambda_weights = trapezoid_weights(lambda_values);

    cie_y_samples.reserve(lambda_values.size());
    for (double lam : lambda_values) {
      cie_y_samples.push_back(interpolate_cie_xyz(lam).y);
    }
  }

  const size_t N_LAMBDA = lambda_values.size();
  double Y = 0.0;

  if (model == "hosek") {
    if (elevation < 0.0 || elevation > 90.0) {
      Rcpp::stop("elevation must be in [0,90]");
    }
    if (turbidity < 1.7 || turbidity > 10.0) {
      Rcpp::stop("turbidity must be in [1.7,10]");
    }

    ArHosekSkyModelState *hosek = arhosekskymodelstate_alloc_init(
      elev_rad, turbidity, albedo);
    if (!hosek) {
      Rcpp::stop("Hosek–Wilkie initialisation failed");
    }

    const double theta = M_PI_2 - elev_rad;
    const double gamma = 0.0;
    for (size_t c = 0; c < N_LAMBDA; ++c) {
      const double lam = lambda_values[c];
      const double L = arhosekskymodel_solar_radiance(hosek, theta, gamma, lam);
      const double d_lambda = lambda_weights[c];
      Y += L * cie_y_samples[c] * d_lambda;
    }
    arhosekskymodelstate_free(hosek);
  } else {
    if (elevation < -4.2 || elevation > 90.0) {
      Rcpp::stop("elevation must be in [-4.2, 90]");
    }

    // Sample sun radiance along the sun direction.
    Vec3<double> sun_dir(std::cos(az_rad) * std::cos(elev_rad),
                         std::sin(az_rad) * std::cos(elev_rad),
                         std::sin(elev_rad));
    auto P = prague_model.computeParameters(
      { 0.0, 0.0, altitude },
      toPrague(sun_dir),
      elev_rad,
      az_rad,
      visibility,
      albedo
    );

    for (size_t i = 0; i < N_LAMBDA; ++i) {
      const double Li = prague_model.sunRadiance(P, lambda_values[i]);
      const double d_lambda = lambda_weights[i];
      Y += Li * cie_y_samples[i] * d_lambda;
    }
  }

  if (!std::isfinite(Y)) {
    Y = 0.0;
  }

  return Y;
}

// [[Rcpp::export]]
double calculate_sun_radiance_band_rcpp(
    double      albedo       = 0.5,
    double      turbidity    = 3.0,
    double      elevation    = 10.0,
    double      azimuth_deg  = 90.0,
    std::string model        = "hosek",
    std::string prg_dataset  = "",
    double      altitude     = 0.0,
    double      visibility   = 50.0,
    Rcpp::Nullable<Rcpp::NumericVector> lambda_nm = R_NilValue
) {
  if (albedo < 0.0 || albedo > 1.0) {
    Rcpp::stop("albedo must be in [0,1]");
  }

  const double elev_rad = elevation * M_PI / 180.0;
  const double az_rad = azimuth_deg * M_PI / 180.0;

  // Prague model initialisation
  PragueSkyModel::AvailableData avail{};
  if (model == "prague") {
    if (prg_dataset.empty()) {
      Rcpp::stop("prg_dataset must be supplied when model=\"prague\"");
    }

    if (prague_loaded_file != prg_dataset) {
      prague_model.initialize(prg_dataset);
      prague_loaded_file = prg_dataset;
    }
    avail = prague_model.getAvailableData();
    if (albedo < avail.albedoMin || albedo > avail.albedoMax ||
        visibility < avail.visibilityMin || visibility > avail.visibilityMax) {
      Rcpp::stop("albedo or visibility outside dataset range");
    }
  } else if (model != "hosek") {
    Rcpp::stop("unknown model \"%s\"", model.c_str());
  }

  // Spectral sampling setup
  std::vector<double> lambda_values;
  std::vector<double> lambda_weights;

  if (model == "prague") {
    lambda_values.reserve(avail.channels);
    for (int i = 0; i < avail.channels; ++i) {
      const double lam = avail.channelStart + (i + 0.5) * avail.channelWidth;
      if (lam >= 380.0 && lam <= 740.0) {
        lambda_values.push_back(lam);
      }
    }
    if (lambda_values.empty()) {
      Rcpp::stop("No Prague spectral channels fall within [380, 740] nm.");
    }
    lambda_weights.assign(lambda_values.size(), avail.channelWidth);
  } else {
    const auto &default_values = default_lambda_low();

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

    lambda_weights = trapezoid_weights(lambda_values);
  }

  const size_t N_LAMBDA = lambda_values.size();
  double L_band = 0.0;

  if (model == "hosek") {
    if (elevation < 0.0 || elevation > 90.0) {
      Rcpp::stop("elevation must be in [0,90]");
    }
    if (turbidity < 1.7 || turbidity > 10.0) {
      Rcpp::stop("turbidity must be in [1.7,10]");
    }

    ArHosekSkyModelState *hosek = arhosekskymodelstate_alloc_init(
      elev_rad, turbidity, albedo);
    if (!hosek) {
      Rcpp::stop("Hosek–Wilkie initialisation failed");
    }

    const double theta = M_PI_2 - elev_rad;
    const double gamma = 0.0;
    for (size_t c = 0; c < N_LAMBDA; ++c) {
      const double lam = lambda_values[c];
      const double L = arhosekskymodel_solar_radiance(hosek, theta, gamma, lam);
      const double d_lambda = lambda_weights[c];
      L_band += L * d_lambda;
    }
    arhosekskymodelstate_free(hosek);
  } else {
    if (elevation < -4.2 || elevation > 90.0) {
      Rcpp::stop("elevation must be in [-4.2, 90]");
    }

    Vec3<double> sun_dir(std::cos(az_rad) * std::cos(elev_rad),
                         std::sin(az_rad) * std::cos(elev_rad),
                         std::sin(elev_rad));
    auto P = prague_model.computeParameters(
      { 0.0, 0.0, altitude },
      toPrague(sun_dir),
      elev_rad,
      az_rad,
      visibility,
      albedo
    );

    for (size_t i = 0; i < N_LAMBDA; ++i) {
      const double Li = prague_model.sunRadiance(P, lambda_values[i]);
      const double d_lambda = lambda_weights[i];
      L_band += Li * d_lambda;
    }
  }

  if (!std::isfinite(L_band)) {
    L_band = 0.0;
  }

  return L_band;
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
  auto avail = prague_model.getAvailableData();

  // Spectral sampling: lock to Prague dataset channel centers within CIE range (380–740 nm).
  std::vector<double> lambda_values;
  lambda_values.reserve(avail.channels);
  for (int i = 0; i < avail.channels; ++i) {
    const double lam = avail.channelStart + (i + 0.5) * avail.channelWidth;
    if (lam >= 380.0 && lam <= 740.0) {
      lambda_values.push_back(lam);
    }
  }
  if (lambda_values.empty()) {
    Rcpp::stop("No Prague spectral channels fall within [380, 740] nm.");
  }
  const size_t N_LAMBDA = lambda_values.size();
  // Constant band weights equal to channel width for rectangular integration.
  std::vector<double> lambda_weights(N_LAMBDA, avail.channelWidth);
  std::vector<CIETristimulus> cie_samples;
  cie_samples.reserve(N_LAMBDA);
  for (double lam : lambda_values) {
    cie_samples.push_back(interpolate_cie_xyz(lam));
  }

  // Initialize image buffer
  std::vector<double> img(3 * phi.size(), 0.0);
  std::vector<double> band(phi.size(), 0.0);

  // Render loop (parallel over rows)
  RcppThread::parallelFor(
    0u, size_t(phi.size()),
    [&](size_t t) {
      double theta_rad = theta[t] * M_PI/180;
      double phi_rad = phi[t] * M_PI/180;
      double elev_rad = elevation[t] * M_PI/180;
  	  double az_rad = azimuth[t] * M_PI / 180;

      if (theta_rad > M_PI_2) {
        return;
      }

      Vec3<double> v_zup(std::cos(phi_rad) * std::cos(theta_rad),
                         std::sin(phi_rad) * std::cos(theta_rad),
                         std::sin(theta_rad));
      auto P = prague_model.computeParameters(
        /*viewPoint*/   { 0.0, 0.0, altitude[t] },
        /*viewDir  */   toPrague(v_zup),
        /*sunElev  */   elev_rad,
        /*sunAzim  */   az_rad,
        /*visibility*/  visibility[t],
        /*albedo   */   albedo[t]);
        double X = 0.0, Y = 0.0, Z = 0.0;
        double L_band = 0.0;
        const bool view_above_horizon = theta_rad <= M_PI_2;
        for (size_t c = 0; c < N_LAMBDA; ++c) {
          const double lam = lambda_values[c];
          double Li = prague_model.skyRadiance(P, lam);
          if (render_solar_disk && view_above_horizon) Li += prague_model.sunRadiance(P, lam);
          const double d_lambda = lambda_weights[c];
          X += Li * cie_samples[c].x * d_lambda;
          Y += Li * cie_samples[c].y * d_lambda;
          Z += Li * cie_samples[c].z * d_lambda;
          L_band += Li * d_lambda;
        }
        if (!std::isfinite(X) || !std::isfinite(Y) || !std::isfinite(Z)) {
          X = 0.0;
          Y = 0.0;
          Z = 0.0;
        }
        double R, G, B;
        xyz_to_srgb(X, Y, Z, R, G, B);
        if (!std::isfinite(R) || !std::isfinite(G) || !std::isfinite(B)) {
          R = 0.0;
          G = 0.0;
          B = 0.0;
        }
		img[3 * t + 0] = R;
		img[3 * t + 1] = G;
		img[3 * t + 2] = B;
        band[t] = L_band;
    },
    numbercores);

  Rcpp::NumericMatrix rgb(phi.size(),3);
  for (int p = 0; p < phi.size(); p++) {
    rgb(p,0) = img[3*p];
    rgb(p,1) = img[3*p + 1];
    rgb(p,2) = img[3*p + 2];
  }
  Rcpp::NumericVector band_vals(band.begin(), band.end());
  rgb.attr("L_band") = band_vals;
  return(rgb);
}
