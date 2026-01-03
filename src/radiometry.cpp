#include <Rcpp.h>

#include "CIE1931Data.h"
#include "hosek/cie.h"
#include "math/units.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <string>

static inline std::string to_lower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return value;
}

static inline double spd_blackbody(double lambda_nm, double temperature_k) {
  const double lambda_m = lambda_nm * 1e-9;
  const double c2 = 1.4387769e-2; // hc/k
  const double expo = c2 / (lambda_m * temperature_k);
  const double denom = std::exp(expo) - 1.0;
  if (!std::isfinite(denom) || denom <= 0.0) {
    return 0.0;
  }
  return 1.0 / (std::pow(lambda_m, 5.0) * denom);
}

static inline double spd_d65(double lambda_nm) {
  constexpr double x_d = 0.31271;
  constexpr double y_d = 0.32902;
  const double denom = 0.0241 + 0.2562 * x_d - 0.7341 * y_d;
  const double m1 = (-1.3515 - 1.7703 * x_d + 5.9114 * y_d) / denom;
  const double m2 = (0.03 - 31.4424 * x_d + 30.0717 * y_d) / denom;
  const auto &s0 = S0_function();
  const auto &s1 = S1_function();
  const auto &s2 = S2_function();
  const auto wl = lambda_nm * nm;
  const double s0_val = s0(wl)();
  const double s1_val = s1(wl)();
  const double s2_val = s2(wl)();
  return s0_val + m1 * s1_val + m2 * s2_val;
}

// [[Rcpp::export]]
Rcpp::List cie_1931_2deg_rcpp() {
  Rcpp::NumericVector lambda_nm(CIE_N);
  Rcpp::NumericVector x_bar(CIE_N);
  Rcpp::NumericVector y_bar(CIE_N);
  Rcpp::NumericVector z_bar(CIE_N);
  for (int i = 0; i < CIE_N; ++i) {
    lambda_nm[i] = cie_lambda_nm[i];
    x_bar[i] = cie_xbar[i];
    y_bar[i] = cie_ybar[i];
    z_bar[i] = cie_zbar[i];
  }
  return Rcpp::List::create(
    Rcpp::Named("lambda_nm") = lambda_nm,
    Rcpp::Named("x_bar") = x_bar,
    Rcpp::Named("y_bar") = y_bar,
    Rcpp::Named("z_bar") = z_bar
  );
}

// [[Rcpp::export]]
Rcpp::NumericVector spd_values_rcpp(Rcpp::NumericVector lambda_nm,
                                    std::string spd_type = "D65") {
  const std::string key = to_lower(spd_type);
  Rcpp::NumericVector spd(lambda_nm.size());
  if (key == "d65") {
    for (R_xlen_t i = 0; i < lambda_nm.size(); ++i) {
      spd[i] = spd_d65(lambda_nm[i]);
    }
  } else if (key == "bb5778") {
    for (R_xlen_t i = 0; i < lambda_nm.size(); ++i) {
      spd[i] = spd_blackbody(lambda_nm[i], 5778.0);
    }
  } else {
    Rcpp::stop("Unknown spd_type \"%s\"", spd_type.c_str());
  }
  return spd;
}
