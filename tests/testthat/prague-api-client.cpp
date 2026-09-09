// [[Rcpp::depends(skymodelr, RcppThread)]]
#include <Rcpp.h>
#include <skymodelr/prague.hpp>
#include <future>
#include <array>
#include <fstream>
using Model = skymodelr::PragueSkyModel;

// [[Rcpp::export]]
Rcpp::List prague_api_errors(std::string missing, std::string malformed) {
  Model model;
  auto fails = [](auto f) { try { f(); return false; } catch (const std::exception &) { return true; } };
  auto get = reinterpret_cast<skymodelr_prague_get_api_v1_fn>(R_GetCCallable("skymodelr", "prague_get_api_v1"));
  const auto *api = get();
  char message[4] = {'x','x','x','x'};
  bool missing_file = !api->create(missing.c_str(), 50, message, sizeof(message));
  return Rcpp::List::create(
    Rcpp::Named("version") = api->abi_version,
    Rcpp::Named("null_query") = fails([&] { model.getAvailableData(); }),
    Rcpp::Named("missing_file") = missing_file,
    Rcpp::Named("terminated_error") = message[3] == '\0',
    Rcpp::Named("malformed_file") = fails([&] { model.initialize(malformed, 50); }));
}

// [[Rcpp::export]]
Rcpp::List prague_api_probe(std::string filename, std::string missing) {
  Model model(filename, 50);
  auto metadata = model.getAvailableData();
  std::array<double, 9> wavelengths{100, 380, 420, 500, 560, 620, 700, 740, 3000};
  std::array<double, 9> radiance, transmission;
  auto params = model.computeParameters({0, 0, 0}, {0, 0, 1}, .3, .7, 50, .5);
  model.skyRadianceSpectrum(params, wavelengths.data(), wavelengths.size(), radiance.data());
  model.transmittanceSpectrum(params, wavelengths.data(), wavelengths.size(), 10000, transmission.data());
  double batch_error = 0;
  for (size_t i = 0; i < wavelengths.size(); ++i) {
    batch_error = std::max(batch_error, std::abs(radiance[i] - model.skyRadiance(params, wavelengths[i])));
    batch_error = std::max(batch_error, std::abs(transmission[i] - model.transmittance(params, wavelengths[i], 10000)));
  }
  // Resolve symbols and initialize coefficients before entering native workers.
  std::array<std::future<double>, 4> futures;
  for (auto &f : futures) f = std::async(std::launch::async, [&] {
    double error = 0;
    std::array<double, 9> r, t;
    for (int i = 0; i < 100; ++i) {
      model.skyRadianceSpectrum(params, wavelengths.data(), wavelengths.size(), r.data());
      model.transmittanceSpectrum(params, wavelengths.data(), wavelengths.size(), 10000, t.data());
      for (size_t j = 0; j < wavelengths.size(); ++j)
        error = std::max({error, std::abs(r[j] - radiance[j]), std::abs(t[j] - transmission[j])});
    }
    return error;
  });
  double thread_error = 0;
  for (auto &f : futures) thread_error = std::max(thread_error, f.get());
  bool failed_initialize = false;
  try { model.initialize(missing, 50); } catch (const std::exception &) { failed_initialize = true; }
  bool preserved = model.getAvailableData().channels == metadata.channels;
  params.gamma = 0;
  double intrinsic = model.sunRadiance(params, 560, false), attenuated = model.sunRadiance(params, 560);
  auto get = reinterpret_cast<skymodelr_prague_get_api_v1_fn>(R_GetCCallable("skymodelr", "prague_get_api_v1"));
  const auto *api = get(); char error[512];
  auto *handle = api->create(filename.c_str(), 50, error, sizeof(error));
  if (!handle) Rcpp::stop(error);
  // Provider contains errors: malformed public C inputs do not escape as C++ exceptions.
  double p[3] = {0, 0, 0}, w[3] = {0, 0, 0}, result;
  bool invalid_direction = !api->parameters(handle, p, w, .3, .7, 50, .5, &params, error, sizeof(error));
  bool invalid_kind = !api->spectrum(handle, 99, &params, wavelengths.data(), 1, 1, 1, &result, error, sizeof(error));
  auto nan = params; nan.gamma = std::numeric_limits<double>::quiet_NaN();
  bool invalid_parameters = !api->spectrum(handle, SKYMODELR_PRAGUE_SKY, &nan, wavelengths.data(), 1, 1, 1, &result, error, sizeof(error));
  api->destroy(handle);
  return Rcpp::List::create(
    Rcpp::Named("batch_error") = batch_error, Rcpp::Named("thread_error") = thread_error,
    Rcpp::Named("radiance") = Rcpp::NumericVector(radiance.begin(), radiance.end()),
    Rcpp::Named("transmission") = Rcpp::NumericVector(transmission.begin(), transmission.end()),
    Rcpp::Named("intrinsic_sun") = intrinsic, Rcpp::Named("attenuated_sun") = attenuated,
    Rcpp::Named("elevation_max") = metadata.elevationMax,
    Rcpp::Named("channels") = metadata.channels, Rcpp::Named("memory") = double(model.memoryUsage()),
    Rcpp::Named("failed_initialize") = failed_initialize, Rcpp::Named("preserved") = preserved,
    Rcpp::Named("invalid_direction") = invalid_direction, Rcpp::Named("invalid_kind") = invalid_kind,
    Rcpp::Named("invalid_parameters") = invalid_parameters);
}
