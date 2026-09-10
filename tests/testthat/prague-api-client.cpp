// [[Rcpp::depends(skymodelr, RcppThread)]]
#include <Rcpp.h>
#include <skymodelr/prague.hpp>
#include <future>
#include <array>
#include <fstream>
#include <memory>
using Model = skymodelr::PragueSkyModel;

// [[Rcpp::export]]
Rcpp::List prague_api_errors(std::string missing, std::string malformed) {
  Model model;
  auto fails = [](auto f) { try { f(); return false; } catch (const std::exception &) { return true; } };
  auto get = reinterpret_cast<skymodelr_prague_get_api_fn>(R_GetCCallable("skymodelr", "prague_get_api"));
  const auto *api = get();
  char message[4] = {'x','x','x','x'};
  bool missing_file = !api->create(missing.c_str(), 50, 0, 0, 512, message, sizeof(message));
  return Rcpp::List::create(
    Rcpp::Named("version") = api->abi_version,
    Rcpp::Named("valid_table_size") = api->struct_size >= sizeof(skymodelr_prague_api),
    Rcpp::Named("null_query") = fails([&] { model.getAvailableData(); }),
    Rcpp::Named("missing_file") = missing_file,
    Rcpp::Named("terminated_error") = message[3] == '\0',
    Rcpp::Named("malformed_file") = fails([&] { model.initialize(malformed, 50); }));
}

// [[Rcpp::export]]
Rcpp::List prague_api_acceleration(std::string filename) {
  auto get = reinterpret_cast<skymodelr_prague_get_api_fn>(
    R_GetCCallable("skymodelr", "prague_get_api"));
  const auto *api = get();
  char error[512];
  using Handle = std::unique_ptr<skymodelr_prague_handle,
                                decltype(api->destroy)>;
  Handle reference(api->create(filename.c_str(), 50, 0, 0, 512, error, sizeof(error)), api->destroy);
  if (!reference) Rcpp::stop(error);
  bool invalid_cache = !api->create(filename.c_str(), 50, 2, 1, 512, error, sizeof(error));
  bool invalid_table = !api->create(filename.c_str(), 50, 1, -1, 512, error, sizeof(error));
  double largest_error = 0;
  size_t differing_values = 0;
  Rcpp::NumericVector memory(4);
  const std::array<double, 17> wavelengths{100, 380, 400, 420, 450, 480, 500, 520,
                                         550, 560, 600, 620, 650, 680, 700, 740, 3000};
  auto compare = [&](const Model &model, int query) {
    double height = (query % 5) * 3000.;
    double direction[3] = {1, .1, (query % 11 - 5) * .01};
    double position[3] = {0, 0, height};
    double visibility = query % 3 == 0 ? 20 : query % 3 == 1 ? 50 : 131.8;
    Model::Parameters p;
    if (!api->parameters(reference.get(), position, direction, -.03, 1.2,
                           visibility, .5, &p, error, sizeof(error))) Rcpp::stop(error);
    // Include scalar, normal visible batches, and batches larger than the cache.
    for (size_t count : {size_t(1), size_t(13), wavelengths.size()}) {
      std::array<double, 17> expected{}, actual{};
      for (double distance : {-1., 0., .01, 30., 100., 10000., 1e30}) {
        int kind = distance < 0 ? SKYMODELR_PRAGUE_SKY : SKYMODELR_PRAGUE_TRANSMITTANCE;
        if (!api->spectrum(reference.get(), kind, &p, wavelengths.data(), count,
                              std::max(0., distance), 1, expected.data(), error, sizeof(error)))
          Rcpp::stop(error);
        // Repeat the sky query to exercise an exact cache hit as well as misses.
        for (int repeat = 0; repeat < 2; ++repeat) {
          if (distance < 0) model.skyRadianceSpectrum(p, wavelengths.data(), count, actual.data());
          else model.transmittanceSpectrum(p, wavelengths.data(), count, distance, actual.data());
          for (size_t i = 0; i < count; ++i) {
            differing_values += actual[i] != expected[i];
            largest_error = std::max(largest_error, std::abs(actual[i] - expected[i]));
          }
        }
      }
    }
  };
  for (int mode = 0; mode < 4; ++mode) {
    Model model(filename, 50, mode & 1, mode & 2);
    memory[mode] = double(model.memoryUsage());
    for (int query = 0; query < 165; ++query) compare(model, query);
  }
  // Reinitialize a populated cache with different coefficients. The same query
  // must now agree with a new uncached handle, including when addresses recur.
  Model reloaded(filename, 50, true, false);
  compare(reloaded, 1);
  reloaded.initialize(filename, 100, true, false);
  reference.reset(api->create(filename.c_str(), 100, 0, 0, 512, error, sizeof(error)));
  if (!reference) Rcpp::stop(error);
  compare(reloaded, 1);
  return Rcpp::List::create(Rcpp::Named("version") = api->abi_version,
    Rcpp::Named("invalid_cache") = invalid_cache, Rcpp::Named("invalid_table") = invalid_table,
    Rcpp::Named("max_error") = largest_error, Rcpp::Named("differing_values") = double(differing_values),
    Rcpp::Named("memory") = memory);
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
  auto get = reinterpret_cast<skymodelr_prague_get_api_fn>(R_GetCCallable("skymodelr", "prague_get_api"));
  const auto *api = get(); char error[512];
  auto *handle = api->create(filename.c_str(), 50, 0, 0, 512, error, sizeof(error));
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

// [[Rcpp::export]]
Rcpp::List prague_api_table_limits(std::string filename) {
  Model reference(filename, 50, false, false);
  const double wavelengths[] = {380, 440, 550, 620, 700, 740};
  auto p = reference.computeParameters({0, 0, 5000}, {1, 0, .01}, .3, .7, 50, .5);
  std::array<double, 6> expected{};
  reference.transmittanceSpectrum(p, wavelengths, 6, 10000, expected.data());
  const double baseline = double(reference.memoryUsage());
  double required;
  {
    Model unlimited(filename, 50, false, true, std::numeric_limits<double>::infinity());
    required = double(unlimited.memoryUsage()) - baseline;
  }
  double required_mib = required / (1024 * 1024);
  Rcpp::NumericVector caps = Rcpp::NumericVector::create(
    0., 1., required_mib - 1. / (1024 * 1024), required_mib,
    512., 1024., std::numeric_limits<double>::infinity());
  Rcpp::NumericVector allocated(caps.size());
  double max_error = 0;
  for (R_xlen_t i = 0; i < caps.size(); ++i) {
    Model model(filename, 50, false, true, caps[i]);
    allocated[i] = double(model.memoryUsage()) - baseline;
    std::array<double, 6> actual{};
    model.transmittanceSpectrum(p, wavelengths, 6, 10000, actual.data());
    for (size_t j = 0; j < actual.size(); ++j)
      max_error = std::max(max_error, std::abs(actual[j] - expected[j]));
  }
  auto get = reinterpret_cast<skymodelr_prague_get_api_fn>(
    R_GetCCallable("skymodelr", "prague_get_api"));
  const auto *api = get();
  char error[512];
  bool invalid = true;
  for (double cap : {-1., -std::numeric_limits<double>::infinity(),
                     std::numeric_limits<double>::quiet_NaN()}) {
    auto *handle = api->create(filename.c_str(), 50, 0, 1, cap, error, sizeof(error));
    invalid &= handle == nullptr;
    if (handle) api->destroy(handle);
  }
  return Rcpp::List::create(Rcpp::Named("required_mib") = required_mib,
    Rcpp::Named("caps") = caps, Rcpp::Named("allocated_bytes") = allocated,
    Rcpp::Named("max_error") = max_error, Rcpp::Named("invalid_rejected") = invalid);
}
