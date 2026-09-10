test_that("installed Prague headers support a separate native consumer", {
  skip_on_cran()
  skip_if_not_installed("Rcpp")
  # Compile a client against installed public headers, with no private source
  # includes, model objects, library search paths, or direct C++ DLL symbols.
  client = new.env(parent = globalenv())
  Rcpp::sourceCpp(test_path("prague-api-client.cpp"), env = client)
  missing = tempfile(fileext = ".dat")
  malformed = tempfile(fileext = ".dat")
  writeBin(charToRaw("invalid"), malformed)
  on.exit(unlink(malformed), add = TRUE)
  errors = client$prague_api_errors(missing, malformed)
  expect_equal(errors$version, 3)
  expect_true(all(unlist(errors[-1])))
  filename = tryCatch(
    resolve_prague_coef_file(0, allow_download = FALSE),
    error = function(e) ""
  )
  skip_if(!file.exists(filename), "Prague ground dataset is not installed")
  x = client$prague_api_probe(filename, missing)
  expect_equal(x$batch_error, 0)
  expect_equal(x$thread_error, 0)
  expect_true(all(is.finite(x$radiance)))
  expect_true(all(is.finite(x$transmission)))
  expect_true(all(x$transmission >= 0 & x$transmission <= 1))
  expect_gt(x$intrinsic_sun, x$attenuated_sun)
  expect_gt(x$attenuated_sun, 0)
  expect_equal(x$elevation_max, 90) # Metadata uses degrees; query angles use radians.
  expect_gt(x$channels, 0)
  expect_gt(x$memory, 0)
  expect_true(all(unlist(x[c(
    "failed_initialize",
    "preserved",
    "invalid_direction",
    "invalid_kind",
    "invalid_parameters"
  )])))

  full = tryCatch(
    resolve_prague_coef_file(5000, allow_download = FALSE),
    error = function(e) ""
  )
  skip_if(!file.exists(full), "Prague full-altitude dataset is not installed")
  accelerated = client$prague_api_acceleration(full)
  expect_equal(accelerated$version, 3)
  expect_true(accelerated$invalid_cache)
  expect_true(accelerated$invalid_table)
  expect_equal(accelerated$max_error, 0)
  expect_equal(accelerated$differing_values, 0)
  expect_equal(accelerated$memory[1], accelerated$memory[2])
  expect_gt(accelerated$memory[3], accelerated$memory[1])
  expect_equal(accelerated$memory[3], accelerated$memory[4])
  expect_lte(accelerated$memory[3] - accelerated$memory[1], 512 * 1024^2)

  limits = client$prague_api_table_limits(full)
  expect_gt(limits$required_mib, 0)
  expect_equal(
    limits$allocated_bytes,
    ifelse(limits$caps >= limits$required_mib, limits$required_mib * 1024^2, 0)
  )
  expect_equal(limits$max_error, 0)
  expect_true(limits$invalid_rejected)

  wide_file = tryCatch(
    resolve_prague_coef_file(0, wide_spectrum = TRUE, allow_download = FALSE),
    error = function(e) ""
  )
  skip_if(
    !file.exists(wide_file),
    "Prague extended-spectrum dataset is not installed"
  )
  wide = client$prague_api_table_limits(wide_file)
  expect_gt(wide$required_mib, 512)
  expect_lt(wide$required_mib, 1024)
  expect_equal(
    wide$allocated_bytes,
    ifelse(wide$caps >= wide$required_mib, wide$required_mib * 1024^2, 0)
  )
  expect_equal(wide$max_error, 0)
  expect_true(wide$invalid_rejected)
})
