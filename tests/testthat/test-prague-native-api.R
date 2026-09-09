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
  expect_equal(errors$version, 1)
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
})
