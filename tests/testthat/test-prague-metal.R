test_that('Prague Metal exports own their data and share immutable snapshots', {
  skip_on_cran()
  skip_if_not_installed('Rcpp')
  filename = tryCatch(
    resolve_prague_coef_file(5000, allow_download = FALSE),
    error = function(e) ''
  )
  skip_if(!file.exists(filename), 'Full-altitude Prague data not installed')
  client = new.env(parent = globalenv())
  Rcpp::sourceCpp(test_path('prague-metal-client.cpp'), env = client)
  result = client$prague_metal_contract(filename)
  expect_true(all(unlist(result[c(
    'shared',
    'bounds',
    'finite',
    'source',
    'parameters'
  )])))
  expect_gt(result$bytes, 0)
  transmission = client$prague_short_transmission(filename)
  expect_true(all(is.finite(transmission)))
  expect_lt(diff(range(transmission)), 1e-5)
})
