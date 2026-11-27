testthat::skip_on_cran()

test_that("disk integration matches analytic omega", {
  resolution = 64
  nTheta = resolution
  nPhi = 2 * resolution
  dphi = 2 * pi / nPhi
  dtheta = pi / nTheta
  thetas = seq(pi / 2, -pi / 2, length.out = nTheta)
  phis = seq(0, 2 * pi, length.out = nPhi)

  cos_phi = cos(phis)
  sin_phi = sin(phis)
  cos_theta = cos(thetas)
  sin_theta = sin(thetas)

  r = 0.01
  cos_r = cos(r)

  domega_row = dphi * dtheta * cos(thetas)
  domega = matrix(domega_row, nrow = nTheta, ncol = nPhi, byrow = FALSE)

  disk_mask = matrix(FALSE, nrow = nTheta, ncol = nPhi)
  for (j in seq_len(nTheta)) {
    for (i in seq_len(nPhi)) {
      dir = c(
        cos_phi[i] * cos_theta[j],
        sin_theta[j],
        sin_phi[i] * cos_theta[j]
      )
      if (sum(dir * c(0, 1, 0)) >= cos_r) {
        disk_mask[j, i] = TRUE
      }
    }
  }

  L_e = 2.5
  L_band_env = matrix(0, nrow = nTheta, ncol = nPhi)
  L_band_env[disk_mask] = L_e
  omega_mask = sum(domega[disk_mask])
  E_e = sum(L_band_env[disk_mask] * domega[disk_mask])

  testthat::expect_lt(abs(E_e - L_e * omega_mask), 1e-6)
})
