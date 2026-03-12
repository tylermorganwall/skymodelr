test_that("star env integration is resolution invariant", {
	star = data.frame(
		ra_rad = 0,
		dec_rad = 0,
		v_mag = 0,
		b_v = 0.65,
		r = 1,
		g = 1,
		b = 1
	)
	make_env = function(resolution) {
		make_starfield_rcpp(
			stars = star,
			resolution = resolution,
			lon_deg = 0,
			lat_deg = 0,
			jd = 2451545.0,
			star_width = 1,
			use_rgb = TRUE,
			atmosphere_effects = FALSE,
			upper_hemisphere_only = FALSE,
			number_cores = 1
		)
	}
	integrate_env = function(env) {
		nTheta = dim(env)[1]
		nPhi = dim(env)[2]
		dtheta = pi / nTheta
		dphi = 2 * pi / nPhi
		theta_c = (seq_len(nTheta) - 0.5) * dtheta
		omega_row = dtheta * dphi * sin(theta_c)
		domega = matrix(omega_row, nrow = nTheta, ncol = nPhi, byrow = FALSE)
		L = env[,, 1] + env[,, 2] + env[,, 3]
		sum(L * domega)
	}
	env_low = make_env(64)
	env_high = make_env(128)
	E_low = integrate_env(env_low)
	E_high = integrate_env(env_high)
	rel_err = abs(E_low - E_high) / max(E_low, E_high, 1e-12)
	testthat::expect_lt(rel_err, 1e-4)
})

test_that("star radiometric scale matches moon conversion order", {
	m_moon = -12.73
	moon_lux = 10.76 * 10^(-0.4 * (m_moon + 16.57))
	K_eff = luminous_efficacy_blackbody(5778)
	moon_E_e = lux_to_radiometric_irradiance(moon_lux, K_eff)
	star_info = star_radiometric_amplitude(mV = m_moon, bv = 0.65)
	star_E_e = star_info$E_e[1]
	ratio = star_E_e / moon_E_e
	testthat::expect_true(is.finite(ratio))
	testthat::expect_lt(abs(log(ratio)), log(2.5))
})
