prague_ground_dataset = function() {
	file.path(
		tools::R_user_dir("skymodelr", "data"),
		"SkyModelDatasetGround.dat"
	)
}

test_that("calculate_sky_radiance splits Prague atmosphere and sun", {
	coef_file = prague_ground_dataset()
	testthat::skip_if_not(file.exists(coef_file))

	args = list(
		phi = c(90, 180),
		theta = c(20, 20),
		lambda_nm = 550,
		altitude = 0,
		elevation = 20,
		visibility = 50,
		albedo = 0.2,
		azimuth = 90,
		number_cores = 1
	)

	L_all = do.call(calculate_sky_radiance, c(args, list(render_mode = "all")))
	L_atm = do.call(calculate_sky_radiance, c(args, list(render_mode = "atmosphere")))
	L_sun = do.call(calculate_sky_radiance, c(args, list(render_mode = "sun")))

	testthat::expect_equal(L_all, L_atm + L_sun, tolerance = 1e-8)
	testthat::expect_true(all(is.finite(L_all)))
	testthat::expect_gt(L_sun[1], 0)
	testthat::expect_gte(L_sun[1], L_sun[2])
})
