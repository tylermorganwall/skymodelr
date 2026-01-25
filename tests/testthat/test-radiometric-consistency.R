test_that("constant env is resolution invariant", {
	testthat::skip_if_not_installed("rayrender")
	testthat::skip_if_not_installed("rayimage")

	make_env = function(resolution) {
		env = array(1, dim = c(resolution, resolution * 2, 3))
		env_path = tempfile(fileext = ".exr")
		rayimage::ray_write_image(env, env_path, clamp = FALSE)
		env_path
	}

	render_plane = function(env_path) {
		scene = rayrender::xz_rect(
			x = 0,
			xwidth = 2,
			z = 0,
			zwidth = 2,
			y = 0,
			material = rayrender::diffuse(color = "#ffffff")
		)
		set.seed(1)
		render = rayrender::render_scene(
			scene,
			width = 32,
			height = 32,
			fov = 0,
			ortho_dimensions = c(2, 2),
			lookfrom = c(0, 1, 0),
			lookat = c(0, 0, 0),
			camera_up = c(0, 0, 1),
			samples = 64,
			ambient_light = FALSE,
			backgroundhigh = "#000000",
			backgroundlow = "#000000",
			environment_light = env_path,
			intensity_env = 1,
			tonemap = "raw",
			bloom = FALSE,
			preview = FALSE,
			interactive = FALSE,
			plot_scene = FALSE,
			progress = FALSE
		)
		mean(render[,, 1:3])
	}

	env_low = make_env(32)
	env_high = make_env(256)

	mean_low = render_plane(env_low)
	mean_high = render_plane(env_high)

	testthat::expect_lt(abs(mean_low - mean_high), 0.05)
})

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

test_that("Prague env irradiance matches rayrender plane", {
	testthat::skip_if_not_installed("rayrender")
	testthat::skip_if_not_installed("rayimage")

	coef_file = file.path(
		tools::R_user_dir("skymodelr", "data"),
		"SkyModelDatasetGround.dat"
	)
	testthat::skip_if_not(file.exists(coef_file))

	sky = generate_sky(
		resolution = 64,
		hosek = FALSE,
		altitude = 0,
		visibility = 50,
		albedo = 0.2,
		elevation = 20,
		azimuth = 180,
		render_mode = "atmosphere",
		numbercores = 1
	)
	L_band = attr(sky, "L_band")
	testthat::skip_if(is.null(L_band))

	nTheta = dim(sky)[1]
	nPhi = dim(sky)[2]
	dphi = 2 * pi / nPhi
	dtheta = pi / nTheta
	thetas = seq(pi / 2, -pi / 2, length.out = nTheta)
	domega_row = dphi * dtheta * cos(thetas)
	domega = matrix(domega_row, nrow = nTheta, ncol = nPhi, byrow = FALSE)
	cos_incidence = pmax(0, sin(thetas))
	cos_incidence = matrix(cos_incidence, nrow = nTheta, ncol = nPhi, byrow = FALSE)
	E_env = sum(L_band * cos_incidence * domega)

	env_path = tempfile(fileext = ".exr")
	rayimage::ray_write_image(sky[,, 1:3], env_path, clamp = FALSE)
	scene = rayrender::xz_rect(
		x = 0,
		xwidth = 2,
		z = 0,
		zwidth = 2,
		y = 0,
		material = rayrender::diffuse(color = "#ffffff")
	)
	set.seed(1)
	render = rayrender::render_scene(
		scene,
		width = 32,
		height = 32,
		fov = 0,
		ortho_dimensions = c(2, 2),
		lookfrom = c(0, 1, 0),
		lookat = c(0, 0, 0),
		camera_up = c(0, 0, 1),
		samples = 128,
		ambient_light = FALSE,
		backgroundhigh = "#000000",
		backgroundlow = "#000000",
		environment_light = env_path,
		intensity_env = 1,
		tonemap = "raw",
		bloom = FALSE,
		preview = FALSE,
		interactive = FALSE,
		plot_scene = FALSE,
		progress = FALSE
	)
	mean_radiance = mean(render[,, 1:3])
	E_render = pi * mean_radiance
	rel_err = abs(E_env - E_render) / max(E_env, 1e-6)
	testthat::expect_lt(rel_err, 0.15)
})
