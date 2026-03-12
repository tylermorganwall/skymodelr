sun_solid_angle_sr = function(diameter_deg = 0.533) {
	r = (diameter_deg * 0.5) * pi / 180
	2 * pi * (1 - cos(r))
}

#' Generate the atmosphere with the moon
#'
#' @description Note that this is just a scaled version of `generate_sky()`, scaled down by the luminance
#' of the moon as compared to the sun. This function takes the phase of the moon into account,
#' along with the increase in luminosity around a full moon (known as opposition surge). Moonlight
#' attenuation uses the Rozenberg/Krisciunas-Schaefer airmass approximation with configurable `moon_extinction_kV`.
#'
#' @param datetime           POSIX-compatible date-time.
#' @param lat                Observer latitude (degrees N).
#' @param lon                Observer longitude (degrees E; west < 0).
#' @param filename           Default `NA`. Path to the image file to write.
#' @param albedo             Default `0.5`. Ground albedo, range 0 to 1.
#' @param turbidity          Default `3`. Atmospheric turbidity, range
#'   1.7 to 10 (*Hosek only*).
#' @param altitude           Default `0`. Observer altitude (m), range
#'   0 to 15000 (*Prague only*).
#' @param resolution         Default `2048`. Image height in pixels (width = 2 * height).
#' @param number_cores        Default `1`. CPU threads to use.
#' @param moon_atmosphere    Default `FALSE`. If `TRUE`, this generates atmospheric scattering from light from the moon.
#' @param earthshine         Default `TRUE`. If `FALSE`, disable earthshine contribution on the moon.
#' @param earthshine_albedo  Default `0.19`. Effective Earth albedo term used
#'   in the earthshine irradiance approximation.
#' @param solar_irradiance_w_m2 Default `1300`. Reference solar irradiance at
#'   1 AU (W/m^2) used to normalize earthshine emission intensity.
#' @param moon_extinction_kV Default `0.172`. V-band atmospheric extinction
#'   coefficient for direct moonlight attenuation.
#' @param hosek              Default `TRUE`. `FALSE` selects the Prague model.
#' @param wide_spectrum      Default `FALSE`. 55-channel Prague coefficients (altitude = 0m only).
#' @param visibility         Default `50`. Meteorological range (km); *Prague only*.
#' @param moon_texture_width Default `801`. Internal moon-texture render width.
#' @param moon_texture_height Default `801`. Internal moon-texture render height.
#' @param verbose            Default `FALSE`. Whether to print progress bars/diagnostic info.
#'
#' @return Either the image array, or the array is invisibly returned if a file
#'   is written. The array has dimensions `(resolution, 2 * resolution, 4)`.
#' @note Writing to non-EXR formats will introduce precision loss because HDR
#'   data are quantised to the destination format, and low dynamic range outputs like PNG
#'   and JPEG files will not represent the true luminosity values encoded in the array.
#' @export
#'
#' @examples
#' # Moonlit sky (Hosek), mid-evening in DC
#' if(run_documentation()) {
#' generate_moon_latlong(
#'   datetime   = as.POSIXct("2025-09-05 19:30:00",tz="America/New_York"),
#'   lat        = 38.9072,
#'   lon        = -77.0369,
#'   resolution = 400,
#'   turbidity  = 3,
#'   verbose    = TRUE
#' ) |>
#'   rayimage::render_exposure(15) |>
#'   rayimage::plot_image()
#' }
generate_moon_latlong = function(
	datetime,
	lat,
	lon,
	filename = NA,
	albedo = 0.5,
	turbidity = 3,
	altitude = 0,
	resolution = 2048,
	number_cores = 1,
	moon_atmosphere = FALSE,
	earthshine = TRUE,
	earthshine_albedo = 0.19,
	solar_irradiance_w_m2 = 1300,
	moon_extinction_kV = 0.172,
	hosek = TRUE,
	wide_spectrum = FALSE,
	visibility = 50,
	moon_texture_width = 801,
	moon_texture_height = 801,
	verbose = FALSE
) {
	moon_sun_data = swe_dirs_topo_moon_sun(
		datetime = datetime,
		lat = lat,
		lon = lon,
		elev_m = altitude,
		moon_extinction_kV = moon_extinction_kV
	)
	moon_dir = moon_sun_data$moon_dir_topo
	moon_elevation = asin(moon_dir[3]) * 180 / pi
	moon_azimuth = (180 + atan2(moon_dir[1], moon_dir[2]) * 180 / pi + 360) %% 360
	if (hosek) {
		if (moon_elevation < 0.0) {
			if (verbose) {
				message(
					"Drawing black image as Hosek model does not produce valid output for elevation < 0."
				)
			}
			black_sky = array(0, dim = c(resolution, resolution * 2, 4))
			black_sky[,, 4] = 1
			if (!is.na(filename)) {
				warn_precision_loss(filename)
				rayimage::ray_write_image(black_sky, filename)
				return(invisible(black_sky))
			} else {
				return(black_sky)
			}
		}
	} else {
		if (moon_elevation < -4.2) {
			if (verbose) {
				message(
					"Drawing black image as Prague model does not produce valid output for elevation < -4.2."
				)
			}
			black_sky = array(0, dim = c(resolution, resolution * 2, 4))
			black_sky[,, 4] = 1
			if (!is.na(filename)) {
				warn_precision_loss(filename)
				rayimage::ray_write_image(black_sky, filename)
				return(invisible(black_sky))
			} else {
				return(black_sky)
			}
		}
	}
	# Sample sun colour so the moon inherits the same atmospheric tint.
	sun_rgb_ratio = c(1, 1, 1)
	tint_min_elevation = if (hosek) 0 else -4.2
	if (is.finite(moon_elevation) && moon_elevation >= tint_min_elevation) {
		sample_resolution = max(256, min(1024, resolution))
		sun_sample = tryCatch(
			{
				generate_sky(
					albedo = albedo,
					turbidity = turbidity,
					altitude = altitude,
					elevation = moon_elevation,
					visibility = visibility,
					azimuth = moon_azimuth,
					resolution = sample_resolution,
					number_cores = number_cores,
					hosek = hosek,
					wide_spectrum = wide_spectrum,
					render_mode = "sun",
					verbose = FALSE
				)
			},
			error = function(e) {
				NULL
			}
		)
		if (!is.null(sun_sample)) {
			sun_lum = sun_sample[,, 1] + sun_sample[,, 2] + sun_sample[,, 3]
			peak_idx = which.max(sun_lum)
			sun_rgb = c(
				sun_sample[,, 1][peak_idx],
				sun_sample[,, 2][peak_idx],
				sun_sample[,, 3][peak_idx]
			)
			sun_mean = mean(sun_rgb)
			if (all(is.finite(sun_rgb)) && sun_mean > 0) {
				sun_rgb_ratio = sun_rgb / sun_mean
			}
		}
	}

	moon_lux = moon_sun_data$moon_brightness_lux
	moon_phase = moon_sun_data$moon_phase
	spd_type = "BB5778"
	rgb_unit = compute_spd_rgb_unit(spd_type)
	K_eff = compute_K_eff(spd_type)
	moon_irradiance = lux_to_radiometric_irradiance(moon_lux, K_eff)
	if (!is.finite(moon_irradiance)) {
		moon_irradiance = 0
	}

	moon_info_list = generate_moon_image_latlong(
		datetime = datetime,
		lat = lat,
		lon = lon,
		elev_m = altitude,
		width = moon_texture_width,
		height = moon_texture_height,
		earthshine = earthshine,
		earthshine_albedo = earthshine_albedo,
		solar_irradiance_w_m2 = solar_irradiance_w_m2,
		moon_extinction_kV = moon_extinction_kV,
		verbose = verbose
	)
	moon_luminance_array = moon_info_list$moon_luminance_array
	moon_angular_diameter_deg = moon_info_list$moon_angular_diameter_deg
	moon_angular_diameter_rad = moon_angular_diameter_deg * pi / 180

	if (verbose) {
		message(sprintf(
			"Moon: %0.1f elevation, %0.1f azimuth, %0.3f phase, %f lux",
			moon_elevation,
			moon_azimuth,
			moon_phase,
			moon_lux
		))
	}
	if (moon_atmosphere) {
		sea_level = altitude == 0
		filesize = ""
		if (sea_level & !wide_spectrum) {
			filesize = "107MB"
		} else if (sea_level & wide_spectrum) {
			filesize = "574MB"
		} else if (!sea_level & !wide_spectrum) {
			filesize = "2.4GB"
		}
		check_coef_file = function(filename) {
			coef_file = file.path(tools::R_user_dir("skymodelr", "data"), filename)
			if (!file.exists(coef_file)) {
				response = readline(
					prompt = sprintf(
						" Coefficient file for this setting not yet present: this is a large file (%s), download? [y/n] ",
						filesize
					)
				)
				if (response == "y") {
					download_sky_data(sea_level, wide_spectrum)
				} else if (response == "n") {
					return("")
				} else {
					stop("Input not recognized.")
				}
			}
			return(coef_file)
		}
		coef_file = ""
		stopifnot(all(altitude >= 0 & altitude <= 15000))
		if (!hosek) {
			if (!wide_spectrum) {
				if (sea_level) {
					coef_file = check_coef_file("SkyModelDatasetGround.dat")
				} else {
					coef_file = check_coef_file("SkyModelDataset.dat")
				}
			} else {
				if (sea_level) {
					coef_file = check_coef_file("PragueSkyModelDatasetGroundInfra.dat")
				} else {
					stop("`wide_spectrum = TRUE` is only valid when `altitude == 0`.")
				}
			}
			if (coef_file == "") {
				stop("No coefficient file downloaded for this set of inputs.")
			}
		}
		moon_array = generate_sky(
			albedo = albedo,
			turbidity = turbidity,
			altitude = altitude,
			elevation = moon_elevation,
			azimuth = moon_azimuth,
			resolution = resolution,
			number_cores = number_cores,
			hosek = hosek,
			wide_spectrum = wide_spectrum,
			visibility = visibility,
			verbose = verbose,
			render_mode = "atmosphere"
		)
		lambda_values = seq(360, 720, by = 40)
		model = if (hosek) "hosek" else "prague"
		E_unit = calculate_sun_radiance_band_rcpp(
			albedo = albedo,
			turbidity = turbidity,
			elevation = moon_elevation,
			azimuth_deg = moon_azimuth,
			model = model,
			prg_dataset = coef_file,
			altitude = altitude,
			visibility = visibility,
			lambda_nm = if (hosek) lambda_values else NULL
		) *
			sun_solid_angle_sr(0.533)
		scale_sun_to_moon = if (is.finite(E_unit) && E_unit > 0) {
			moon_irradiance / E_unit
		} else {
			0
		}
		moon_luminance_array[,, 1:3] = sweep(
			moon_luminance_array[,, 1:3] * scale_sun_to_moon,
			3,
			sun_rgb_ratio,
			"*"
		)
		# Apply radiometric scaling and atmospheric tint
		moon_array[,, 1:3] = sweep(
			moon_array[,, 1:3] * scale_sun_to_moon,
			3,
			sun_rgb_ratio,
			"*"
		)
	} else {
		moon_array = array(
			0,
			dim = c(resolution, resolution * 2, 4)
		)
		moon_array[,, 4] = 1
	}

	moon_dir_vec = c(
		cospi(moon_azimuth / 180) * cospi(moon_elevation / 180),
		sinpi(moon_elevation / 180),
		sinpi(moon_azimuth / 180) * cospi(moon_elevation / 180)
	)

	moon_diameter_pixels = resolution / 180 * moon_angular_diameter_deg
	resize_moon_dim = ceiling(moon_diameter_pixels)
	if (resize_moon_dim %% 2 == 0) {
		resize_moon_dim = resize_moon_dim + 1
	}
	resized_moon_luminance_array = rayimage::render_resized(
		moon_luminance_array,
		dims = c(resize_moon_dim, resize_moon_dim),
		method = "box"
	)
	patch_mask = resized_moon_luminance_array[,, 4] > 0
	if (!any(patch_mask)) {
		patch_mask[,] = TRUE
	}
	alpha_mask = resized_moon_luminance_array[,, 4] > 0
	if (any(alpha_mask)) {
		alpha_inds = which(alpha_mask, arr.ind = TRUE)
		min_row = min(alpha_inds[, 1])
		max_row = max(alpha_inds[, 1])
		min_col = min(alpha_inds[, 2])
		max_col = max(alpha_inds[, 2])
		center_x = (min_col + max_col) / 2
		center_y = (min_row + max_row) / 2
		radius_x = max(1, (max_col - min_col + 1) / 2)
		radius_y = max(1, (max_row - min_row + 1) / 2)
	} else {
		half_dim = (resize_moon_dim - 1) / 2
		center_x = half_dim + 1
		center_y = half_dim + 1
		radius_x = half_dim
		radius_y = half_dim
	}
	rgb_target = rgb_unit * sun_rgb_ratio
	if (sum(rgb_target) > 0) {
		rgb_target = rgb_target / sum(rgb_target)
	} else {
		rgb_target = rgb_unit
	}
	mean_rgb = c(
		mean(resized_moon_luminance_array[,, 1][patch_mask], na.rm = TRUE),
		mean(resized_moon_luminance_array[,, 2][patch_mask], na.rm = TRUE),
		mean(resized_moon_luminance_array[,, 3][patch_mask], na.rm = TRUE)
	)
	mean_rgb[!is.finite(mean_rgb)] = 0
	rgb_scale = ifelse(mean_rgb > 0, rgb_target / mean_rgb, 0)
	resized_moon_luminance_array[,, 1:3] = sweep(
		resized_moon_luminance_array[,, 1:3],
		3,
		rgb_scale,
		"*"
	)
	mean_sum = mean(
		rowSums(resized_moon_luminance_array[,, 1:3])[patch_mask],
		na.rm = TRUE
	)
	if (is.finite(mean_sum) && mean_sum > 0) {
		resized_moon_luminance_array[,, 1:3] =
			resized_moon_luminance_array[,, 1:3] / mean_sum
	}

	nTheta = resolution
	nPhi = 2 * resolution
	disk_mask = matrix(FALSE, nrow = nTheta, ncol = nPhi)

	r = moon_angular_diameter_rad / 2
	cos_r = cos(r)

	phi_c = moon_elevation * pi / 180
	j_min = max(1, floor(((pi / 2) - (phi_c + r)) / pi * (nTheta - 1)) + 1)
	j_max = min(
		nTheta,
		ceiling(((pi / 2) - (phi_c - r)) / pi * (nTheta - 1)) + 1
	)

	thetas = seq(pi / 2, -pi / 2, length.out = nTheta)
	phis = seq(0, 2 * pi, length.out = nPhi)

	cos_phi = cos(phis)
	sin_phi = sin(phis)
	cos_theta = cos(thetas)
	sin_theta = sin(thetas)

	# Tangent frame at moon center (z = moon_dir_vec)
	lam_c = moon_azimuth * pi / 180
	z_hat = moon_dir_vec
	e_hat = c(-sin(lam_c), 0, cos(lam_c))
	e_hat = e_hat / sqrt(sum(e_hat * e_hat)) # local east
	n_hat = c(
		z_hat[2] * e_hat[3] - z_hat[3] * e_hat[2],
		z_hat[3] * e_hat[1] - z_hat[1] * e_hat[3],
		z_hat[1] * e_hat[2] - z_hat[2] * e_hat[1]
	)
	n_hat = n_hat / sqrt(sum(n_hat * n_hat)) # local north

	tan_r = tan(r) # r = moon_angular_diameter_rad/2

	for (j in seq(j_min, j_max)) {
		for (i in seq_len(nPhi)) {
			sample_dir_vec = c(
				cos_phi[i] * cos_theta[j],
				sin_theta[j],
				sin_phi[i] * cos_theta[j]
			)
			if (sample_dir_vec[2] < 0) {
				next
			}
			dot = sum(moon_dir_vec * sample_dir_vec)
			if (dot < cos_r) {
				next
			}

			disk_mask[j, i] = TRUE

			# Camera (tangent) components of the sample direction
			vx = sum(sample_dir_vec * e_hat)
			vy = sum(sample_dir_vec * n_hat)
			vz = sum(sample_dir_vec * z_hat) # == dot

			# Inverse rectilinear mapping (square FOV 2*r):
			u_px = center_x + radius_x * (vx / (tan_r * vz))
			v_px = center_y + radius_y * (vy / (tan_r * vz))

			u_i = rayrender:::clamp(round(u_px), 1, resize_moon_dim)
			v_i = rayrender:::clamp(round(v_px), 1, resize_moon_dim)

			moon_array[j, i, 1] = resized_moon_luminance_array[v_i, u_i, 1]
			moon_array[j, i, 2] = resized_moon_luminance_array[v_i, u_i, 2]
			moon_array[j, i, 3] = resized_moon_luminance_array[v_i, u_i, 3]
		}
	}

	dphi = 2 * pi / nPhi
	dtheta = pi / nTheta
	thetas = seq(pi / 2, -pi / 2, length.out = nTheta) # you already have this
	domega_row = dphi * dtheta * cos(thetas) # cos(latitude) = sin(colatitude)
	domega = matrix(domega_row, nrow = nTheta, ncol = nPhi, byrow = FALSE)
	omega_disk = sum(domega[disk_mask])
	if (is.finite(omega_disk) && omega_disk > 0) {
		texture_sum = moon_array[,, 1] + moon_array[,, 2] + moon_array[,, 3]
		texture_mean_env = sum(
			texture_sum[disk_mask] * domega[disk_mask],
			na.rm = TRUE
		) /
			omega_disk
		if (is.finite(texture_mean_env) && texture_mean_env > 0) {
			moon_array[,, 1][disk_mask] = moon_array[,, 1][disk_mask] /
				texture_mean_env
			moon_array[,, 2][disk_mask] = moon_array[,, 2][disk_mask] /
				texture_mean_env
			moon_array[,, 3][disk_mask] = moon_array[,, 3][disk_mask] /
				texture_mean_env
		}
	}
	omega_moon = moon_solid_angle_sr(moon_angular_diameter_deg)
	omega_used = if (is.finite(omega_disk) && omega_disk > 0) {
		omega_disk
	} else {
		omega_moon
	}
	L_e = if (is.finite(omega_used) && omega_used > 0) {
		moon_irradiance / omega_used
	} else {
		0
	}
	if (L_e > 0) {
		moon_array[,, 1][disk_mask] = moon_array[,, 1][disk_mask] * L_e
		moon_array[,, 2][disk_mask] = moon_array[,, 2][disk_mask] * L_e
		moon_array[,, 3][disk_mask] = moon_array[,, 3][disk_mask] * L_e
	}
	moon_band = attr(moon_array, "L_band")
	if (is.null(moon_band)) {
		moon_band = moon_array[,, 1] + moon_array[,, 2] + moon_array[,, 3]
	} else if (any(disk_mask)) {
		moon_band[disk_mask] = moon_array[,, 1][disk_mask] +
			moon_array[,, 2][disk_mask] +
			moon_array[,, 3][disk_mask]
	}
	moon_array = rayimage::ray_read_image(
		moon_array,
		assume_white = "D65",
		assume_colorspace = rayimage::CS_SRGB
	)
	attr(moon_array, "L_band") = moon_band
	if (!is.na(filename)) {
		warn_precision_loss(filename)
		rayimage::ray_write_image(moon_array, filename, clamp = FALSE)
		return(invisible(moon_array))
	} else {
		return(moon_array)
	}
}
