angle_diff = function(a, b) {
	diff = (a - b + 180) %% 360 - 180
	abs(diff)
}

column_azimuth = function(col_index, degrees_per_pixel) {
	((col_index - 0.5) * degrees_per_pixel) %% 360
}

row_altitude = function(row_index, degrees_per_pixel) {
	90 - ((row_index - 0.5) * degrees_per_pixel)
}

test_that("sun is correct at equinox", {
	resolution_height = 2000
	resolution_width = 2 * resolution_height
	degrees_per_pixel = 360 / (resolution_width)
	sun_width_pixels = round(0.5 / degrees_per_pixel)
	half_sun_width = sun_width_pixels / 2 + 1
	east_pixel_center = 0.25 * resolution_width
	generate_sky_latlong(
		datetime = as.POSIXct("2025-03-21 06:15:00", tz = "EST"),
		lat = 38.9072,
		lon = -77.0369,
		resolution = resolution_height,
		number_cores = 2,
		hosek = FALSE
	) -> due_east
	testthat::expect_true(
		abs(which.max(apply(due_east[,, 3], 2, sum)) - east_pixel_center) <=
			half_sun_width
	)

	west_pixel_center = 0.75 * resolution_width
	generate_sky_latlong(
		datetime = as.POSIXct("2025-03-21 18:15:00", tz = "EST"),
		lat = 38.9072,
		lon = -77.0369,
		resolution = resolution_height,
		number_cores = 2,
		hosek = FALSE
	) -> due_west
	testthat::expect_true(
		abs(which.max(apply(due_west[,, 3], 2, sum)) - west_pixel_center) <=
			half_sun_width
	)
})

# Additional test cases
# 1) Eclipse at 44.906907, -72.804280, April 8, 2024, approx 3:29 PM EST
# 2) On Tuesday, September 16, 2025 at 04:21:11 UTC the Sun is at its zenith at Latitude: 2 33' North, Longitude: 113 26' East
# 3) On Tuesday, September 16, 2025 at 04:21:11 UTC the Moon is at its zenith at Latitude: 27 08' North, Longitude: 47 16' East

test_that("sun position matches ephemeris during eclipse", {
	resolution_height = 512
	resolution_width = 2 * resolution_height
	degrees_per_pixel = 360 / resolution_width
	date_time = as.POSIXct("2024-04-08 15:29:00", tz = "America/New_York")
	lat = 44.906907
	lon = -72.804280
	sky = generate_sky_latlong(
		datetime = date_time,
		lat = lat,
		lon = lon,
		resolution = resolution_height,
		number_cores = 1
	)
	col_luminance = apply(sky[,, 3], 2, sum)
	brightest_col = which.max(col_luminance)
	sunmoon_data = swe_dirs_topo_moon_sun(date_time, lat = lat, lon = lon)
	sun_dir = sunmoon_data$sun_dir_topo
	expected_az = (180 + atan2(sun_dir[1], sun_dir[2]) * 180 / pi + 360) %% 360
	actual_az = column_azimuth(brightest_col, degrees_per_pixel)
	testthat::expect_lt(
		angle_diff(actual_az, expected_az),
		max(2, degrees_per_pixel * 3)
	)
})

test_that("sun is overhead at Borneo zenith time", {
	resolution_height = 2000
	resolution_width = 2 * resolution_height
	degrees_per_pixel = 360 / resolution_width
	date_time = as.POSIXct("2025-09-16 04:21:11", tz = "UTC")
	lat = 2 + 33 / 60
	lon = 113 + 26 / 60
	sky = generate_sky_latlong(
		datetime = date_time,
		lat = lat,
		lon = lon,
		resolution = resolution_height,
		number_cores = 1
	)
	col_luminance = apply(sky[,, 3], 2, sum)
	row_luminance = apply(sky[,, 3], 1, sum)

	brightest_col = which.max(col_luminance)
	brightest_row = which.max(row_luminance)

	sunmoon_data = swe_dirs_topo_moon_sun(date_time, lat = lat, lon = lon)
	sun_dir = sunmoon_data$sun_dir_topo
	expected_az = (180 + atan2(sun_dir[1], sun_dir[2]) * 180 / pi + 360) %% 360
	actual_az = column_azimuth(brightest_col, degrees_per_pixel)
	testthat::expect_lt(
		angle_diff(actual_az, expected_az),
		max(2, degrees_per_pixel * 3)
	)
	expected_alt = asin(sun_dir[3]) * 180 / pi
	testthat::expect_gt(expected_alt, 85)
	actual_alt = row_altitude(brightest_row, degrees_per_pixel)
	testthat::expect_lt(
		angle_diff(expected_alt, actual_alt),
		max(2, degrees_per_pixel * 3)
	)
})

test_that("moon aligns with zenith ephemeris", {
	resolution_height = 2000
	resolution_width = 2 * resolution_height
	degrees_per_pixel = 360 / resolution_width
	date_time = as.POSIXct("2025-09-16 04:21:11", tz = "UTC")
	lat = 27 + 8 / 60
	lon = 47 + 16 / 60
	moon = generate_moon_latlong(
		datetime = date_time,
		lat = lat,
		lon = lon,
		resolution = resolution_height,
		number_cores = 1,
		atmospheric_scattering = FALSE
	)
	col_luminance = apply(moon[,, 3], 2, sum)
	row_luminance = apply(moon[,, 3], 1, sum)

	brightest_col = which.max(col_luminance)
	brightest_row = which.max(row_luminance)

	sunmoon_data = swe_dirs_topo_moon_sun(date_time, lat = lat, lon = lon)
	moon_dir = sunmoon_data$moon_dir_topo
	#Azimuth doesn't work because it's spread fairly uniformly across the horizontal image
	expected_alt = asin(moon_dir[3]) * 180 / pi
	testthat::expect_gt(expected_alt, 80)
	actual_alt = row_altitude(brightest_row, degrees_per_pixel)
	testthat::expect_lt(
		angle_diff(expected_alt, actual_alt),
		max(2, degrees_per_pixel * 3)
	)
})

test_that("sun radiance does not drop at 90 deg elevation", {
	resolution_height = 512
	lambda_values = seq(360, 720, by = 40)
	sky_899 = makesky_rcpp(
		elevation = 89.9,
		azimuth_deg = 90,
		resolution = resolution_height,
		number_cores = 1,
		model = "hosek",
		render_mode = "sun",
		lambda_nm = lambda_values
	)
	sky_90 = makesky_rcpp(
		elevation = 90,
		azimuth_deg = 90,
		resolution = resolution_height,
		number_cores = 1,
		model = "hosek",
		render_mode = "sun",
		lambda_nm = lambda_values
	)
	L_899 = attr(sky_899, "L_band")
	L_90 = attr(sky_90, "L_band")
	testthat::skip_if(is.null(L_899) || is.null(L_90))

	max_899 = max(L_899)
	max_90 = max(L_90)
	testthat::expect_gt(max_899, 0)
	testthat::expect_gt(max_90, 0)
	testthat::expect_gt(max_90 / max_899, 0.7)

	row_luminance = rowSums(sky_90[,, 2, drop = FALSE])
	max_row = which.max(row_luminance)
	testthat::expect_lte(max_row, 2)
})

test_that("Prague sun radiance does not drop at 90 deg elevation", {
	coef_file = file.path(
		tools::R_user_dir("skymodelr", "data"),
		"SkyModelDatasetGround.dat"
	)
	testthat::skip_if_not(file.exists(coef_file))

	resolution_height = 512
	sky_899 = makesky_rcpp(
		elevation = 89.9,
		azimuth_deg = 90,
		resolution = resolution_height,
		number_cores = 1,
		model = "prague",
		prg_dataset = coef_file,
		altitude = 0,
		visibility = 50,
		albedo = 0.2,
		render_mode = "sun"
	)
	sky_90 = makesky_rcpp(
		elevation = 90,
		azimuth_deg = 90,
		resolution = resolution_height,
		number_cores = 1,
		model = "prague",
		prg_dataset = coef_file,
		altitude = 0,
		visibility = 50,
		albedo = 0.2,
		render_mode = "sun"
	)
	L_899 = attr(sky_899, "L_band")
	L_90 = attr(sky_90, "L_band")
	testthat::skip_if(is.null(L_899) || is.null(L_90))

	max_899 = max(L_899)
	max_90 = max(L_90)
	testthat::expect_gt(max_899, 0)
	testthat::expect_gt(max_90, 0)
	testthat::expect_gt(max_90 / max_899, 0.7)

	row_luminance = rowSums(sky_90[,, 2, drop = FALSE])
	max_row = which.max(row_luminance)
	testthat::expect_lte(max_row, 2)
})
