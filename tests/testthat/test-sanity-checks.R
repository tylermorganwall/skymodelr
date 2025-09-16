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
		numbercores = 2,
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
		numbercores = 2,
		hosek = FALSE
	) -> due_west
	testthat::expect_true(
		abs(which.max(apply(due_west[,, 3], 2, sum)) - west_pixel_center) <=
			half_sun_width
	)
})

# eclipse at 44.906907, -72.804280, April 8, 2024, approx 3:29 PM EST
# On Tuesday, September 16, 2025 at 04:21:11 UTC the Sun is at its zenith at Latitude: 2° 33' North, Longitude: 113° 26' East
# On Tuesday, September 16, 2025 at 04:21:11 UTC the Moon is at its zenith at Latitude: 27° 08' North, Longitude: 47° 16' East
