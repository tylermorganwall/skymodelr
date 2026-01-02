#' Normalize a vector to unit length
#'
#' @description Scale a numeric vector so it has unit Euclidean length.
#' @param v Numeric vector to normalise.
#' @keywords internal
normalize = function(v) {
	s = sqrt(sum(v * v))
	if (s == 0) {
		stop("Zero-length vector")
	}
	return(v / s)
}

#' 3D cross product helper
#'
#' @description Compute the 3D cross product of two numeric vectors.
#' @param x Numeric vector of length 3.
#' @param y Numeric vector of length 3.
#' @keywords internal
cross_prod = function(x, y) {
	return(c(
		x[2] * y[3] - x[3] * y[2],
		-(x[1] * y[3] - x[3] * y[1]),
		x[1] * y[2] - x[2] * y[1]
	))
}

#' Convert altitude/azimuth to ENU
#'
#' @description Map altitude/azimuth angles (north-referenced) to a unit
#' East-North-Up direction vector.
#' @param alt_rad Altitude in radians.
#' @param azN_rad Azimuth in radians measured from north.
#' @keywords internal
altaz_to_enu = function(alt_rad, azN_rad) {
	east = sin(azN_rad) * cos(alt_rad)
	north = cos(azN_rad) * cos(alt_rad)
	up = sin(alt_rad)
	normalize(c(east, north, up)) # (E, N, U) matches engine XYZ
}


#' Moon V-band magnitude from phase angle
#'
#' @param phase_deg Phase angle ψ in degrees (0 = full Moon, 180 = new).
#' @param dist_km   Default `384400`. Current geocentric distance (km).
#' @param mean_km   Default `384400`. Mean distance (km).
#'
#' @return Apparent V magnitude `m`.
moon_mag = function(
	phase_deg,
	dist_km = 384400,
	mean_km = 384400
) {
	# Allen-1976 magnitude for a disc at mean distance, no surge (phase_deg in degrees):
	m0 = -12.73 + 0.026 * abs(phase_deg) + 4e-9 * phase_deg^4

	# Geometry: inverse-square scaling with distance
	geom = 5 * log10(dist_km / mean_km)

	# Opposition surge term (Allen’s “extra” brightness near full Moon)
	S = max(1, 1.35 - 0.05 * abs(phase_deg)) # Allen’s empirical fit
	surge_corr = -2.5 * log10(S)

	m = m0 + geom + surge_corr
	return(m)
}


# Krisciunas & Schaefer (1991) / Rozenberg (1966) airmass approximation:
# X(Z) = [cos(Z) + 0.025 * exp(-11 * cos(Z))]^(-1)
# where Z is zenith distance in radians; bounded to ~40 at the horizon.
airmass_rozenberg1966 = function(zenith_deg) {
	z = zenith_deg * pi / 180
	c = cos(z)
	1 / (c + 0.025 * exp(-11 * c))
}

apply_airmass_extinction = function(lux_top_atm, alt_deg, kV = 0.172) {
	# No direct-beam moonlight below (or at) the horizon
	if (!is.finite(lux_top_atm) || !is.finite(alt_deg) || alt_deg <= 0) {
		return(0)
	}

	# Rozenberg is for 0 <= Z <= 90 deg; clamp tiny numerical overshoots
	alt_deg = min(alt_deg, 90)
	zenith_deg = 90 - alt_deg

	X = airmass_rozenberg1966(zenith_deg)
	lux_top_atm * 10^(-0.4 * kV * X)
}

#' Compute topocentric moon/sun directions
#'
#' @description Query Swiss Ephemeris for topocentric sun and moon direction
#' vectors and magnitudes for a given observer location.
#' @param datetime POSIXct observation time (UTC).
#' @param lat Latitude in degrees.
#' @param lon Longitude in degrees.
#' @param elev_m Observer elevation above sea level in metres.
#' @param albedo Default `0.5`. Ground albedo, range 0 to 1.
#' @param turbidity Default `3`. Atmospheric turbidity, range 1.7 to 10
#'   (*Hosek only*).
#' @param visibility Default `50`. Meteorological range (km); *Prague only*.
#' @param hosek Default `TRUE`. `FALSE` selects the Prague model.
#' @param wide_spectrum Default `FALSE`. Use wide-spectrum (55-channel)
#'   coefficients for Prague at sea level only.
#' @keywords internal
swe_dirs_topo_moon_sun = function(
	datetime,
	lat,
	lon,
	elev_m = 0,
	albedo = 0.5,
	turbidity = 3,
	visibility = 50,
	hosek = TRUE,
	wide_spectrum = FALSE
) {
	swephR::swe_set_topo(lon, lat, elev_m)
	attr(datetime, "tzone") = "UTC"
	ts = as.POSIXlt(datetime, tz = "UTC")
	jd_ut = swephR::swe_utc_to_jd(
		ts$year + 1900,
		ts$mon + 1,
		ts$mday,
		ts$hour,
		ts$min,
		ts$sec,
		swephR::SE$GREG_CAL
	)$dret[1]

	flg_eq_topo = swephR::SE$FLG_SWIEPH +
		swephR::SE$FLG_EQUATORIAL +
		swephR::SE$FLG_TOPOCTR
	flg_eq_geo = swephR::SE$FLG_SWIEPH + swephR::SE$FLG_EQUATORIAL
	long_rad = lon * pi / 180
	lat_rad = lat * pi / 180

	enu_to_ecef = matrix(
		c(
			-sin(long_rad),
			cos(long_rad),
			0,
			-sin(lat_rad) * cos(long_rad),
			-sin(lat_rad) * sin(long_rad),
			cos(lat_rad),
			cos(lat_rad) * cos(long_rad),
			cos(lat_rad) * sin(long_rad),
			sin(lat_rad)
		),
		ncol = 3,
		nrow = 3,
		byrow = TRUE
	)
	sun_t = swephR::swe_calc_ut(jd_ut, swephR::SE$SUN, flg_eq_topo)$xx
	moon_t = swephR::swe_calc_ut(jd_ut, swephR::SE$MOON, flg_eq_topo)$xx
	moon_size = swephR::swe_pheno_ut(jd_ut, swephR::SE$MOON, flg_eq_topo)$attr
	sun_size = swephR::swe_pheno_ut(jd_ut, swephR::SE$SUN, flg_eq_topo)$attr

	dist_au = moon_t[3]
	moon_distance_km = dist_au * 149597870.7

	moon_phase = moon_size[1]
	moon_diameter_degrees = moon_size[4]
	moon_brightness_magnitude = moon_size[5]

	m = moon_mag(
		moon_phase,
		dist_km = moon_distance_km,
		mean_km = 384400
	)

	ftcnds_to_lux = 10.76
	moon_brightness_lux_raw = ftcnds_to_lux * 10^(-0.4 * (m + 16.57))

	sun_diameter_degrees = sun_size[4]
	sun_brightness_magnitude = sun_size[5]
	sun_brightness_lux = ftcnds_to_lux *
		10^(-0.4 * (sun_brightness_magnitude + 16.57))

	to_local_horizontal_coordinates = function(ra_dec) {
		aa = swephR::swe_azalt(
			jd_ut,
			swephR::SE$EQU2HOR,
			c(lon, lat, elev_m),
			0,
			0,
			xin = c(ra_dec[1], ra_dec[2], 1)
		)$xaz
		aa
	}
	to_rad = pi / 180

	local_horizon_sun = to_local_horizontal_coordinates(sun_t) * to_rad
	local_horizon_moon = to_local_horizontal_coordinates(moon_t) * to_rad

	azN_s = local_horizon_sun[1]
	alt_s = local_horizon_sun[2]
	azN_mt = local_horizon_moon[1]
	alt_mt = local_horizon_moon[2]

	moon_elev_deg = alt_mt * 180 / pi
	moon_azimuth_deg = (azN_mt * 180 / pi + 180) %% 360
	atten_elev_deg = moon_elev_deg
	if (hosek) {
		atten_elev_deg = max(min(atten_elev_deg, 90), 0)
	} else {
		atten_elev_deg = max(min(atten_elev_deg, 90), -4.2)
	}

	# Then, after you compute moon_elev_deg:
	moon_brightness_lux = apply_airmass_extinction(
		lux_top_atm = moon_brightness_lux_raw,
		alt_deg = moon_elev_deg,
		kV = 0.172
	)

	list(
		sun_dir_topo = altaz_to_enu(alt_s, azN_s),
		sun_ecef_geo = enu_to_ecef %*% altaz_to_enu(alt_s, azN_s),
		moon_dir_topo = altaz_to_enu(alt_mt, azN_mt),
		moon_ecef_geo = enu_to_ecef %*% altaz_to_enu(alt_mt, azN_mt),
		local_up_geo = enu_to_ecef %*% c(0, 0, 1),
		moon_diameter_degrees = moon_diameter_degrees,
		moon_brightness_magnitude = moon_brightness_magnitude,
		moon_brightness_lux = moon_brightness_lux, #lux
		moon_phase = moon_phase,
		sun_diameter_degrees = sun_diameter_degrees,
		sun_brightness_magnitude = sun_brightness_magnitude,
		sun_brightness_lux = sun_brightness_lux #lux
	) |>
		lapply(as.numeric)
}


#' Convert Swiss ephemeris data to a starfield frame
#'
#' @description Build a data frame compatible with `make_starfield_rcpp()` from
#' Swiss Ephemeris results.
#' @param sweph_info List of planetary info objects.
#' @keywords internal
sweph_info_to_stars_df = function(sweph_info) {
	sweph_info_df = do.call(rbind, lapply(sweph_info, as.data.frame))
	planets = rownames(sweph_info_df)
	data.frame(
		planet = planets,
		ra_rad = sweph_info_df$ra_rad,
		dec_rad = sweph_info_df$dec_rad,
		v_mag = sweph_info_df$v_mag,
		r = 1,
		g = 1,
		b = 1
	)
}

#' Compute bright-planet positions via Swiss Ephemeris
#'
#' @description Produce a data frame of planetary apparent positions and
#' magnitudes for the observer.
#' @param datetime POSIXct observation time (UTC).
#' @param lat Latitude in degrees.
#' @param lon Longitude in degrees.
#' @param elev_m Observer elevation above sea level in metres.
#' @keywords internal
swe_dirs_topo_planets_df = function(datetime, lat, lon, elev_m = 0) {
	swephR::swe_set_topo(lon, lat, elev_m)
	attr(datetime, "tzone") = "UTC"
	ts = as.POSIXlt(datetime, tz = "UTC")
	jd_ut = swephR::swe_utc_to_jd(
		ts$year + 1900,
		ts$mon + 1,
		ts$mday,
		ts$hour,
		ts$min,
		ts$sec,
		swephR::SE$GREG_CAL
	)$dret[1]

	flg_eq_topo = swephR::SE$FLG_SWIEPH +
		swephR::SE$FLG_EQUATORIAL +
		swephR::SE$FLG_TOPOCTR
	flg_eq_geo = swephR::SE$FLG_SWIEPH + swephR::SE$FLG_EQUATORIAL
	long_rad = lon * pi / 180
	lat_rad = lat * pi / 180

	enu_to_ecef = matrix(
		c(
			-sin(long_rad),
			cos(long_rad),
			0,
			-sin(lat_rad) * cos(long_rad),
			-sin(lat_rad) * sin(long_rad),
			cos(lat_rad),
			cos(lat_rad) * cos(long_rad),
			cos(lat_rad) * sin(long_rad),
			sin(lat_rad)
		),
		ncol = 3,
		nrow = 3,
		byrow = TRUE
	)

	to_rad = pi / 180

	to_local_horizontal_coordinates = function(ra_dec) {
		aa = swephR::swe_azalt(
			jd_ut,
			swephR::SE$EQU2HOR,
			c(lon, lat, elev_m),
			0,
			0,
			xin = c(ra_dec[1], ra_dec[2], 1)
		)$xaz
		aa
	}
	get_enu_magnitude = function(planet_flag) {
		position_rad = #to_local_horizontal_coordinates(
			swephR::swe_calc_ut(jd_ut, planet_flag, flg_eq_topo)$xx * to_rad
		# ) *
		# to_rad
		# position_enu = altaz_to_enu(position_hori[1], position_hori[2])
		info = swephR::swe_pheno_ut(jd_ut, planet_flag, flg_eq_topo)$attr
		list(ra_rad = position_rad[1], dec_rad = position_rad[2], v_mag = info[5])
	}
	planet_info = list(
		mercury = get_enu_magnitude(swephR::SE$MERCURY),
		venus = get_enu_magnitude(swephR::SE$VENUS),
		mars = get_enu_magnitude(swephR::SE$MARS),
		jupiter = get_enu_magnitude(swephR::SE$JUPITER),
		saturn = get_enu_magnitude(swephR::SE$SATURN),
		uranus = get_enu_magnitude(swephR::SE$URANUS),
		neptune = get_enu_magnitude(swephR::SE$NEPTUNE)
	)
	return(sweph_info_to_stars_df(planet_info))
}


#' Build an orthonormal frame from a Z axis
#'
#' @description Construct a 3×3 orientation matrix with the third column aligned
#' to a supplied direction.
#' @param x Numeric vector giving the desired Z axis.
#' @keywords internal
build_from_z = function(x) {
	zz = normalize(x)
	a = if (abs(zz[1]) > 0.9999999) c(0, 1, 0) else c(1, 0, 0)
	yy = normalize(cross_prod(zz, a))
	xx = cross_prod(yy, zz)
	return(matrix(c(xx, yy, zz), ncol = 3, byrow = FALSE))
}

sun_solid_angle_sr = function(diameter_deg = 0.533) {
	r = (diameter_deg * 0.5) * pi / 180
	2 * pi * (1 - cos(r))
}

# Prague direct-beam illuminance (lux) corresponding to 1 "unit" of Prague sun radiance
prague_unit_direct_lux = function(
	elevation_deg,
	azimuth_deg,
	albedo,
	visibility,
	altitude,
	prg_dataset,
	sun_diameter_deg = 0.533
) {
	# returns Y accumulated from prague_model.sunRadiance(...) over channels :contentReference[oaicite:4]{index=4}
	Y_unit = calculate_sun_brightness_rcpp(
		albedo = albedo,
		elevation = elevation_deg,
		azimuth_deg = azimuth_deg,
		model = "prague",
		prg_dataset = prg_dataset,
		altitude = altitude,
		visibility = visibility
	)
	683 * Y_unit * sun_solid_angle_sr(sun_diameter_deg)
}


#' Generate the moon into a lat-long image array patch
#'
#' @description Produce a rasterised moon texture aligned for composition into
#' the sky dome. The Earth phase is inferred from the Sun-Moon geometry, and earthshine is
#' implemented as an emissive term scaled relative to solar irradiance.
#'
#' @param datetime POSIXct observation time (UTC).
#' @param lat Latitude in degrees.
#' @param lon Longitude in degrees.
#' @param elev_m Observer elevation above sea level in metres.
#' @param width Output width in pixels.
#' @param height Output height in pixels.
#' @param earthshine Default `TRUE`. If `FALSE`, skip the earthshine emissive term.
#' @param verbose Default `FALSE`. If `TRUE`, print earthshine diagnostics.
#' @keywords internal
generate_moon_image_latlong = function(
	datetime,
	lat,
	lon,
	elev_m = 0,
	width = 400,
	height = 400,
	earthshine = TRUE,
	verbose = FALSE
) {
	moon_sun_data = swe_dirs_topo_moon_sun(datetime, lat, lon, elev_m)
	dir_moon = normalize(moon_sun_data$moon_ecef_geo)
	dir_sun = -normalize(moon_sun_data$sun_ecef_geo)
	local_up = normalize(moon_sun_data$local_up_geo)
	moon_brightness_lux = moon_sun_data$moon_brightness_lux

	tex = system.file(
		"textures",
		"lroc_color_poles_1k.jpg",
		package = "skymodelr"
	)

	clamp = function(x, a, b) {
		pmax(a, pmin(b, x))
	}
	cot = function(x) {
		1 / tan(x)
	}

	cos_phi = clamp(sum(dir_moon * dir_sun), -1, 1)
	phi = acos(cos_phi)

	earth_phase = pi - phi
	earth_phase = clamp(earth_phase, 1e-6, pi - 1e-6)

	# Paper-based earthshine irradiance model:
	# E_em = 0.19 * 0.5 * [1 - sin(phase/2) * tan(phase/2) * ln(cot(phase/4))]
	E_em = 0.19 *
		0.5 *
		(1 -
			sin(earth_phase / 2) *
				tan(earth_phase / 2) *
				log(cot(earth_phase / 4)))

	if (!is.finite(E_em) || E_em < 0) {
		E_em = 0
	}

	# Convert to unitless emission scale relative to solar irradiance at 1 AU.
	# sunlight ~ 1.3e3 W/m^2.
	solar_irradiance = 1.3e3

	earth_emissive_intensity = 0
	if (earthshine) {
		earth_emissive_intensity = E_em / solar_irradiance
	}

	# Build a base mesh transform shared by both passes
	build_moon_mesh = function(material) {
		rayvertex::sphere_mesh(material = material, radius = 0.5) |>
			rayvertex::subdivide_mesh() |>
			rayvertex::rotate_mesh(c(0, -90, 6.68), order_rotation = c(3, 1, 2)) |>
			rayvertex::transform_mesh(rayvertex::lookat_transform(
				pos = -dir_moon,
				look = c(0, 0, 0),
				up = c(0, 0, 1)
			)) |>
			rayvertex::translate_mesh(dir_moon * 10)
	}

	# Pass 1: sun-only material (no emission) for magnitude calibration
	moon_material_sun = rayvertex::material_list(
		texture_location = tex,
		emissive_texture_location = tex,
		diffuse = c(1, 1, 1),
		emission = c(1, 1, 1),
		diffuse_intensity = 1,
		emission_intensity = 0,
		sigma = 10
	)

	moon_mesh_sun = build_moon_mesh(moon_material_sun)

	moon_image_sun = rayvertex::rasterize_scene(
		moon_mesh_sun,
		lookfrom = c(0, 0, 0),
		lookat = dir_moon * 10,
		camera_up = local_up,
		fov = 0,
		light_info = rayvertex::directional_light(
			direction = -dir_sun,
			intensity = 1
		),
		transparent_background = TRUE,
		plot = FALSE,
		width = width,
		height = height
	)

	# Pass 2: combined sun + earthshine using emissive intensity
	moon_material = rayvertex::material_list(
		texture_location = tex,
		emissive_texture_location = tex,
		diffuse = c(1, 1, 1),
		emission = c(1, 1, 1),
		diffuse_intensity = 1,
		emission_intensity = earth_emissive_intensity,
		sigma = 10
	)

	moon_mesh = build_moon_mesh(moon_material)

	moon_image = rayvertex::rasterize_scene(
		moon_mesh,
		lookfrom = c(0, 0, 0),
		lookat = dir_moon * 10,
		camera_up = local_up,
		fov = 0,
		light_info = rayvertex::directional_light(
			direction = -dir_sun,
			intensity = 1
		),
		transparent_background = TRUE,
		plot = FALSE,
		width = width,
		height = height
	)

	if (verbose) {
		message(sprintf(
			"phi=%0.4f, earth_phase=%0.4f, E_em=%0.4g W/m^2, emission_intensity=%0.4g",
			phi,
			earth_phase,
			E_em,
			earth_emissive_intensity
		))
	}

	moon_image[,, 4][moon_image[,, 4] > 1] = 1

	list(
		moon_luminance_array = moon_image,
		moon_angular_diameter_deg = moon_sun_data$moon_diameter_degrees
	)
}
