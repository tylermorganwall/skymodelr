#' Moon disk solid angle in steradians
#'
#' @param diameter_deg Angular diameter in degrees.
#' @return Solid angle in steradians.
#' @keywords internal
moon_solid_angle_sr = function(diameter_deg = 0.533) {
	r = (diameter_deg * 0.5) * pi / 180
	2 * pi * (1 - cos(r))
}

cie_1931_2deg = local({
	cache = NULL
	function() {
		if (is.null(cache)) {
			cache = cie_1931_2deg_rcpp()
		}
		cache
	}
})

trapezoid_weights = function(lambda_nm) {
	n = length(lambda_nm)
	if (n == 0) {
		return(numeric(0))
	}
	if (n == 1) {
		return(1)
	}
	w = numeric(n)
	w[1] = (lambda_nm[2] - lambda_nm[1]) / 2
	for (i in 2:(n - 1)) {
		w[i] = (lambda_nm[i + 1] - lambda_nm[i - 1]) / 2
	}
	w[n] = (lambda_nm[n] - lambda_nm[n - 1]) / 2
	w
}

xyz_to_srgb = function(XYZ) {
	M = matrix(
		c(
			3.2404542,
			-1.5371385,
			-0.4985314,
			-0.9692660,
			1.8760108,
			0.0415560,
			0.0556434,
			-0.2040259,
			1.0572252
		),
		nrow = 3,
		byrow = TRUE
	)
	as.numeric(M %*% XYZ)
}

#' Compute a unit RGB vector for a chosen SPD
#'
#' @param spd_type SPD shape. Either `"D65"` or `"BB5778"`.
#' @return Numeric vector length 3, normalized so sum equals 1.
#' @keywords internal
compute_spd_rgb_unit = function(spd_type = c("D65", "BB5778")) {
	spd_type = match.arg(spd_type)
	cie = cie_1931_2deg()
	lambda_nm = cie$lambda_nm
	spd = spd_values_rcpp(lambda_nm, spd_type)
	w = trapezoid_weights(lambda_nm)
	X = sum(spd * cie$x_bar * w)
	Y = sum(spd * cie$y_bar * w)
	Z = sum(spd * cie$z_bar * w)
	rgb = xyz_to_srgb(c(X, Y, Z))
	rgb[!is.finite(rgb)] = 0
	rgb = pmax(rgb, 0)
	total = sum(rgb)
	if (!is.finite(total) || total <= 0) {
		return(c(0, 0, 0))
	}
	rgb / total
}

#' Compute effective luminous efficacy for a chosen SPD
#'
#' @param spd_type SPD shape. Either `"D65"` or `"BB5778"`.
#' @param wavelength_grid Optional numeric vector of wavelengths (nm).
#' @param V_lambda Optional numeric vector of photopic V(lambda) samples.
#' @return Effective luminous efficacy in lm/W.
#' @keywords internal
compute_K_eff = function(
	spd_type = c("D65", "BB5778"),
	wavelength_grid = NULL,
	V_lambda = NULL
) {
	spd_type = match.arg(spd_type)
	if (is.null(wavelength_grid) || is.null(V_lambda)) {
		cie = cie_1931_2deg()
		wavelength_grid = cie$lambda_nm
		V_lambda = cie$y_bar
	}
	if (length(wavelength_grid) != length(V_lambda)) {
		stop("wavelength_grid and V_lambda must have the same length.")
	}
	spd = spd_values_rcpp(wavelength_grid, spd_type)
	w = trapezoid_weights(wavelength_grid)
	num = sum(spd * V_lambda * w)
	den = sum(spd * w)
	if (!is.finite(num) || !is.finite(den) || den <= 0) {
		return(NA_real_)
	}
	683 * num / den
}

#' Convert photometric illuminance to radiometric irradiance
#'
#' @param E_v_lux Illuminance in lux.
#' @param K_eff Effective luminous efficacy in lm/W.
#' @return Irradiance in W/m^2.
#' @keywords internal
lux_to_radiometric_irradiance = function(E_v_lux, K_eff) {
	E_v_lux / K_eff
}
