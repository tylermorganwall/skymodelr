#' Convert B-V color index to blackbody temperature
#'
#' @param bv Numeric vector of B-V values.
#' @return Numeric vector of temperatures in Kelvin.
#' @keywords internal
bv_to_temperature_K = function(bv) {
	bv = pmax(pmin(bv, 2.0), -0.4)
	temp = 4600 * (1 / (0.92 * bv + 1.7) + 1 / (0.92 * bv + 0.62))
	temp[!is.finite(temp)] = NA_real_
	pmax(pmin(temp, 40000), 1000)
}

#' Convert V magnitude to illuminance in lux
#'
#' @param mV Numeric vector of V magnitudes.
#' @param zero_point Optional scalar multiplier for catalog scaling.
#' @return Numeric vector of photometric illuminance in lux.
#' @keywords internal
vmag_to_ev_lux = function(mV, zero_point = 1) {
	zero_point * 10^((-14.18 - mV) / 2.5)
}

#' Compute luminous efficacy for a blackbody SPD
#'
#' @param temperature_K Numeric vector of blackbody temperatures in Kelvin.
#' @param wavelengths_nm Optional wavelength grid in nm.
#' @param V_lambda Optional photopic V(lambda) samples matching the grid.
#' @return Numeric vector of K_eff values in lm/W.
#' @keywords internal
luminous_efficacy_blackbody = function(
	temperature_K,
	wavelengths_nm = NULL,
	V_lambda = NULL
) {
	if (is.null(wavelengths_nm) || is.null(V_lambda)) {
		cie = cie_1931_2deg()
		wavelengths_nm = cie$lambda_nm
		V_lambda = cie$y_bar
	}
	if (length(wavelengths_nm) != length(V_lambda)) {
		stop("wavelengths_nm and V_lambda must have the same length.")
	}
	w = trapezoid_weights(wavelengths_nm)
	c2 = 1.4387769e-2
	blackbody_spd = function(T) {
		lambda_m = wavelengths_nm * 1e-9
		denom = exp(c2 / (lambda_m * T)) - 1
		spd = 1 / (lambda_m^5 * denom)
		spd[!is.finite(spd)] = 0
		spd
	}
	vapply(temperature_K, function(T) {
		if (!is.finite(T) || T <= 0) {
			return(NA_real_)
		}
		spd = blackbody_spd(T)
		num = sum(spd * V_lambda * w)
		den = sum(spd * w)
		if (!is.finite(num) || !is.finite(den) || den <= 0) {
			return(NA_real_)
		}
		683 * num / den
	}, numeric(1))
}

#' Compute radiometric star amplitude and RGB unit vector
#'
#' @param mV Numeric vector of V magnitudes.
#' @param bv Numeric vector of B-V values.
#' @param wavelengths_nm Optional wavelength grid in nm.
#' @param V_lambda Optional photopic V(lambda) samples matching the grid.
#' @return List with E_v (lux), E_e (W/m^2), K_eff (lm/W), and rgb_unit.
#' @keywords internal
star_radiometric_amplitude = function(
	mV,
	bv,
	wavelengths_nm = NULL,
	V_lambda = NULL
) {
	if (length(mV) != length(bv)) {
		stop("mV and bv must have the same length.")
	}
	cie = cie_1931_2deg()
	if (is.null(wavelengths_nm)) {
		wavelengths_nm = cie$lambda_nm
	}
	if (is.null(V_lambda)) {
		V_lambda = cie$y_bar
	}
	if (length(wavelengths_nm) != length(V_lambda)) {
		stop("wavelengths_nm and V_lambda must have the same length.")
	}
	w = trapezoid_weights(wavelengths_nm)
	c2 = 1.4387769e-2
	blackbody_spd = function(T) {
		lambda_m = wavelengths_nm * 1e-9
		denom = exp(c2 / (lambda_m * T)) - 1
		spd = 1 / (lambda_m^5 * denom)
		spd[!is.finite(spd)] = 0
		spd
	}
	out = mapply(function(m, b) {
		T = bv_to_temperature_K(b)
		E_v = vmag_to_ev_lux(m)
		K_eff = luminous_efficacy_blackbody(T, wavelengths_nm, V_lambda)
		E_e = E_v / K_eff
		spd = blackbody_spd(T)
		X = sum(spd * cie$x_bar * w)
		Y = sum(spd * cie$y_bar * w)
		Z = sum(spd * cie$z_bar * w)
		rgb = xyz_to_srgb(c(X, Y, Z))
		rgb[!is.finite(rgb)] = 0
		rgb = pmax(rgb, 0)
		sum_rgb = sum(rgb)
		if (!is.finite(sum_rgb) || sum_rgb <= 0) {
			rgb = c(1, 1, 1) / 3
		} else {
			rgb = rgb / sum_rgb
		}
		list(E_v = E_v, E_e = E_e, K_eff = K_eff, rgb_unit = rgb)
	}, mV, bv, SIMPLIFY = FALSE)
	E_v = vapply(out, function(x) x$E_v, numeric(1))
	E_e = vapply(out, function(x) x$E_e, numeric(1))
	K_eff = vapply(out, function(x) x$K_eff, numeric(1))
	rgb_unit = do.call(rbind, lapply(out, function(x) x$rgb_unit))
	list(E_v = E_v, E_e = E_e, K_eff = K_eff, rgb_unit = rgb_unit)
}
