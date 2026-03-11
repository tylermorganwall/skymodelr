prep_cmf = function(lambda_nm, cmf_file = "data-raw/CIE_xyz_1931_2deg.csv") {
	stopifnot(is.numeric(lambda_nm), length(lambda_nm) >= 2)
	cmf_src = read.csv(cmf_file, header = FALSE)
	colnames(cmf_src) = c("lambda", "x_bar", "y_bar", "z_bar")

	f = function(col) {
		approx(
			cmf_src$lambda,
			cmf_src[[col]],
			xout = lambda_nm,
			rule = 2,
			ties = "ordered"
		)$y
	}
	x_bar = f("x_bar")
	y_bar = f("y_bar")
	z_bar = f("z_bar")

	# trapezoidal weights
	d = diff(lambda_nm)
	w = numeric(length(lambda_nm))
	if (length(lambda_nm) == 2) {
		w[] = d[1] / 2
	} else {
		w[1] = d[1] / 2
		w[2:(length(lambda_nm) - 1)] = (d[-1] + d[-length(d)]) / 2
		w[length(lambda_nm)] = d[length(d)] / 2
	}

	list(lambda = lambda_nm, cmf_xyz = cbind(x_bar, y_bar, z_bar), weights = w)
}


plot(
	x = seq(380, 700, by = 10),
	y = prep_cmf(lambda_nm = seq(380, 700, by = 60))$cmf_xyz[, 3],
	col = "blue",
	type = "l"
)
lines(
	x = seq(380, 700, by = 10),
	y = prep_cmf(lambda_nm = seq(380, 700, by = 60))$cmf_xyz[, 2],
	col = "green"
)
lines(
	x = seq(380, 700, by = 10),
	y = prep_cmf(lambda_nm = seq(380, 700, by = 60))$cmf_xyz[, 1],
	col = "red"
)


lines(
	x = seq(380, 700, by = 60),
	prep_cmf(lambda_nm = seq(380, 700, by = 60))$cmf_xyz[, 3],
	col = "blue",
	lty = "dashed"
)
lines(
	x = seq(380, 700, by = 60),
	prep_cmf(lambda_nm = seq(380, 700, by = 60))$cmf_xyz[, 2],
	col = "green",
	lty = "dashed"
)
lines(
	x = seq(380, 700, by = 60),
	prep_cmf(lambda_nm = seq(380, 700, by = 60))$cmf_xyz[, 1],
	col = "red",
	lty = "dashed"
)
