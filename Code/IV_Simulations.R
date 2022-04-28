# ------------------------------------------------------------------------------
# IV Simulations
# Written by Peter Hull, 11/15/21 (v1)
# ------------------------------------------------------------------------------

library(fixest)
library(purrr)

# ------------------------------------------------------------------------------
# Simulation showing weak instruments move IV estimate -> OLS estimate
#
# DGP:
# x = \Pi * z1 + eta, z1-z100 ~ N(0,1) and eta ~ N(0, 0.25)
# y = eps, eps = eta + N(0, 0.25)
#
# True parameter value of X is 0.
# Note that X is correlated with epsilon via eta, so OLS is biased.
# ------------------------------------------------------------------------------

results_weak <- purrr::map_dfr(1:2000, \(i) {
	n_units <- 100
	df <- data.table(z = rnorm(n_units))

	df[,
		eta := 0.5 * rnorm(n_units)
	][,
		eps := eta + 0.5 * rnorm(n_units),
	][,
		y := eps
	]

	# \Pi = 1
	df[,
		 x := 1 * z + eta
	]
	est_strong <- feols(y ~ 1 | 0 | x ~ z, df)

	# \Pi = 0.1
	df[,
		 x := 0.1 * z + eta
	]
	est_weaker <- feols(y ~ 1 | 0 | x ~ z, df)

	# \Pi = 0.01
	df[,
		 x := 0.01 * z + eta
	]
	est_weakest <- feols(y ~ 1 | 0 | x ~ z, df)


	return(data.table(
		name = c("Strong", "Weaker", "Weakest"),
		est  = c(
			coef(est_strong)[["fit_x"]],
			coef(est_weaker)[["fit_x"]],
			coef(est_weakest)[["fit_x"]]
		)
	))
})


ggplot(results_weak) +
	# Plot estimates
	geom_density(
		aes(x = est, group = name, color = name)
	) +
	xlim(-2, 2) +
	# True Parameter
	geom_vline(xintercept = 0) +
	annotate(
		"text", label = "True Parameter",
		x = 0.05, y = 3.75, size = 6, hjust = 0
	) +
	# OLS Estimate
	geom_vline(xintercept = 1) +
	annotate(
		"text", label = "OLS Coefficient",
		x = 1.05, y = 3.75, size = 6, hjust = 0
	)


# ------------------------------------------------------------------------------
# Simulation showing many instruments move IV estimate -> OLS estimate
#
# DGP:
# x = \Pi * z1 + eta, z1-z100 ~ N(0,1) and eta ~ N(0, 0.25)
# y = eps, eps = 5 * eta + N(0, 0.25)
#
# True parameter value of X is 0.
# Note that X is correlated with epsilon via eta, so OLS is biased.
# ------------------------------------------------------------------------------


results_many <- purrr::map_dfr(1:2000, \(i) {
	n_units <- 100

	# Generate many instruments
	df <- data.table(z1 = rnorm(n_units))
	for(i in 2:100) {
		df[[paste0("z", i)]] <- rnorm(n_units)
	}

	df[,
		eta := 0.5 * rnorm(n_units)
	][,
		eps := 5 * eta + 0.5 * rnorm(n_units),
	][,
		y := eps
	][,
		# \Pi = 1
		x := 1 * z1 + eta
	]


	# IV with 1, 10, and 100 instruments
	est_1 <- feols(y ~ 1 | 0 | x ~ z.[1], df)

	est_10 <- feols(y ~ 1 | 0 | x ~ z.[1:10], df)

	est_100 <- feols(y ~ 1 | 0 | x ~ z.[1:100], df)

	return(data.table(
		name = c("1 Instrument", "10 Instruments", "100 Instruments"),
		est  = c(
			coef(est_1)[["fit_x"]],
			coef(est_10)[["fit_x"]],
			coef(est_100)[["fit_x"]]
		)
	))
})



ggplot(results_many) +
	# Plot estimates
	geom_density(
		aes(x = est, group = name, color = name)
	) +
	xlim(-2, 2) +
	# True Parameter
	geom_vline(xintercept = 0) +
	annotate(
		"text", label = "True Parameter",
		x = 0.05, y = 3.75, size = 6, hjust = 0
	) +
	# OLS Estimate
	geom_vline(xintercept = 1) +
	annotate(
		"text", label = "OLS Coefficient",
		x = 1.05, y = 3.75, size = 6, hjust = 0
	)
