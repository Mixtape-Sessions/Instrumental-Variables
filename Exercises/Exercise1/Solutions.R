## Example Solutions for Exercise 1 --------------------------------------------
## Kyle Butts, CU Boulder Economics
##
## Solutions for Peter Hull's Instrumental Variables Mixtape Session

library(haven) # Read .dta files
library(data.table) # For working with data
library(fixest) # For regressions
library(binsreg) # For binscatter
library(ggplot2)

## Load data
# data <- haven::read_dta("https://github.com/Mixtape-Sessions/Instrumental-Variables/blob/main/Exercises/Exercise1/angrist_krueger_91.dta?raw=true")
data <- read_dta("~/Downloads/angrist_krueger_91.dta")
data <- as.data.table(data)

data[, qob_1 := (qob == 1)]
data[, qob_2 := (qob == 2)]
data[, qob_3 := (qob == 3)]
data[, qob_4 := (qob == 4)]



# ---- OLS and Binscatter ------------------------------------------------------

feols(
	lwage ~ educ, # Regression formula
	data,
	vcov = "hc1" # ,r
)

binscatter <- binsreg(data$lwage, data$educ)

# plot and add labels
binscatter$bins_plot +
	labs(y = "Log wages", x = "Years of Completed Schooling")



# ---- Simple (Wald) IV Estimator ----------------------------------------------

# Formula y ~ exogenous | fixed effects | endogenous ~ instrument
# 1 = constant, 0 = no fixed effects
feols(
	lwage ~ 1 | 0 | educ ~ qob_1,
	data,
	vcov = "hc1"
)


data[,
		 .(n = .N, mean = mean(lwage), sd = sd(lwage), min = min(lwage), max = max(lwage)),
		 by = qob_1
]
data[,
		 .(n = .N, mean = mean(educ), sd = sd(educ), min = min(educ), max = max(educ)),
		 by = qob_1
]



# ---- Overidentified IV Estimator ---------------------------------------------

feols(
	lwage ~ 1 | 0 | educ ~ qob_1 + qob_2 + qob_3,
	data,
	vcov = "hc1"
)

# collapse data by qob
collapsed <- data[,
		 .(lwage = mean(lwage), educ = mean(educ)),
		 by = qob
]

# plot means
plot(collapsed$educ, collapsed$lwage)

# add regression line
abline(feols(lwage ~ educ, collapsed))



# ---- Putting the 2S in 2SLS --------------------------------------------------

feols(
	lwage ~ 1 | yob | educ ~ qob_1 + qob_2 + qob_3,
	data,
	vcov = "hc1"
)

first_stage <- feols(educ ~ i(qob) | yob, data)

data[, educ_hat := predict(first_stage)]

feols(
	lwage ~ educ_hat | yob,
	data,
	vcov = "hc1"
)



# ---- Many IV Bias ------------------------------------------------------------

feols(
	lwage ~ 1 | yob | educ ~ i(yob, qob_1) + i(yob, qob_2) + i(yob, qob_3),
	data,
	vcov = "hc1"
)

