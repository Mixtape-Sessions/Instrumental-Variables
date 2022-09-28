## Example Solutions for Exercise 2 --------------------------------------------
## Kyle Butts, CU Boulder Economics
##
## Solutions for Peter Hull's Instrumental Variables Mixtape Session


library(data.table)
library(haven)
library(fixest)
library(SteinIV)

# ---- Setup Data --------------------------------------------------------------
df <- read_dta("https://github.com/Mixtape-Sessions/Instrumental-Variables/blob/main/Data/stevenson.dta?raw=true")
setDT(df)
df <- df[black == 1, ]

# ---- Get dummy IV ------------------------------------------------------------
est <- feols(jail3 ~ judge_pre_.[1:7], df)
df[, leniency := predict(est)]
df[, more_lenient := (leniency > median(leniency))]

# ---- Naive OLS and 2SLS ------------------------------------------------------
feols(guilt ~ jail3, df, vcov="hc1")
feols(guilt ~ 1 | 0 | jail3 ~ i(more_lenient), df, vcov="hc1")


# ---- Balance on prior felonies -----------------------------------------------
feols(prior_felChar ~ i(more_lenient), df, vcov="hc1")
feols(prior_felChar ~ i(more_lenient) | bailDate, df, vcov="hc1")

# ---- Controlled OLS and 2SLS -------------------------------------------------
feols(guilt ~ jail3 | bailDate, df, vcov="hc1")
feols(guilt ~ 1 | bailDate | jail3 ~ more_lenient, df, vcov="hc1")

# ---- Estimate complier Y0 and compare to overall E[Y|D=0] --------------------
df[,
	 Y_omD := guilt * (1 - jail3)
][,
	omD := (1 - jail3)
]

feols(Y_omD ~ 1 | bailDate | omD ~ i(more_lenient), df, vcov="hc1")

summary(df[jail3 == 0, "guilt"])

# ---- 2SLS using all judges ---------------------------------------------------
feols(guilt ~ 1 | 0 | jail3 ~ judge_pre_.[1:7], df, vcov="hc1")
feols(guilt ~ 1 | 0 | jail3 ~ leniency, df, vcov="hc1")


# ---- JIVE using all judges ---------------------------------------------------

jive.est(
	y = df$guilt, X = as.matrix(df[,.(intercept = rep(1, .N), jail3)]),
	Z = as.matrix(df[,.SD, .SDcols=patterns("judge_pre_")])
)

# Make judge ID
df[, judge := (judge_pre_1 == 1) * 1 +
	 	(judge_pre_2 == 1) * 2 +
	 	(judge_pre_3 == 1) * 3 +
	 	(judge_pre_4 == 1) * 4 +
	 	(judge_pre_5 == 1) * 5 +
	 	(judge_pre_6 == 1) * 6 +
	 	(judge_pre_7 == 1) * 7 +
	 	(judge_pre_8 == 1) * 8]

df[, num := .N, by = judge]
df[, lo_leniency := (leniency*num - jail3)/(num - 1)]

feols(guilt ~ 1 | 0 | jail3 ~ lo_leniency, df, vcov="hc1")


# ---- 2SLS / JIVE / UJIVE using all judges  & controls ------------------------

feols(guilt ~ 1 | bailDate | jail3 ~ judge_pre_.[1:7], df, vcov="hc1")

# Jonathan Seward ported UJIVE to R
# Can't use bailDate fixed effects since the matrix is too large to invert
source("https://raw.githubusercontent.com/j-seward/ujive/main/ivreg_001.R")

ivreg(
	y = df$guilt, T = df$jail3,
	Z = as.matrix(df[,.SD, .SDcols=patterns("judge_pre_")]),
	W = as.matrix(df[, .(rep(1, .N))]), # Constant only
	noConstant = FALSE
)
