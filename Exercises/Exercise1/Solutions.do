/****************************************/
/*    Example Solution for Exercise 1   */
/* Written by Peter Hull, 3/6/2022 (v1) */
/****************************************/

/* Install any needed packages */
cap ssc install binscatter
cap ssc install ranktest
cap ssc install ivreg2

/* Setup data */
clear all
use angrist_krueger_91, clear

forval q=1/4 {
	gen qob_`q'=(qob==`q')
}

/* OLS and binscatter */
reg lwage educ, r
binscatter lwage educ, xtitle("Years of Completed Schooling") ytitle("Log Wages")

/* Simple (Wald) IV */
ivreg2 lwage (educ=qob_1), r
foreach var of varlist lwage educ {
	summ `var' if qob_1==1
	summ `var' if qob_1==0
}

/*Overidentified IV */
ivreg2 lwage (educ=qob_1 qob_2 qob_3), r
sto
preserve
collapse (mean) lwage educ, by(qob)
scatter lwage educ || lfit lwage educ
reg lwage educ, r
restore

/*Putting the 2S in 2SLS */
ivregress 2sls lwage (educ=qob_1 qob_2 qob_3) i.yob, r
reg educ qob_1 qob_2 qob_3 i.yob
predict educ_hat, xb
reg lwage educ_hat i.yob, r

/* Many-IV Bias */
ivreg2 lwage (educ=qob_1#yob qob_2#yob qob_3#yob) i.yob, r
