/****************************************/
/*    Example Solution for Exercise 2   */
/* Written by Peter Hull, 3/9/2022 (v1) */
/****************************************/

/* Install any needed packages */
* cap ssc install ranktest
* cap ssc install ivreg2

* To install jive: 
* -findit jive-

* To install manyiv from local file:
* type -sysdir-, and copy all the manyiv files into PERSONAL folder

/* Setup data */
clear all
use "https://github.com/Mixtape-Sessions/Instrumental-Variables/blob/main/Data/stevenson.dta?raw=true", clear
keep if black==1

/* Get dummy IV */
reg jail3 judge_pre_1-judge_pre_7
predict leniency, xb
summ leniency, d
gen more_lenient=leniency>r(p50)

/* Naive OLS and 2SLS */
reg guilt jail3, r
ivreg2 guilt (jail3=more_lenient), r

/* Balance on prior felonies */
reg prior_felChar more_lenient, r
reghdfe prior_felChar more_lenient, vce(robust) absorb(bailDate)

/* Controlled OLS and 2SLS */
reghdfe guilt jail3, absorb(bailDate) vce(robust)
ivreghdfe guilt (jail3=more_lenient), absorb(bailDate) r
 
/* Estimate complier Y0 and compare to overall E[Y|D=0] */
gen Y_omD=guilt*(1-jail3)
gen omD=1-jail3
ivreghdfe Y_omD (omD=more_lenient), absorb(bailDate) r
summ guilt if jail3==0

/* 2SLS using all judges */
ivreghdfe guilt (jail3=judge_pre_1-judge_pre_7), r 
ivreghdfe guilt (jail3=leniency), r 

/* JIVE using all judges */
jive guilt (jail3=judge_pre_1-judge_pre_7), r 
egen judge=group(judge_pre_*)
bys judge: egen num=count(jail3)
gen lo_leniency=(leniency*num-jail3)/(num-1)
ivreghdfe guilt (jail3=lo_leniency), r 

/* 2SLS / JIVE / UJIVE using all judges  & controls */
ivreghdfe guilt (jail3=judge_pre_1-judge_pre_7), absorb(bailDate) r 
manyiv guilt (jail3 = judge_pre_1-judge_pre_7), absorb(bailDate) 
