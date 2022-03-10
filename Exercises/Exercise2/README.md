# Mixtape IV Workshop: Coding Lab 2

This lab will walk through some standard judge IV analyses using data from Stevenson (2018): a paper you can find in our Readings folder. Stevenson leverages the quasi-random assignment of 8 judges (magistrates) in Philadelphia to study the effects pretrial detention on several outcomes, including whether or not a defendant subsequently pleads guilty. To start, load the `stevenson.dta` data file into Stata or R. This is the same dataset from Scott’s Mixtape textbook.

1. Limit the sample to Black defendants (`black==1`) and compute the “leniency” of a defendant’s assigned judge. This is the average pretrial detention rate (indicated by `jail3`) by the eight judges (indicated by `judge_pre_1-judge_pre_7`). Generate an indicator for a defendant having an above-median lenient judge; call this variable `more_lenient`. Estimate OLS and 2SLS regressions of the guilty plea outcome (`guilt`) on pretrial detention, with the latter instrumenting by `more_lenient`. Report your coefficients and robust standard errors. What do you find interesting here?

  - OLS Estimate:
  - Standard Error:
  
  - 2SLS Estimate:
  - Standard Error:

  - Interesting Things:

2. Show that this instrument is correlated with a defendant having a prior felony charge (indicated by `prior_felChar==1`). What assumption of the LATE theorem appears to be violated, given such a correlation? Show that this correlation mostly goes away when controlling for date fixed effects (coded in `bailDate`). Re-run your OLS and 2SLS estimates with these controls and comment on the change.

  - Controlled OLS Estimate:
  - Standard Error:

  - Controlled 2SLS Estimate:
  - Standard Error:

  - Interesting Things:

3. Estimate the average untreated potential outcome for compliers in the controlled 2SLS specification, using the simple trick we saw in lecture. Compare this to the average outcome among `jail3==0` defendants in the full population. Interpret the difference

  - E[Y(0)|D(1)>D(0)]:
  - E[Y(0)|D=0]:

4. Returning to the 2SLS specification without controls, replace the single `more_lenient` instrument with indicators for seven of the eight judges. Show that you get an identical 2SLS estimate and standard error if you instrument by the assigned judge leniency.

  - Overidentified 2SLS Estimate:
  - Standard Error:

5. Estimate the coefficient of interest by JIVE, _without_ controls, using the judge indicators as seven instruments. Show that you get a (basically) identical 2SLS estimate and standard error if you instrument by the _leave-out_ assigned judge leniency. That is, the average pretrial detention rate among other defendants assigned to each defendant’s judge. 

  - JIVE Estimate:
  - Standard Error:

6. Add date fixed effects to the above 2SLS and JIVE specifications, and comment on how the coefficients change. Finally, estimate the coefficient by Kolesar’s UJIVE with controls and judge indicators as seven instruments. You should use the `manyiv` Stata command for this; unzip the `manyiv-0.5.0.zip` file included in this repository and put its contents in the folder you’re running Stata from (and check out the help file). Comment on the different estimates. *Note:* R doesn't really have a `manyiv` function

  - Overidentified+Controlled 2SLS Estimate:
  - Standard Error:

  - Overidentified+Controlled JIVE Estimate:
  - Standard Error:

  - Overidentified+Controlled UJIVE Estimate:
  - Standard Error:

