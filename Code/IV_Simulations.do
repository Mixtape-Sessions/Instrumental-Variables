/*------------------------------------------------------------------------------
IV Simulations
Written by Peter Hull, 11/15/21 (v1)
------------------------------------------------------------------------------*/ 


clear all
set matsize 2000
set seed 42


/*------------------------------------------------------------------------------
Simulation showing weak instruments move IV estimate -> OLS estimate

DGP:
x = \Pi * z1 + eta, z1-z100 ~ N(0,1) and eta ~ N(0, 0.25)
y = eps, eps = eta + N(0, 0.25)

True parameter value of X is 0. 
Note that X is correlated with epsilon via eta, so OLS is biased.
------------------------------------------------------------------------------*/ 

matrix sims=J(2000,3,.)
forval s=1/200 { /* set to 2000 to get the graphs from class */
	disp `s'
	qui {
	clear
	set obs 100
	gen z = rnormal()
	gen eta = 0.5*rnormal()
	gen eps = eta+0.5*rnormal()
	gen y = eps /* Note: no treatment effect */
	gen x = 1*z + eta
	ivregress 2sls y (x=z)
		matrix sims[`s',1]=_b[x]
	replace x = 0.1*z + eta
	ivregress 2sls y (x=z)
		matrix sims[`s',2]=_b[x]
	replace x = 0.01*z + eta
	ivregress 2sls y (x=z)
		matrix sims[`s',3]=_b[x]
	}
}
clear
svmat sims
twoway (kdensity sims1 if sims1>-2 & sims1<2, xline(0, lpattern(dash) lcolor(black)) xline(1, lcolor(black) lpattern(shortdash)) ///
	 graphregion(color(white)) ytitle("Density") xtitle("2SLS Estimate") xlab(-2(1)2) ///
	 text(6.2 -0.4 "True Parameter") text(6.2 1.6 "Regression Coefficient"))
graph export strongpi.png, replace

twoway (kdensity sims2 if sims2>-2 & sims2<2, xline(0, lpattern(dash) lcolor(black)) xline(1, lcolor(black) lpattern(shortdash)) ///
	 graphregion(color(white)) ytitle("Density") xtitle("2SLS Estimate") xlab(-2(1)2) ///
	 text(0.81 -0.4 "True Parameter") text(0.81 1.6 "Regression Coefficient"))
graph export medpi.png, replace

twoway (kdensity sims3 if sims3>-2 & sims3<2, xline(0, lpattern(dash) lcolor(black)) xline(1, lcolor(black) lpattern(shortdash)) ///
	 graphregion(color(white)) ytitle("Density") xtitle("2SLS Estimate") xlab(-2(1)2) ///
	 text(0.51 -0.4 "True Parameter") text(0.51 1.6 "Regression Coefficient"))
graph export weakpi.png, replace



/*------------------------------------------------------------------------------
Simulation showing many instruments move IV estimate -> OLS estimate

DGP:
x = z1 + eta, z1-z100 ~ N(0,1) and eta ~ N(0, 0.25)
y = eps, eps = 5 * eta + N(0, 0.25)

True parameter value of X is 0. 
Note that X is correlated with epsilon via eta, so OLS is biased.
------------------------------------------------------------------------------*/ 

matrix sims=J(2000,3,.)
forval s=1/200 { /* set to 2000 to get the graphs from class */
	disp `s'
	qui {
	clear
	set obs 100
	forval j=1/100 {
		gen z_`j' = rnormal()
	}
	gen eta = 0.5*rnormal()
	gen eps = 5*eta+0.5*rnormal()
	gen y = eps /* Note: no treatment effect */
	gen x = 1*z_1 + eta
	ivregress 2sls y (x=z_1)
		matrix sims[`s',1]=_b[x]
	ivregress 2sls y (x=z_1-z_10)
		matrix sims[`s',2]=_b[x]
	ivregress 2sls y (x=z_*)
		matrix sims[`s',3]=_b[x]
	}
}
clear
svmat sims
twoway (kdensity sims1 if sims1>-2 & sims1<2, xline(0, lpattern(dash) lcolor(black)) xline(1, lcolor(black) lpattern(shortdash)) ///
	 graphregion(color(white)) ytitle("Density") xtitle("2SLS Estimate") xlab(-2(1)2) ///
	 text(1.55 -0.4 "True Parameter") text(1.55 1.6 "Regression Coefficient"))
graph export fewz.png, replace

twoway (kdensity sims2 if sims2>-2 & sims2<2, xline(0, lpattern(dash) lcolor(black)) xline(1, lcolor(black) lpattern(shortdash)) ///
	 graphregion(color(white)) ytitle("Density") xtitle("2SLS Estimate") xlab(-2(1)2) ///
	 text(1.55 -0.4 "True Parameter") text(1.55 1.6 "Regression Coefficient"))
graph export somez.png, replace

twoway (kdensity sims3 if sims3>-2 & sims3<2, xline(0, lpattern(dash) lcolor(black)) xline(1, lcolor(black) lpattern(shortdash)) ///
	 graphregion(color(white)) ytitle("Density") xtitle("2SLS Estimate") xlab(-2(1)2) ///
	 text(2 -0.4 "True Parameter") text(2 1.6 "Regression Coefficient"))
graph export manyz.png, replace
