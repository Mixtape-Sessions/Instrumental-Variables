{smcl}
{* *! version 0.5.0 14Nov2021}{...}
{viewerdialog manyiv "dialog manyiv"}{...}
{vieweralsosee "[R] manyiv" "mansection R manyiv"}{...}
{viewerjumpto "Syntax" "manyiv##syntax"}{...}
{viewerjumpto "Description" "manyiv##description"}{...}
{viewerjumpto "Options" "manyiv##options"}{...}
{viewerjumpto "Examples" "manyiv##examples"}{...}
{title:Title}

{p2colset 5 17 17 2}{...}
{p2col :{cmd:manyiv} {hline 2}}Many IV regressions (OLS, TSLS, LIML, MBTSLS, JIVE, UJIVE, RTSLS){p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{pstd}
Run multiple IV regressions:

{p 8 15 2}
{cmd:manyiv}
{depvar}
{cmd:(}{it:endogenous}{cmd:=}{it:instruments}{cmd:)}
[{it:exogenous}]
{ifin}
[{cmd:,} {it:{help manyiv##table_options:options}}]

{synoptset 18 tabbed}{...}
{marker table_options}{...}
{synopthdr}
{synoptline}
{synopt :{opth absorb(varlist)}} Controls to absorb as fixed effects.
{p_end}
{synopt :{opth absorbiv(varlist)}} Instruments to absorb as fixed effects.
{p_end}
{synopt :{opth cluster(varname)}} ClusterSEs by variable.
{p_end}
{synopt :{opt skipsingletons}} Skip singleton absorb groups.
{p_end}
{synopt :{opt keepsingletons}} Keep singleton absorb groups (jive/ujive not estimated with this option).
{p_end}
{synopt :{opth save:results(str)}} Save results into mata object.
{p_end}
{synopt :{opt noc:onstant}} Do not include constant.
{p_end}
{synopt :{opt nosmall}} Small-sample adjustmnt.
{p_end}
{synopt :{opt noprint}} Do not print table.
{p_end}
{synopt :{opt nose}} Do not compute se.
{p_end}
{synopt :{opt nostats}} Do not compute stats.
{p_end}
{synopt :{opt nosquarem}} No SQUAREM acceeration (multiple absorb only).
{p_end}

{p2colreset}{...}
{p 4 6 2}

{marker description}{...}
{title:Description}

{pstd}
This package computes several instrumental variables estimators as well as
many-instrument robust SEs for selected estimators. In particular, the
package computes OLS, TSLS, LIML, MBTSLS, JIVE, UJIVE, RTSLS.

{pstd}
LIML is Limited information maximum likelihood, MBTSLS is the Kolesar,
Chetty, Friedman, Glaeser and Imbens (2011) modified bias-corrected
two-stage least squares, JIVE is the JIVE1 estimator in Angrist, Imbens,
and Krueger (1999), and UJIVE is the Kolesar (2012) version of the JIVE
estimator. RTSLS is the reverse two-stage least squares estimator.

{pstd}
The Standard errors are returned in a 5x7 matrix. The first row gives
standard errors valid under homoscedasticity and standard asymptotics
(classic standard errors). The second row computes robust standard
errors that are valid under heteroscedasticity.  If {opt cluster()}
is specified, the third row computes clustered standard errors.  The fourth
row computes many-instrument robust standard errors that are
valid under homoscedasticity and Bekker asymptotics, allowing also for
the presence of many covariates. The fifth row computes estimates of
the standard errors that are valid under the many invalid instruments
sequence of Kolesar, Chetty, Friedman, Glaeser and Imbens (2011)

{pstd}
{cmd:manyiv} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(F)}}the first-stage F-statistic{p_end}
{synopt:{cmd:e(small)}}small-sample adjustment{p_end}
{synopt:{cmd:e(jive)}}jive/ujive estimated{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)            }}Coefficient vector{p_end}
{synopt:{cmd:e(se)           }}Matrix with SEs{p_end}
{synopt:{cmd:e(rf)           }}Vector with reduced form (if non-absorb instruments){p_end}
{synopt:{cmd:e(fs)           }}Vector with first stage (if non-absorb instruments){p_end}
{synopt:{cmd:e(Omega)        }}Estimate of the reduced-form covariance matrix{p_end}
{synopt:{cmd:e(Xi)           }}Estimate of XI{p_end}
{synopt:{cmd:e(Sargan)       }}2-by-1 vector, with the first element equal to the Sargan test statistic and the second element equal to the p-value. (Missing for single instrument).{p_end}
{synopt:{cmd:e(CD)           }}2-by-1 vector, with the first element equal to the Cragg-Donald test statistic and the second element equal to the p-value. The p-value contains a size-correction derived in Kolesar (2012) that ensures correct coverage under many-instrument asymptotics. (Missing for single instrument).{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}

{marker example}{...}
{title:Examples}

{phang2}{cmd:. clear                                            }{p_end}
{phang2}{cmd:. set seed 1729                                    }{p_end}
{phang2}{cmd:. set obs 1000                                     }{p_end}
{phang2}{cmd:. gen u  = rnormal()                               }{p_end}
{phang2}{cmd:. gen z1 = rnormal()                               }{p_end}
{phang2}{cmd:. gen z2 = rnormal()                               }{p_end}
{phang2}{cmd:. gen e  = rnormal() + u                           }{p_end}
{phang2}{cmd:. gen c  = int(runiform() * 10)                    }{p_end}
{phang2}{cmd:. gen fe = int(runiform() * 15)                    }{p_end}
{phang2}{cmd:. gen iv = int(runiform() * 8)                     }{p_end}
{phang2}{cmd:. gen w  = rnormal()                               }{p_end}
{phang2}{cmd:. gen x  = 1 + 0.1 * z1 - 0.2 * z2 - 1/(1 + iv) + u}{p_end}
{phang2}{cmd:. gen y  = 1 + x + w + 1/(1 + fe) + e              }{p_end}

{phang2}{cmd:. manyiv y (x = z1 z2) w                                    }{p_end}
{phang2}{cmd:. manyiv y (x = z1 z2) w, cluster(c)                        }{p_end}
{phang2}{cmd:. manyiv y (x = z1 z2) w, absorb(fe) cluster(c)             }{p_end}
{phang2}{cmd:. manyiv y (x = z1 z2) w, absorbiv(iv) cluster(c)           }{p_end}
{phang2}{cmd:. manyiv y (x = z1 z2) w, absorb(fe) absorbiv(iv) cluster(c)}{p_end}
{phang2}{cmd:. manyiv y (x = .)     w, absorb(fe) absorbiv(iv) cluster(c)}{p_end}
