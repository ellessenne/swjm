version 18

// Log
capture log close
log using "02b-data-analysis-log.txt", text replace

// Constant intervention model

use "02-dt-ci.dta", clear
compress
save, replace

// Linear mixed-effects model
gsem (yobs <- ibn.j i.x M1[id]@1, noconstant family(gaussian))
estimates store mod_ci_lmm

// The equivalent model with -mixed- is:
mixed yobs ibn.j i.x, noconstant || id:

// Joint model
gsem ///
	(yobs <- ibn.j i.x M1[id]@1, noconstant family(gaussian)) ///
	(t <- i.x M1[id], family(weibull, failure(d) lt(t0)))
estimates store mod_ci_jm

// Compare the estimates
estimates table mod_ci_lmm mod_ci_jm

// General time on treatment model

use "02-dt-gi.dta", clear
compress
save, replace

// Linear mixed-effects model
gsem (yobs <- ibn.j i.cumx M1[id]@1, noconstant family(gaussian))
estimates store mod_gi_lmm

// The equivalent model with -mixed- is:
mixed yobs ibn.j i.cumx, noconstant || id:

// Joint model
gsem ///
	(yobs <- ibn.j i.cumx M1[id]@1, noconstant family(gaussian)) ///
	(t <- i.x M1[id], family(weibull, failure(d) lt(t0)))
estimates store mod_gi_jm

// Compare the estimates
estimates table mod_gi_lmm mod_gi_jm

// Close log
log close
