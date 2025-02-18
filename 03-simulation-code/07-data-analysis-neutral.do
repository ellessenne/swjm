version 18

cd "~/Documents/stepped-wedge-jm" // ... or wherever

// Local values
local scenarios = 1
local B = 1000

// Constant intervention model

// Create postfile
capture postclose pf
postfile pf double(i model scenario parameter b stderr converged error) str20(name) using "simulation/data/07-full-results-ci.dta", replace

// Loop over scenarios and repetitions
forvalues s = 1/`scenarios' {
	nois _dots 0, title(Scenario: `s') reps(`B')
	forvalues i = 1/`B' {
		quietly {
			use "simulation/data/simdata/06-full-ci-`s'-`i'", clear
			* use "simulation/data/simdata/06-full-ci-1-2", clear
			compress
			save, replace
			// LMM
			capture gsem ///
				(yobs <- ibn.j i.x M1[i]@1 M2[i>id]@1, noconstant family(gaussian)) ///
				, covstructure(M1[i] M2[i>id], diag) startgrid technique(bfgs)
			// Post results
			if (_rc > 0 | e(converged) == 0) {
				// gsem failed or did not converge
				post pf (`i') (2) (`s') (1) (.) (.) (.) (99) ("beta1")
				post pf (`i') (2) (`s') (2) (.) (.) (.) (99) ("beta2")
				post pf (`i') (2) (`s') (3) (.) (.) (.) (99) ("beta3")
				post pf (`i') (2) (`s') (4) (.) (.) (.) (99) ("beta4")
				post pf (`i') (2) (`s') (5) (.) (.) (.) (99) ("beta5")
				post pf (`i') (2) (`s') (6) (.) (.) (.) (99) ("delta")
				post pf (`i') (2) (`s') (7) (.) (.) (.) (99) ("nu")
				post pf (`i') (2) (`s') (8) (.) (.) (.) (99) ("omega1")
				post pf (`i') (2) (`s') (9) (.) (.) (.) (99) ("omega3")
				post pf (`i') (2) (`s') (10) (.) (.) (.) (99) ("ln_lambda")
				post pf (`i') (2) (`s') (11) (.) (.) (.) (99) ("ln_p")
				post pf (`i') (2) (`s') (12) (.) (.) (.) (99) ("sigma2_alpha")
				post pf (`i') (2) (`s') (13) (.) (.) (.) (99) ("sigma2_phi")
				post pf (`i') (2) (`s') (14) (.) (.) (.) (99) ("sigma2_epsilon")
				// Then, parameter = 20 is ICCa, parameter = 21 is ICCw
				post pf (`i') (2) (`s') (20) (.) (.) (.) (99) ("icca")
				post pf (`i') (2) (`s') (21) (.) (.) (.) (99) ("iccw")
			}
			else {
				// gsem ok
				post pf (`i') (2) (`s') (1) (_b[yobs:1.j]) (_se[yobs:1.j]) (e(converged)) (0) ("beta1")
				post pf (`i') (2) (`s') (2) (_b[yobs:2.j]) (_se[yobs:2.j]) (e(converged)) (0) ("beta2")
				post pf (`i') (2) (`s') (3) (_b[yobs:3.j]) (_se[yobs:3.j]) (e(converged)) (0) ("beta3")
				post pf (`i') (2) (`s') (4) (_b[yobs:4.j]) (_se[yobs:4.j]) (e(converged)) (0) ("beta4")
				post pf (`i') (2) (`s') (5) (_b[yobs:5.j]) (_se[yobs:5.j]) (e(converged)) (0) ("beta5")
				post pf (`i') (2) (`s') (6) (_b[yobs:1.x]) (_se[yobs:1.x]) (e(converged)) (0) ("delta")
				post pf (`i') (2) (`s') (7) (.) (.) (e(converged)) (0) ("nu")
				post pf (`i') (2) (`s') (8) (.) (.) (e(converged)) (0) ("omega1")
				post pf (`i') (2) (`s') (9) (.) (.) (e(converged)) (0) ("omega3")
				post pf (`i') (2) (`s') (10) (.) (.) (e(converged)) (0) ("ln_lambda")
				post pf (`i') (2) (`s') (11) (.) (.) (e(converged)) (0) ("ln_p")
				post pf (`i') (2) (`s') (12) (_b[/:var(M1[i])]) (_se[/:var(M1[i])]) (e(converged)) (0) ("sigma2_alpha")
				post pf (`i') (2) (`s') (13) (_b[/:var(M2[i>id])]) (_se[/:var(M2[i>id])]) (e(converged)) (0) ("sigma2_phi")
				post pf (`i') (2) (`s') (14) (_b[/:var(e.yobs)]) (_se[/:var(e.yobs)]) (e(converged)) (0) ("sigma2_epsilon")
				// Calculate ICC and post that too, with SE using numerical delta method
				nlcom (_b[/:var(M1[i])] + _b[/:var(M2[i>id])]) / (_b[/:var(M1[i])] + _b[/:var(M2[i>id])] + _b[/:var(e.yobs)]), iterate(100000)
				post pf (`i') (2) (`s') (20) (r(b)[1,1]) (sqrt(r(V)[1,1])) (e(converged)) (0) ("icca")
				nlcom _b[/:var(M1[i])] / (_b[/:var(M1[i])] + _b[/:var(M2[i>id])] + _b[/:var(e.yobs)]), iterate(100000)
				post pf (`i') (2) (`s') (21) (r(b)[1,1]) (sqrt(r(V)[1,1])) (e(converged)) (0) ("iccw")
			}
			// JM
			capture gsem ///
				(yobs <- ibn.j i.x M1[i]@1 M2[i>id]@1, noconstant family(gaussian)) ///
				(t <- i.x M1[i] M2[i>id], family(weibull, failure(d) lt(t0))) ///
				, covstructure(M1[i] M2[i>id], diag) startgrid technique(bfgs)
			// Post results
			if (_rc > 0 | e(converged) == 0) {
				// gsem failed or did not converge
				post pf (`i') (1) (`s') (1) (.) (.) (.) (99) ("beta1")
				post pf (`i') (1) (`s') (2) (.) (.) (.) (99) ("beta2")
				post pf (`i') (1) (`s') (3) (.) (.) (.) (99) ("beta3")
				post pf (`i') (1) (`s') (4) (.) (.) (.) (99) ("beta4")
				post pf (`i') (1) (`s') (5) (.) (.) (.) (99) ("beta5")
				post pf (`i') (1) (`s') (6) (.) (.) (.) (99) ("delta")
				post pf (`i') (1) (`s') (7) (.) (.) (.) (99) ("nu")
				post pf (`i') (1) (`s') (8) (.) (.) (.) (99) ("omega1")
				post pf (`i') (1) (`s') (9) (.) (.) (.) (99) ("omega3")
				post pf (`i') (1) (`s') (10) (.) (.) (.) (99) ("ln_lambda")
				post pf (`i') (1) (`s') (11) (.) (.) (.) (99) ("ln_p")
				post pf (`i') (1) (`s') (12) (.) (.) (.) (99) ("sigma2_alpha")
				post pf (`i') (1) (`s') (13) (.) (.) (.) (99) ("sigma2_phi")
				post pf (`i') (1) (`s') (14) (.) (.) (.) (99) ("sigma2_epsilon")
				// Then, parameter = 20 is ICCa, parameter = 21 is ICCw
				post pf (`i') (1) (`s') (20) (.) (.) (.) (99) ("icca")
				post pf (`i') (1) (`s') (21) (.) (.) (.) (99) ("iccw")
			}
			else {
				// gsem ok
				post pf (`i') (1) (`s') (1) (_b[yobs:1.j]) (_se[yobs:1.j]) (e(converged)) (0) ("beta1")
				post pf (`i') (1) (`s') (2) (_b[yobs:2.j]) (_se[yobs:2.j]) (e(converged)) (0) ("beta2")
				post pf (`i') (1) (`s') (3) (_b[yobs:3.j]) (_se[yobs:3.j]) (e(converged)) (0) ("beta3")
				post pf (`i') (1) (`s') (4) (_b[yobs:4.j]) (_se[yobs:4.j]) (e(converged)) (0) ("beta4")
				post pf (`i') (1) (`s') (5) (_b[yobs:5.j]) (_se[yobs:5.j]) (e(converged)) (0) ("beta5")
				post pf (`i') (1) (`s') (6) (_b[yobs:1.x]) (_se[yobs:1.x]) (e(converged)) (0) ("delta")
				post pf (`i') (1) (`s') (7) (_b[t:1.x]) (_se[t:1.x]) (e(converged)) (0) ("nu")
				post pf (`i') (1) (`s') (8) (_b[t:M1[i]]) (_se[t:M1[i]]) (e(converged)) (0) ("omega1")
				post pf (`i') (1) (`s') (9) (_b[t:M2[i>id]]) (_se[t:M2[i>id]]) (e(converged)) (0) ("omega3")
				post pf (`i') (1) (`s') (10) (_b[t:_cons]) (_se[t:_cons]) (e(converged)) (0) ("ln_lambda")
				post pf (`i') (1) (`s') (11) (_b[/t:ln_p]) (_se[/t:ln_p]) (e(converged)) (0) ("ln_p")
				post pf (`i') (1) (`s') (12) (_b[/:var(M1[i])]) (_se[/:var(M1[i])]) (e(converged)) (0) ("sigma2_alpha")
				post pf (`i') (1) (`s') (13) (_b[/:var(M2[i>id])]) (_se[/:var(M2[i>id])]) (e(converged)) (0) ("sigma2_phi")
				post pf (`i') (1) (`s') (14) (_b[/:var(e.yobs)]) (_se[/:var(e.yobs)]) (e(converged)) (0) ("sigma2_epsilon")
				// Calculate ICCs and post that too, with SE using numerical delta method
				nlcom (_b[/:var(M1[i])] + _b[/:var(M2[i>id])]) / (_b[/:var(M1[i])] + _b[/:var(M2[i>id])] + _b[/:var(e.yobs)]), iterate(100000)
				post pf (`i') (1) (`s') (20) (r(b)[1,1]) (sqrt(r(V)[1,1])) (e(converged)) (0) ("icca")
				nlcom _b[/:var(M1[i])] / (_b[/:var(M1[i])] + _b[/:var(M2[i>id])] + _b[/:var(e.yobs)]), iterate(100000)
				post pf (`i') (1) (`s') (21) (r(b)[1,1]) (sqrt(r(V)[1,1])) (e(converged)) (0) ("iccw")
			}
		}
		nois _dots `i' 0
	}
}
postclose pf

// Compress to get a smaller file
use "simulation/data/07-full-results-ci.dta", replace
compress
save, replace

// Then, move onto the general time on treatment model

// Create postfile
capture postclose pf
postfile pf double(i model scenario parameter b stderr converged error) str20(name) using "simulation/data/07-full-results-gi.dta", replace

// Loop over scenarios and repetitions
forvalues s = 1/`scenarios' {
	nois _dots 0, title(Scenario: `s') reps(`B')
	forvalues i = 1/`B' {
		quietly {
			use "simulation/data/simdata/06-full-gi-`s'-`i'", clear
			* use "simulation/data/simdata/06-full-gi-1-1", clear
			compress
			save, replace
			// LMM
			capture gsem ///
				(yobs <- ibn.j i.cumx M1[i]@1 M2[i>id]@1, noconstant family(gaussian)) ///
				, covstructure(M1[i] M2[i>id], diag) startgrid technique(bfgs)
			// Post results
			if (_rc > 0 | e(converged) == 0) {
				// gsem failed or did not converge
				post pf (`i') (2) (`s') (1) (.) (.) (.) (99) ("beta1")
				post pf (`i') (2) (`s') (2) (.) (.) (.) (99) ("beta2")
				post pf (`i') (2) (`s') (3) (.) (.) (.) (99) ("beta3")
				post pf (`i') (2) (`s') (4) (.) (.) (.) (99) ("beta4")
				post pf (`i') (2) (`s') (5) (.) (.) (.) (99) ("beta5")
				post pf (`i') (2) (`s') (6) (.) (.) (.) (99) ("delta0")
				post pf (`i') (2) (`s') (7) (.) (.) (.) (99) ("delta1")
				post pf (`i') (2) (`s') (8) (.) (.) (.) (99) ("delta2")
				post pf (`i') (2) (`s') (9) (.) (.) (.) (99) ("delta3")
				post pf (`i') (2) (`s') (10) (.) (.) (.) (99) ("nu")
				post pf (`i') (2) (`s') (11) (.) (.) (.) (99) ("omega1")
				post pf (`i') (2) (`s') (12) (.) (.) (.) (99) ("omega3")
				post pf (`i') (2) (`s') (13) (.) (.) (.) (99) ("ln_lambda")
				post pf (`i') (2) (`s') (14) (.) (.) (.) (99) ("ln_p")
				post pf (`i') (2) (`s') (15) (.) (.) (.) (99) ("sigma2_alpha")
				post pf (`i') (2) (`s') (16) (.) (.) (.) (99) ("sigma2_phi")
				post pf (`i') (2) (`s') (17) (.) (.) (.) (99) ("sigma2_epsilon")
				// Then, parameter = 20 is ICCa, parameter = 21 is ICCw
				post pf (`i') (2) (`s') (20) (.) (.) (.) (99) ("icca")
				post pf (`i') (2) (`s') (21) (.) (.) (.) (99) ("iccw")
			}
			else {
				// gsem ok
				post pf (`i') (2) (`s') (1) (_b[yobs:1.j]) (_se[yobs:1.j]) (e(converged)) (0) ("beta1")
				post pf (`i') (2) (`s') (2) (_b[yobs:2.j]) (_se[yobs:2.j]) (e(converged)) (0) ("beta2")
				post pf (`i') (2) (`s') (3) (_b[yobs:3.j]) (_se[yobs:3.j]) (e(converged)) (0) ("beta3")
				post pf (`i') (2) (`s') (4) (_b[yobs:4.j]) (_se[yobs:4.j]) (e(converged)) (0) ("beta4")
				post pf (`i') (2) (`s') (5) (_b[yobs:5.j]) (_se[yobs:5.j]) (e(converged)) (0) ("beta5")
				post pf (`i') (2) (`s') (6) (_b[yobs:1.cumx]) (_se[yobs:1.cumx]) (e(converged)) (0) ("delta0")
				post pf (`i') (2) (`s') (7) (_b[yobs:2.cumx]) (_se[yobs:2.cumx]) (e(converged)) (0) ("delta1")
				post pf (`i') (2) (`s') (8) (_b[yobs:3.cumx]) (_se[yobs:3.cumx]) (e(converged)) (0) ("delta2")
				post pf (`i') (2) (`s') (9) (_b[yobs:4.cumx]) (_se[yobs:4.cumx]) (e(converged)) (0) ("delta3")
				post pf (`i') (2) (`s') (10) (.) (.) (e(converged)) (0) ("nu")
				post pf (`i') (2) (`s') (11) (.) (.) (e(converged)) (0) ("omega1")
				post pf (`i') (2) (`s') (12) (.) (.) (e(converged)) (0) ("omega3")
				post pf (`i') (2) (`s') (13) (.) (.) (e(converged)) (0) ("ln_lambda")
				post pf (`i') (2) (`s') (14) (.) (.) (e(converged)) (0) ("ln_p")
				post pf (`i') (2) (`s') (15) (_b[/:var(M1[i])]) (_se[/:var(M1[i])]) (e(converged)) (0) ("sigma2_alpha")
				post pf (`i') (2) (`s') (16) (_b[/:var(M2[i>id])]) (_se[/:var(M2[i>id])]) (e(converged)) (0) ("sigma2_phi")
				post pf (`i') (2) (`s') (17) (_b[/:var(e.yobs)]) (_se[/:var(e.yobs)]) (e(converged)) (0) ("sigma2_epsilon")
				// Calculate ICCs and post that too, with SE using numerical delta method
				nlcom (_b[/:var(M1[i])] + _b[/:var(M2[i>id])]) / (_b[/:var(M1[i])] + _b[/:var(M2[i>id])] + _b[/:var(e.yobs)]), iterate(100000)
				post pf (`i') (2) (`s') (20) (r(b)[1,1]) (sqrt(r(V)[1,1])) (e(converged)) (0) ("icca")
				nlcom _b[/:var(M1[i])] / (_b[/:var(M1[i])] + _b[/:var(M2[i>id])] + _b[/:var(e.yobs)]), iterate(100000)
				post pf (`i') (2) (`s') (21) (r(b)[1,1]) (sqrt(r(V)[1,1])) (e(converged)) (0) ("iccw")
			}
			// JM
			capture gsem ///
				(yobs <- ibn.j i.cumx M1[i]@1 M2[i>id]@1, noconstant family(gaussian)) ///
				(t <- i.x M1[i] M2[i>id], family(weibull, failure(d) lt(t0))) ///
				, covstructure(M1[i] M2[i>id], diag) startgrid technique(bfgs)
			// Post results
			if (_rc > 0 | e(converged) == 0) {
				// gsem failed or did not converge
				post pf (`i') (1) (`s') (1) (.) (.) (.) (99) ("beta1")
				post pf (`i') (1) (`s') (2) (.) (.) (.) (99) ("beta2")
				post pf (`i') (1) (`s') (3) (.) (.) (.) (99) ("beta3")
				post pf (`i') (1) (`s') (4) (.) (.) (.) (99) ("beta4")
				post pf (`i') (1) (`s') (5) (.) (.) (.) (99) ("beta5")
				post pf (`i') (1) (`s') (6) (.) (.) (.) (99) ("delta0")
				post pf (`i') (1) (`s') (7) (.) (.) (.) (99) ("delta1")
				post pf (`i') (1) (`s') (8) (.) (.) (.) (99) ("delta2")
				post pf (`i') (1) (`s') (9) (.) (.) (.) (99) ("delta3")
				post pf (`i') (1) (`s') (10) (.) (.) (.) (99) ("nu")
				post pf (`i') (1) (`s') (11) (.) (.) (.) (99) ("omega1")
				post pf (`i') (1) (`s') (12) (.) (.) (.) (99) ("omega3")
				post pf (`i') (1) (`s') (13) (.) (.) (.) (99) ("ln_lambda")
				post pf (`i') (1) (`s') (14) (.) (.) (.) (99) ("ln_p")
				post pf (`i') (1) (`s') (15) (.) (.) (.) (99) ("sigma2_alpha")
				post pf (`i') (1) (`s') (16) (.) (.) (.) (99) ("sigma2_phi")
				post pf (`i') (1) (`s') (17) (.) (.) (.) (99) ("sigma2_epsilon")
				// Then, parameter = 20 is ICCa, parameter = 21 is ICCw
				post pf (`i') (1) (`s') (20) (.) (.) (.) (99) ("icca")
				post pf (`i') (1) (`s') (21) (.) (.) (.) (99) ("iccw")
			}
			else {
				// gsem ok
				post pf (`i') (1) (`s') (1) (_b[yobs:1.j]) (_se[yobs:1.j]) (e(converged)) (0) ("beta1")
				post pf (`i') (1) (`s') (2) (_b[yobs:2.j]) (_se[yobs:2.j]) (e(converged)) (0) ("beta2")
				post pf (`i') (1) (`s') (3) (_b[yobs:3.j]) (_se[yobs:3.j]) (e(converged)) (0) ("beta3")
				post pf (`i') (1) (`s') (4) (_b[yobs:4.j]) (_se[yobs:4.j]) (e(converged)) (0) ("beta4")
				post pf (`i') (1) (`s') (5) (_b[yobs:5.j]) (_se[yobs:5.j]) (e(converged)) (0) ("beta5")
				post pf (`i') (1) (`s') (6) (_b[yobs:1.cumx]) (_se[yobs:1.cumx]) (e(converged)) (0) ("delta0")
				post pf (`i') (1) (`s') (7) (_b[yobs:2.cumx]) (_se[yobs:2.cumx]) (e(converged)) (0) ("delta1")
				post pf (`i') (1) (`s') (8) (_b[yobs:3.cumx]) (_se[yobs:3.cumx]) (e(converged)) (0) ("delta2")
				post pf (`i') (1) (`s') (9) (_b[yobs:4.cumx]) (_se[yobs:4.cumx]) (e(converged)) (0) ("delta3")
				post pf (`i') (1) (`s') (10) (_b[t:1.x]) (_se[t:1.x]) (e(converged)) (0) ("nu")
				post pf (`i') (1) (`s') (11) (_b[t:M1[i]]) (_se[t:M1[i]]) (e(converged)) (0) ("omega1")
				post pf (`i') (1) (`s') (12) (_b[t:M2[i>id]]) (_se[t:M2[i>id]]) (e(converged)) (0) ("omega3")
				post pf (`i') (1) (`s') (13) (_b[t:_cons]) (_se[t:_cons]) (e(converged)) (0) ("ln_lambda")
				post pf (`i') (1) (`s') (14) (_b[/t:ln_p]) (_se[/t:ln_p]) (e(converged)) (0) ("ln_p")
				post pf (`i') (1) (`s') (15) (_b[/:var(M1[i])]) (_se[/:var(M1[i])]) (e(converged)) (0) ("sigma2_alpha")
				post pf (`i') (1) (`s') (16) (_b[/:var(M2[i>id])]) (_se[/:var(M2[i>id])]) (e(converged)) (0) ("sigma2_phi")
				post pf (`i') (1) (`s') (17) (_b[/:var(e.yobs)]) (_se[/:var(e.yobs)]) (e(converged)) (0) ("sigma2_epsilon")
				// Calculate ICCs and post that too, with SE using numerical delta method
				nlcom (_b[/:var(M1[i])] + _b[/:var(M2[i>id])]) / (_b[/:var(M1[i])] + _b[/:var(M2[i>id])] + _b[/:var(e.yobs)]), iterate(100000)
				post pf (`i') (1) (`s') (20) (r(b)[1,1]) (sqrt(r(V)[1,1])) (e(converged)) (0) ("icca")
				nlcom _b[/:var(M1[i])] / (_b[/:var(M1[i])] + _b[/:var(M2[i>id])] + _b[/:var(e.yobs)]), iterate(100000)
				post pf (`i') (1) (`s') (21) (r(b)[1,1]) (sqrt(r(V)[1,1])) (e(converged)) (0) ("iccw")
			}
		}
		nois _dots `i' 0
	}
}
postclose pf

// Compress to get a smaller file
use "simulation/data/07-full-results-gi.dta", replace
compress
save, replace
