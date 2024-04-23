This folder includes example code for an individually-randomised stepped wedge trial with a continuous outcome and informative dropout.

[`02a-data-simulation.R`](02a-data-simulation.R) includes the R code required to simulate the data, with information on the R session in [`02a-data-simulation-sessioninfo.txt`](02a-data-simulation-sessioninfo.txt).

The simulated data in included as [`02-dt-ci.dta`](02-dt-ci.dta) (for the constant intervention parametrisation) and [`02-dt-gi.dta`](02-dt-gi.dta) (for the general time on treatment parametrisation).

[`02b-data-analysis.do`](02b-data-analysis.do) includes the Stata code required to fit the linear mixed model and joint model under both the constant intervention and generalised time on treatment parametrisations.

Finally, a Stata log with the results of the analysis is included in [`02b-data-analysis-log.txt`](02b-data-analysis-log.txt).
