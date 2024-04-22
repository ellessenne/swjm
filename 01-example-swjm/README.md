This folder includes example code for a stepped wedge trial with a continuous outcome and informative dropout.

[`01a-data-simulation.R`](01a-data-simulation.R) includes the R code required to simulate the data, with information on the R session in [`01a-data-simulation-sessioninfo.txt`](01a-data-simulation-sessioninfo.txt).

The simulated data in included as [`01-dt-ci.dta`](01-dt-ci.dta) (for the constant intervention parametrisation) and [`01-dt-gi.dta`](01-dt-gi.dta) (for the general time on treatment parametrisation).

[`01b-data-analysis.do`](01b-data-analysis.do) includes the Stata code required to fit the linear mixed model and joint model under both the constant intervention and generalised time on treatment parametrisations.

Finally, a Stata log in included in [`01b-data-analysis-log.txt`](01b-data-analysis-log.txt).
