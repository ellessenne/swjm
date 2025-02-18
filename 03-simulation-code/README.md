This folder includes code in R and Stata that could be used to replicate the Monte Carlo simulations reported in the manuscript.
Specifically:

1. The `01-data-simulation-preliminary.R` script is used to simulate data for 50 preliminary repetitions of the simulation study;

1. `02-data-analysis-preliminary.do` includes the Stata code used to fit models using the preliminary repetitions;

1. Then, `03-estimating-nsim.R` estimates the required number of repetitions to control Monte Carlo standard errors within an acceptable limit based on the preliminary repetitions of the simulation study;

1. `04-data-simulation-ms.R` and `05-data-analysis-ms.do` are used to simulate data and fit models for the entire simulation study, for the data-generating mechanisms based on the joint model;

1. `06-data-simulation-neutral.R` and `07-data-analysis-neutral.do` are used to simulate data and fit models for the neutral simulation scenarios;

1. Finally, `10-simulation-summary-ms.R` and `11-simulation-summary-neutral.R` are used to summarise the results of all simulations and produce tables and plots that are included in the manuscript and supplementary material.
