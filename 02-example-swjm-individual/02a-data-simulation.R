### Required packages
library(simswjm)
library(haven)
library(glue)
library(sessioninfo)
library(here)
library(tidyverse)

### Seed, for reproducibility
set.seed(945837)

### Data-generating parameters
scenarios <- crossing(
  k = 100,
  tribble(
    ~omega1, ~omega2, ~omega3, # omega1 is for alpha, omega2 is for gamma, omega3 is for phi
    0.0, 0.0, log(1.5),
  ),
  delta = 20,
  nu = -1.0,
  ln_lambda = -1.5,
  ln_p = 0.0,
  sigma2_alpha = 0.0,
  sigma2_gamma = 0.0,
  sigma2_phi = 55.0,
  sigma2_epsilon = 40.0,
  beta1 = 30, beta2 = 30, beta3 = 30, beta4 = 30, beta5 = 30
) |>
  mutate(delta0 = delta * 0.0, delta1 = delta * 0.5, delta2 = delta * 1.0, delta3 = delta * 1.5)

### Constant Intervention Model
dt_ci <- swtrial_inf(
  repn = 1,
  k = scenarios$k,
  omega1 = scenarios$omega1, # alpha
  omega2 = scenarios$omega2, # gamma
  omega3 = scenarios$omega3, # phi
  deltas = rep(scenarios$delta, 4),
  betas = c(
    scenarios$beta1,
    scenarios$beta2,
    scenarios$beta3,
    scenarios$beta4,
    scenarios$beta5
  ),
  nu = scenarios$nu,
  lambda = exp(scenarios$ln_lambda),
  p = exp(scenarios$ln_p),
  i = 4 * 8,
  j = 5,
  intervention_seq = 4,
  sigma_alpha = sqrt(scenarios$sigma2_alpha),
  sigma_gamma = sqrt(scenarios$sigma2_gamma),
  sigma_phi = sqrt(scenarios$sigma2_phi),
  sigma_epsilon = sqrt(scenarios$sigma2_epsilon),
  family = "gaussian"
)

### General Time on Treatment Model
dt_gi <- swtrial_inf(
  repn = 1,
  k = scenarios$k,
  omega1 = scenarios$omega1, # alpha
  omega2 = scenarios$omega2, # gamma
  omega3 = scenarios$omega3, # phi
  deltas = c(
    scenarios$delta0,
    scenarios$delta1,
    scenarios$delta2,
    scenarios$delta3
  ),
  betas = c(
    scenarios$beta1,
    scenarios$beta2,
    scenarios$beta3,
    scenarios$beta4,
    scenarios$beta5
  ),
  nu = scenarios$nu,
  lambda = exp(scenarios$ln_lambda),
  p = exp(scenarios$ln_p),
  i = 4 * 8,
  j = 5,
  intervention_seq = 4,
  sigma_alpha = sqrt(scenarios$sigma2_alpha),
  sigma_gamma = sqrt(scenarios$sigma2_gamma),
  sigma_phi = sqrt(scenarios$sigma2_phi),
  sigma_epsilon = sqrt(scenarios$sigma2_epsilon),
  family = "gaussian"
)

### Export datasets
write_dta(data = dt_ci, path = here("02-example-swjm-individual/02-dt-ci.dta"))
write_dta(data = dt_gi, path = here("02-example-swjm-individual/02-dt-gi.dta"))

### Export session info
sink(file = here("02-example-swjm-individual/02a-data-simulation-sessioninfo.txt"))
cat("Last run:\n")
Sys.time()
cat("\n")
session_info()
sink()
