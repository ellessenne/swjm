# With this script, we run B iterations of the simulation study.
# More details in the ms.
library(tidyverse)
library(glue)
library(simswjm)
library(haven)

### Seed
set.seed(28734678)

# Scenarios:
scenarios <- readRDS(file = "simulation/data/01-scenarios.RDS")

### Constant Intervention Model

# Simulate B iterations per scenario, N patients per cluster
B <- 1000
pb <- txtProgressBar(max = B * nrow(scenarios), style = 3)
for (scn in seq(nrow(scenarios))) {
  for (i in seq(B)) {
    file_name <- glue("simulation/data/simdata/04-full-ci-{scn}-{i}.dta")
    tmp <- swtrial_inf(
      repn = i,
      k = scenarios$N[scn],
      omega1 = scenarios$omega1[scn], # alpha
      omega2 = scenarios$omega2[scn], # gamma
      omega3 = scenarios$omega3[scn], # phi
      deltas = rep(scenarios$delta[scn], 4),
      betas = c(
        scenarios$beta1[scn],
        scenarios$beta2[scn],
        scenarios$beta3[scn],
        scenarios$beta4[scn],
        scenarios$beta5[scn]
      ),
      nu = scenarios$nu[scn],
      lambda = exp(scenarios$ln_lambda[scn]),
      p = exp(scenarios$ln_p[scn]),
      i = scenarios$i[scn],
      j = scenarios$j[scn],
      intervention_seq = scenarios$intervention_seq[scn],
      sigma_alpha = sqrt(scenarios$sigma2_alpha[scn]),
      sigma_gamma = sqrt(scenarios$sigma2_gamma[scn]),
      sigma_phi = sqrt(scenarios$sigma2_phi[scn]),
      sigma_epsilon = sqrt(scenarios$sigma2_epsilon[scn]),
      family = "gaussian"
    )
    write_dta(data = tmp, path = file_name)
    setTxtProgressBar(pb = pb, value = pb$getVal() + 1)
  }
}
close(pb)

### General Time on Treatment Model

# Simulate B iterations per scenario, N patients per cluster
pb <- txtProgressBar(max = B * nrow(scenarios), style = 3)
for (scn in seq(nrow(scenarios))) {
  for (i in seq(B)) {
    file_name <- glue("simulation/data/simdata/04-full-gi-{scn}-{i}.dta")
    tmp <- swtrial_inf(
      repn = i,
      k = scenarios$N[scn],
      omega1 = scenarios$omega1[scn], # alpha
      omega2 = scenarios$omega2[scn], # gamma
      omega3 = scenarios$omega3[scn], # phi
      deltas = c(
        scenarios$delta0[scn],
        scenarios$delta1[scn],
        scenarios$delta2[scn],
        scenarios$delta3[scn]
      ),
      betas = c(
        scenarios$beta1[scn],
        scenarios$beta2[scn],
        scenarios$beta3[scn],
        scenarios$beta4[scn],
        scenarios$beta5[scn]
      ),
      nu = scenarios$nu[scn],
      lambda = exp(scenarios$ln_lambda[scn]),
      p = exp(scenarios$ln_p[scn]),
      i = scenarios$i[scn],
      j = scenarios$j[scn],
      intervention_seq = scenarios$intervention_seq[scn],
      sigma_alpha = sqrt(scenarios$sigma2_alpha[scn]),
      sigma_gamma = sqrt(scenarios$sigma2_gamma[scn]),
      sigma_phi = sqrt(scenarios$sigma2_phi[scn]),
      sigma_epsilon = sqrt(scenarios$sigma2_epsilon[scn]),
      family = "gaussian"
    )
    write_dta(data = tmp, path = file_name)
    setTxtProgressBar(pb = pb, value = pb$getVal() + 1)
  }
}
close(pb)
