# With this script, we run B iterations of the simulation study.
# More details in the ms.
library(tidyverse)
library(glue)
library(simswjm)
library(haven)

### Seed
set.seed(898968904)

# Scenarios:
scenarios_neutral <- readRDS(file = "simulation/data/01-scenarios.RDS") |>
  filter(scenario == 9) |>
  mutate(Intercept = -1.5)
saveRDS(object = scenarios_neutral, file = "simulation/data/06-scenarios-neutral.RDS")

### Constant Intervention Model

# Simulate B iterations per scenario, N patients per cluster
B <- 1000
pb <- txtProgressBar(max = B * nrow(scenarios_neutral), style = 3)
for (scn in seq(nrow(scenarios_neutral))) {
  for (i in seq(B)) {
    file_name <- glue("simulation/data/simdata/06-full-ci-{scn}-{i}.dta")
    tmp <- swtrial_inf(
      repn = i,
      k = scenarios_neutral$N[scn],
      omega1 = scenarios_neutral$omega1[scn], # alpha
      omega2 = scenarios_neutral$omega2[scn], # gamma
      omega3 = scenarios_neutral$omega3[scn], # phi
      deltas = rep(scenarios_neutral$delta[scn], 4),
      betas = c(
        scenarios_neutral$beta1[scn],
        scenarios_neutral$beta2[scn],
        scenarios_neutral$beta3[scn],
        scenarios_neutral$beta4[scn],
        scenarios_neutral$beta5[scn]
      ),
      nu = scenarios_neutral$nu[scn],
      i = scenarios_neutral$i[scn],
      j = scenarios_neutral$j[scn],
      intervention_seq = scenarios_neutral$intervention_seq[scn],
      sigma_alpha = sqrt(scenarios_neutral$sigma2_alpha[scn]),
      sigma_gamma = sqrt(scenarios_neutral$sigma2_gamma[scn]),
      sigma_phi = sqrt(scenarios_neutral$sigma2_phi[scn]),
      sigma_epsilon = sqrt(scenarios_neutral$sigma2_epsilon[scn]),
      logistic_dropout = TRUE,
      logistic_intercept = scenarios_neutral$Intercept[scn],
      family = "gaussian"
    )
    write_dta(data = tmp, path = file_name)
    setTxtProgressBar(pb = pb, value = pb$getVal() + 1)
  }
}
close(pb)

### General Time on Treatment Model

# Simulate B iterations per scenario, N patients per cluster
pb <- txtProgressBar(max = B * nrow(scenarios_neutral), style = 3)
for (scn in seq(nrow(scenarios_neutral))) {
  for (i in seq(B)) {
    file_name <- glue("simulation/data/simdata/06-full-gi-{scn}-{i}.dta")
    tmp <- swtrial_inf(
      repn = i,
      k = scenarios_neutral$N[scn],
      omega1 = scenarios_neutral$omega1[scn], # alpha
      omega2 = scenarios_neutral$omega2[scn], # gamma
      omega3 = scenarios_neutral$omega3[scn], # phi
      deltas = c(
        scenarios_neutral$delta0[scn],
        scenarios_neutral$delta1[scn],
        scenarios_neutral$delta2[scn],
        scenarios_neutral$delta3[scn]
      ),
      betas = c(
        scenarios_neutral$beta1[scn],
        scenarios_neutral$beta2[scn],
        scenarios_neutral$beta3[scn],
        scenarios_neutral$beta4[scn],
        scenarios_neutral$beta5[scn]
      ),
      nu = scenarios_neutral$nu[scn],
      i = scenarios_neutral$i[scn],
      j = scenarios_neutral$j[scn],
      intervention_seq = scenarios_neutral$intervention_seq[scn],
      sigma_alpha = sqrt(scenarios_neutral$sigma2_alpha[scn]),
      sigma_gamma = sqrt(scenarios_neutral$sigma2_gamma[scn]),
      sigma_phi = sqrt(scenarios_neutral$sigma2_phi[scn]),
      sigma_epsilon = sqrt(scenarios_neutral$sigma2_epsilon[scn]),
      logistic_dropout = TRUE,
      logistic_intercept = scenarios_neutral$Intercept[scn],
      family = "gaussian"
    )
    write_dta(data = tmp, path = file_name)
    setTxtProgressBar(pb = pb, value = pb$getVal() + 1)
  }
}
close(pb)
