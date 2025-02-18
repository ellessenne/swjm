### With this script, we run a small number of iterations per scenario to estimate the
#   total number of replications for the full simulation study
library(tidyverse)
library(glue)
library(simswjm)
library(haven)

### Seed
set.seed(9462831)

### Constant Intervention Model

# Scenarios:
scenarios <- crossing(
  N = 50,
  tribble(
    ~omega1, ~omega2, ~omega3,
    log(0.9), 0.0, log(0.9),
    0.0, 0.0, 0.0,
    log(2.0), 0.0, log(2.0),
    log(0.5), 0.0, log(0.5)
  ),
  delta = c(0.0, 5, 25),
  nu = c(-0.2, 0.0),
  ln_lambda = -1.5,
  ln_p = 0.0,
  tribble(
    ~sigma2_alpha, ~sigma2_gamma, ~sigma2_phi, ~sigma2_epsilon,
    2.0, 0.0, 55.0, 40.0,
  ),
  j = 5, # periods
  beta1 = 30, beta2 = 30, beta3 = 30, beta4 = 30, beta5 = 30, # period effects
  intervention_seq = 4, # number of intervention sequences
  i = 4 * 8 # number of clusters (8) per intervention_seq (4)
) |>
  mutate(delta0 = delta * 0.0, delta1 = delta * 0.5, delta2 = delta * 1.0, delta3 = delta * 1.25)
scenarios <- bind_rows(
  scenarios,
  scenarios |> filter(nu == -0.2 & delta == 5 & omega1 == log(0.9) & omega3 == log(0.9)) |> mutate(i = 4 * 3),
  scenarios |> filter(nu == -0.2 & delta == 5 & omega1 == log(0.9) & omega3 == log(0.9)) |> mutate(i = 4 * 3, N = 100),
  scenarios |> filter(nu == -0.2 & delta == 5 & omega1 == log(0.9) & omega3 == log(0.9)) |> mutate(sigma2_alpha = sigma2_alpha / 2, sigma2_gamma = sigma2_gamma / 2, sigma2_phi = sigma2_phi / 2),
  scenarios |> filter(nu == -0.2 & delta == 5 & omega1 == log(0.9) & omega3 == log(0.9)) |> mutate(sigma2_alpha = sigma2_alpha * 2, sigma2_gamma = sigma2_gamma * 2, sigma2_phi = sigma2_phi * 2)
) |>
  mutate(scenario = row_number())
scenarios$icca <- map_dbl(.x = 1:nrow(scenarios), .f = function(i) icc(sigma_alpha = sqrt(scenarios$sigma2_alpha[i]), sigma_gamma = sqrt(scenarios$sigma2_gamma[i]), sigma_phi = sqrt(scenarios$sigma2_phi[i]), sigma_epsilon = sqrt(scenarios$sigma2_epsilon[i]))$ICC_a)
scenarios$iccw <- map_dbl(.x = 1:nrow(scenarios), .f = function(i) icc(sigma_alpha = sqrt(scenarios$sigma2_alpha[i]), sigma_gamma = sqrt(scenarios$sigma2_gamma[i]), sigma_phi = sqrt(scenarios$sigma2_phi[i]), sigma_epsilon = sqrt(scenarios$sigma2_epsilon[i]))$ICC_w)
saveRDS(object = scenarios, file = "simulation/data/01-scenarios.RDS")

# Simulate B iterations per scenario, N patients per cluster
B <- 50
pb <- txtProgressBar(max = B * nrow(scenarios), style = 3)
for (scn in seq(nrow(scenarios))) {
  for (i in seq(B)) {
    file_name <- glue("simulation/data/simdata/01-preliminary-ci-{scn}-{i}.dta")
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
    file_name <- glue("simulation/data/simdata/01-preliminary-gi-{scn}-{i}.dta")
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
