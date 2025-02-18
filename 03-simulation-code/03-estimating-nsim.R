# With this script, we summarise results from 50 iterations to estimate how many iterations to run
# This is based on minimising Monte Carlo error for bias, across simulation scenarios

library(tidyverse)
library(haven)
library(rsimsum)
library(formattable)
library(glue)
library(ggrepel)
library(simswjm)

### Constant intervention model

# Preliminary data
simres_ci <- read_dta(file = "simulation/data/02-preliminary-results-ci.dta") %>%
  zap_label() %>%
  zap_labels() %>%
  zap_formats() |>
  mutate(converged = ifelse(error != 0, 0, converged))
# Those with converged = 0 did not converge (according to Stata's gsem checks, or captured errors)
simres_ci %>%
  distinct(scenario, i, model, converged) %>%
  group_by(scenario, model) %>%
  summarise(sum(converged), mean(converged)) |>
  View(title = "Convergence, CI")
#
scenarios_wide <- readRDS(file = "simulation/data/01-scenarios.RDS")
scenarios <- scenarios_wide |>
  pivot_longer(cols = -scenario) |>
  rename(true = value)
#
simres_ci <- left_join(simres_ci, scenarios, by = c("scenario", "name")) |>
  mutate(model = factor(model, levels = seq(2), labels = c("JM", "LMM")))
# Then, run multisimsum for each component of the model separately
# We base our calculations on the coefficients for the longitudinal part, which are the main estimands of interest for this study
simres_ci_summdf <- map_dfr(.x = levels(simres_ci$model), .f = function(m) {
  s <- rsimsum::multisimsum(data = filter(simres_ci, error == 0 & converged == 1 & model == m & grepl("delta", name)), par = "name", estvarname = "b", se = "stderr", by = "scenario", true = "true")
  s <- tidy(summary(s))
  s[["model"]] <- m
  return(s)
}) |>
  mutate(scenario = as.numeric(as.character(scenario))) |>
  left_join(scenarios_wide, by = "scenario")

### General time on treatment model

# Preliminary data
simres_gi <- read_dta(file = "simulation/data/02-preliminary-results-gi.dta") %>%
  zap_label() %>%
  zap_labels() %>%
  zap_formats() |>
  mutate(converged = ifelse(error != 0, 0, converged))
# Those with converged = 0 did not converge (according to Stata's gsem checks, or captured errors)
simres_gi %>%
  distinct(scenario, i, model, converged) %>%
  group_by(scenario, model) %>%
  summarise(sum(converged), mean(converged)) |>
  View(title = "Convergence, GI")
#
simres_gi <- left_join(simres_gi, scenarios, by = c("scenario", "name")) |>
  mutate(model = factor(model, levels = seq(2), labels = c("JM", "LMM")))
# Then, run multisimsum for each component of the model separately
# We base our calculations on the treatment effect coefficients (for the longitudinal part)
simres_gi_summdf <- map_dfr(.x = levels(simres_gi$model), .f = function(m) {
  s <- rsimsum::multisimsum(data = filter(simres_gi, error == 0 & converged == 1 & model == m & grepl("delta", name)), par = "name", estvarname = "b", se = "stderr", by = "scenario", true = "true")
  s <- tidy(summary(s))
  s[["model"]] <- m
  return(s)
}) |>
  mutate(scenario = as.numeric(as.character(scenario))) |>
  left_join(scenarios_wide, by = "scenario")

### Largest EmpSE, ModelSE
mSE_ci <- max(simres_ci_summdf$est[simres_ci_summdf$stat %in% c("empse", "modelse")], na.rm = TRUE)
mVar_ci <- mSE_ci^2
mSE_gi <- max(simres_gi_summdf$est[simres_gi_summdf$stat %in% c("empse", "modelse")], na.rm = TRUE)
mVar_gi <- mSE_gi^2
largest_mVar <- max(c(mVar_ci, mVar_gi))
# Expected MCSE with 500 repetitions:
expected_mcse_500 <- sqrt(largest_mVar / 500)
# Expected MCSE with 1000 repetitions:
expected_mcse_1000 <- sqrt(largest_mVar / 1000)
# Expected MCSE with 1500 repetitions:
expected_mcse_1500 <- sqrt(largest_mVar / 1500)
# Expected ratio to delta values
expected_ratio <- tibble(
  delta = unique(scenarios$true[grepl("delta", scenarios$name)]),
  ratio_500 = expected_mcse_500 / delta,
  ratio_1000 = expected_mcse_1000 / delta,
  ratio_1500 = expected_mcse_1500 / delta
) |>
  arrange(delta)
expected_ratio

### How many repetitions?
sink(file = "simulation/data/03-mcse.md")
cat("CONSTANT INTERVENTION MODEL\n\n")
cat("Based on a preliminary, largest (Model SE, Empirical SE) of:\n")
print(sqrt(mVar_ci))
cat("... i.e., variance of:\n")
print(mVar_ci)
cat("GENERAL TIME ON TREATMENT MODEL\n\n")
cat("Based on a preliminary, largest (Model SE, Empirical SE) of:\n")
print(sqrt(mVar_gi))
cat("... i.e., variance of:\n")
print(mVar_gi)
cat("COMBINED LARGEST mVAR:\n")
print(largest_mVar)
cat("EXPECTED MCSE WITH 500 REPETITIONS:\n")
print(expected_mcse_500)
cat("EXPECTED MCSE WITH 1000 REPETITIONS:\n")
print(expected_mcse_1000)
cat("EXPECTED MCSE WITH 1500 REPETITIONS:\n")
print(expected_mcse_1500)
cat("EXPECTED MCSE / DELTA RATIO:\n")
print(expected_ratio)
sink()
