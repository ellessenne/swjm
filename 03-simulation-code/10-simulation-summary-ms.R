# With this script, we summarise results from the full, 1000 + 50 repetitions of the simulation study.

library(tidyverse)
library(haven)
library(rsimsum)
library(formattable)
library(glue)
library(ggrepel)
library(simswjm)
library(conflicted)
library(here)
library(kableExtra)
library(survival)
library(KMunicate)

#
options(scipen = 100)

# Conflicts
conflicts_prefer(dplyr::filter)

# Custom scales, theme
scale_color_rda <- scale_color_manual(values = c("#152D49", "#E72D3F"))
scale_fill_rda <- scale_fill_manual(values = c("#152D49", "#E72D3F"))
theme_rda <- theme_bw(base_size = 12, base_family = "Atkinson Hyperlegible")

# Read in data
simres_ci_prel <- read_dta(file = here("simulation/data/02-preliminary-results-ci.dta")) |>
  mutate(trt = "ci")
simres_ci <- read_dta(file = here("simulation/data/05-full-results-ci.dta")) |>
  mutate(i = i + max(simres_ci_prel$i), trt = "ci")
simres_gi_prel <- read_dta(file = here("simulation/data/02-preliminary-results-gi.dta")) |>
  mutate(trt = "gi")
simres_gi <- read_dta(file = "simulation/data/05-full-results-gi.dta") |>
  mutate(i = i + max(simres_gi_prel$i), trt = "gi")

# Combine
simres <- bind_rows(simres_ci_prel, simres_ci, simres_gi_prel, simres_gi) |>
  zap_label() %>%
  zap_labels() %>%
  zap_formats() |>
  mutate(converged = ifelse(error != 0, 0, converged))
# Those with converged = 0 did not converge (according to Stata's gsem checks, or captured errors)
simres %>%
  distinct(trt, scenario, i, model, converged) %>%
  group_by(trt, scenario, model) %>%
  summarise(sum(converged), mean(converged)) |>
  View(title = "Convergence")

# Read-in scenarios and true values
scenarios_wide <- readRDS(file = here("simulation/data/01-scenarios.RDS"))
scenarios_labels <- scenarios_wide |>
  select(scenario) |>
  mutate(scenario_lbl = factor(
    scenario,
    levels = 1:28,
    labels = c(
      "1: list(list(omega[1], omega[2]) == log(0.5), nu == -0.2, sigma[alpha]^2 == 2, sigma[phi]^2 == 55, sigma[epsilon]^2 == 40, N == 50, i == 32, delta == 0)",
      "2: list(list(omega[1], omega[2]) == log(0.5), nu ==  0.0, sigma[alpha]^2 == 2, sigma[phi]^2 == 55, sigma[epsilon]^2 == 40, N == 50, i == 32, delta == 0)",
      "3: list(list(omega[1], omega[2]) == log(0.5), nu == -0.2, sigma[alpha]^2 == 2, sigma[phi]^2 == 55, sigma[epsilon]^2 == 40, N == 50, i == 32, delta == 5)",
      "4: list(list(omega[1], omega[2]) == log(0.5), nu ==  0.0, sigma[alpha]^2 == 2, sigma[phi]^2 == 55, sigma[epsilon]^2 == 40, N == 50, i == 32, delta == 5)",
      "5: list(list(omega[1], omega[2]) == log(0.5), nu == -0.2, sigma[alpha]^2 == 2, sigma[phi]^2 == 55, sigma[epsilon]^2 == 40, N == 50, i == 32, delta == 25)",
      "6: list(list(omega[1], omega[2]) == log(0.5), nu ==  0.0, sigma[alpha]^2 == 2, sigma[phi]^2 == 55, sigma[epsilon]^2 == 40, N == 50, i == 32, delta == 25)",
      "7: list(list(omega[1], omega[2]) == log(0.9), nu == -0.2, sigma[alpha]^2 == 2, sigma[phi]^2 == 55, sigma[epsilon]^2 == 40, N == 50, i == 32, delta == 0)",
      "8: list(list(omega[1], omega[2]) == log(0.9), nu ==  0.0, sigma[alpha]^2 == 2, sigma[phi]^2 == 55, sigma[epsilon]^2 == 40, N == 50, i == 32, delta == 0)",
      "9: list(list(omega[1], omega[2]) == log(0.9), nu == -0.2, sigma[alpha]^2 == 2, sigma[phi]^2 == 55, sigma[epsilon]^2 == 40, N == 50, i == 32, delta == 5)",
      "10: list(list(omega[1], omega[2]) == log(0.9), nu ==  0.0, sigma[alpha]^2 == 2, sigma[phi]^2 == 55, sigma[epsilon]^2 == 40, N == 50, i == 32, delta == 5)",
      "11: list(list(omega[1], omega[2]) == log(0.9), nu == -0.2, sigma[alpha]^2 == 2, sigma[phi]^2 == 55, sigma[epsilon]^2 == 40, N == 50, i == 32, delta == 25)",
      "12: list(list(omega[1], omega[2]) == log(0.9), nu ==  0.0, sigma[alpha]^2 == 2, sigma[phi]^2 == 55, sigma[epsilon]^2 == 40, N == 50, i == 32, delta == 25)",
      "13: list(list(omega[1], omega[2]) == log(1.0), nu == -0.2, sigma[alpha]^2 == 2, sigma[phi]^2 == 55, sigma[epsilon]^2 == 40, N == 50, i == 32, delta == 0)",
      "14: list(list(omega[1], omega[2]) == log(1.0), nu ==  0.0, sigma[alpha]^2 == 2, sigma[phi]^2 == 55, sigma[epsilon]^2 == 40, N == 50, i == 32, delta == 0)",
      "15: list(list(omega[1], omega[2]) == log(1.0), nu == -0.2, sigma[alpha]^2 == 2, sigma[phi]^2 == 55, sigma[epsilon]^2 == 40, N == 50, i == 32, delta == 5)",
      "16: list(list(omega[1], omega[2]) == log(1.0), nu ==  0.0, sigma[alpha]^2 == 2, sigma[phi]^2 == 55, sigma[epsilon]^2 == 40, N == 50, i == 32, delta == 5)",
      "17: list(list(omega[1], omega[2]) == log(1.0), nu == -0.2, sigma[alpha]^2 == 2, sigma[phi]^2 == 55, sigma[epsilon]^2 == 40, N == 50, i == 32, delta == 25)",
      "18: list(list(omega[1], omega[2]) == log(1.0), nu ==  0.0, sigma[alpha]^2 == 2, sigma[phi]^2 == 55, sigma[epsilon]^2 == 40, N == 50, i == 32, delta == 25)",
      "19: list(list(omega[1], omega[2]) == log(2.0), nu == -0.2, sigma[alpha]^2 == 2, sigma[phi]^2 == 55, sigma[epsilon]^2 == 40, N == 50, i == 32, delta == 0)",
      "20: list(list(omega[1], omega[2]) == log(2.0), nu ==  0.0, sigma[alpha]^2 == 2, sigma[phi]^2 == 55, sigma[epsilon]^2 == 40, N == 50, i == 32, delta == 0)",
      "21: list(list(omega[1], omega[2]) == log(2.0), nu == -0.2, sigma[alpha]^2 == 2, sigma[phi]^2 == 55, sigma[epsilon]^2 == 40, N == 50, i == 32, delta == 5)",
      "22: list(list(omega[1], omega[2]) == log(2.0), nu ==  0.0, sigma[alpha]^2 == 2, sigma[phi]^2 == 55, sigma[epsilon]^2 == 40, N == 50, i == 32, delta == 5)",
      "23: list(list(omega[1], omega[2]) == log(2.0), nu == -0.2, sigma[alpha]^2 == 2, sigma[phi]^2 == 55, sigma[epsilon]^2 == 40, N == 50, i == 32, delta == 25)",
      "24: list(list(omega[1], omega[2]) == log(2.0), nu ==  0.0, sigma[alpha]^2 == 2, sigma[phi]^2 == 55, sigma[epsilon]^2 == 40, N == 50, i == 32, delta == 25)",
      "25: list(list(omega[1], omega[2]) == log(0.9), nu == -0.2, sigma[alpha]^2 == 2, sigma[phi]^2 == 55, sigma[epsilon]^2 == 40, N == 50, i == 12, delta == 5)",
      "26: list(list(omega[1], omega[2]) == log(0.9), nu == -0.2, sigma[alpha]^2 == 2, sigma[phi]^2 == 55, sigma[epsilon]^2 == 40, N == 100, i == 12, delta == 5)",
      "27: list(list(omega[1], omega[2]) == log(0.9), nu == -0.2, sigma[alpha]^2 == 1, sigma[phi]^2 == 27.5, sigma[epsilon]^2 == 40, N == 50, i == 32, delta == 5)",
      "28: list(list(omega[1], omega[2]) == log(0.9), nu == -0.2, sigma[alpha]^2 == 4, sigma[phi]^2 == 110, sigma[epsilon]^2 == 40, N == 50, i == 32, delta == 5)"
    )
  )) |>
  mutate(scenario_tex = factor(
    scenario,
    levels = 1:28,
    labels = c(
      "1: $\\omega_1 = \\omega_2 = \\log(0.5), \\nu = -0.2, \\sigma_{\\alpha}^2 = 2, \\sigma_{\\phi}^2 = 55, \\sigma_{\\varepsilon}^2 = 40, N = 50, i = 32, \\delta = 0$",
      "2: $\\omega_1 = \\omega_2 = \\log(0.5), \\nu =  0.0, \\sigma_{\\alpha}^2 = 2, \\sigma_{\\phi}^2 = 55, \\sigma_{\\varepsilon}^2 = 40, N = 50, i = 32, \\delta = 0$",
      "3: $\\omega_1 = \\omega_2 = \\log(0.5), \\nu = -0.2, \\sigma_{\\alpha}^2 = 2, \\sigma_{\\phi}^2 = 55, \\sigma_{\\varepsilon}^2 = 40, N = 50, i = 32, \\delta = 5$",
      "4: $\\omega_1 = \\omega_2 = \\log(0.5), \\nu =  0.0, \\sigma_{\\alpha}^2 = 2, \\sigma_{\\phi}^2 = 55, \\sigma_{\\varepsilon}^2 = 40, N = 50, i = 32, \\delta = 5$",
      "5: $\\omega_1 = \\omega_2 = \\log(0.5), \\nu = -0.2, \\sigma_{\\alpha}^2 = 2, \\sigma_{\\phi}^2 = 55, \\sigma_{\\varepsilon}^2 = 40, N = 50, i = 32, \\delta = 25$",
      "6: $\\omega_1 = \\omega_2 = \\log(0.5), \\nu =  0.0, \\sigma_{\\alpha}^2 = 2, \\sigma_{\\phi}^2 = 55, \\sigma_{\\varepsilon}^2 = 40, N = 50, i = 32, \\delta = 25$",
      "7: $\\omega_1 = \\omega_2 = \\log(0.9), \\nu = -0.2, \\sigma_{\\alpha}^2 = 2, \\sigma_{\\phi}^2 = 55, \\sigma_{\\varepsilon}^2 = 40, N = 50, i = 32, \\delta = 0$",
      "8: $\\omega_1 = \\omega_2 = \\log(0.9), \\nu =  0.0, \\sigma_{\\alpha}^2 = 2, \\sigma_{\\phi}^2 = 55, \\sigma_{\\varepsilon}^2 = 40, N = 50, i = 32, \\delta = 0$",
      "9: $\\omega_1 = \\omega_2 = \\log(0.9), \\nu = -0.2, \\sigma_{\\alpha}^2 = 2, \\sigma_{\\phi}^2 = 55, \\sigma_{\\varepsilon}^2 = 40, N = 50, i = 32, \\delta = 5$",
      "10: $\\omega_1 = \\omega_2 = \\log(0.9), \\nu =  0.0, \\sigma_{\\alpha}^2 = 2, \\sigma_{\\phi}^2 = 55, \\sigma_{\\varepsilon}^2 = 40, N = 50, i = 32, \\delta = 5$",
      "11: $\\omega_1 = \\omega_2 = \\log(0.9), \\nu = -0.2, \\sigma_{\\alpha}^2 = 2, \\sigma_{\\phi}^2 = 55, \\sigma_{\\varepsilon}^2 = 40, N = 50, i = 32, \\delta = 25$",
      "12: $\\omega_1 = \\omega_2 = \\log(0.9), \\nu =  0.0, \\sigma_{\\alpha}^2 = 2, \\sigma_{\\phi}^2 = 55, \\sigma_{\\varepsilon}^2 = 40, N = 50, i = 32, \\delta = 25$",
      "13: $\\omega_1 = \\omega_2 = \\log(1.0), \\nu = -0.2, \\sigma_{\\alpha}^2 = 2, \\sigma_{\\phi}^2 = 55, \\sigma_{\\varepsilon}^2 = 40, N = 50, i = 32, \\delta = 0$",
      "14: $\\omega_1 = \\omega_2 = \\log(1.0), \\nu =  0.0, \\sigma_{\\alpha}^2 = 2, \\sigma_{\\phi}^2 = 55, \\sigma_{\\varepsilon}^2 = 40, N = 50, i = 32, \\delta = 0$",
      "15: $\\omega_1 = \\omega_2 = \\log(1.0), \\nu = -0.2, \\sigma_{\\alpha}^2 = 2, \\sigma_{\\phi}^2 = 55, \\sigma_{\\varepsilon}^2 = 40, N = 50, i = 32, \\delta = 5$",
      "16: $\\omega_1 = \\omega_2 = \\log(1.0), \\nu =  0.0, \\sigma_{\\alpha}^2 = 2, \\sigma_{\\phi}^2 = 55, \\sigma_{\\varepsilon}^2 = 40, N = 50, i = 32, \\delta = 5$",
      "17: $\\omega_1 = \\omega_2 = \\log(1.0), \\nu = -0.2, \\sigma_{\\alpha}^2 = 2, \\sigma_{\\phi}^2 = 55, \\sigma_{\\varepsilon}^2 = 40, N = 50, i = 32, \\delta = 25$",
      "18: $\\omega_1 = \\omega_2 = \\log(1.0), \\nu =  0.0, \\sigma_{\\alpha}^2 = 2, \\sigma_{\\phi}^2 = 55, \\sigma_{\\varepsilon}^2 = 40, N = 50, i = 32, \\delta = 25$",
      "19: $\\omega_1 = \\omega_2 = \\log(2.0), \\nu = -0.2, \\sigma_{\\alpha}^2 = 2, \\sigma_{\\phi}^2 = 55, \\sigma_{\\varepsilon}^2 = 40, N = 50, i = 32, \\delta = 0$",
      "20: $\\omega_1 = \\omega_2 = \\log(2.0), \\nu =  0.0, \\sigma_{\\alpha}^2 = 2, \\sigma_{\\phi}^2 = 55, \\sigma_{\\varepsilon}^2 = 40, N = 50, i = 32, \\delta = 0$",
      "21: $\\omega_1 = \\omega_2 = \\log(2.0), \\nu = -0.2, \\sigma_{\\alpha}^2 = 2, \\sigma_{\\phi}^2 = 55, \\sigma_{\\varepsilon}^2 = 40, N = 50, i = 32, \\delta = 5$",
      "22: $\\omega_1 = \\omega_2 = \\log(2.0), \\nu =  0.0, \\sigma_{\\alpha}^2 = 2, \\sigma_{\\phi}^2 = 55, \\sigma_{\\varepsilon}^2 = 40, N = 50, i = 32, \\delta = 5$",
      "23: $\\omega_1 = \\omega_2 = \\log(2.0), \\nu = -0.2, \\sigma_{\\alpha}^2 = 2, \\sigma_{\\phi}^2 = 55, \\sigma_{\\varepsilon}^2 = 40, N = 50, i = 32, \\delta = 25$",
      "24: $\\omega_1 = \\omega_2 = \\log(2.0), \\nu =  0.0, \\sigma_{\\alpha}^2 = 2, \\sigma_{\\phi}^2 = 55, \\sigma_{\\varepsilon}^2 = 40, N = 50, i = 32, \\delta = 25$",
      "25: $\\omega_1 = \\omega_2 = \\log(0.9), \\nu = -0.2, \\sigma_{\\alpha}^2 = 2, \\sigma_{\\phi}^2 = 55, \\sigma_{\\varepsilon}^2 = 40, N = 50, i = 12, \\delta = 5$",
      "26: $\\omega_1 = \\omega_2 = \\log(0.9), \\nu = -0.2, \\sigma_{\\alpha}^2 = 2, \\sigma_{\\phi}^2 = 55, \\sigma_{\\varepsilon}^2 = 40, N = 100, i = 12, \\delta = 5$",
      "27: $\\omega_1 = \\omega_2 = \\log(0.9), \\nu = -0.2, \\sigma_{\\alpha}^2 = 1, \\sigma_{\\phi}^2 = 27.5, \\sigma_{\\varepsilon}^2 = 40, N = 50, i = 32, \\delta = 5$",
      "28: $\\omega_1 = \\omega_2 = \\log(0.9), \\nu = -0.2, \\sigma_{\\alpha}^2 = 4, \\sigma_{\\phi}^2 = 110, \\sigma_{\\varepsilon}^2 = 40, N = 50, i = 32, \\delta = 5$"
    )
  )) |>
  mutate(scenario_extra_lbl = factor(
    scenario,
    levels = c(9, 25:28),
    labels = c(
      "Ref. ~ scenario",
      "i == 3 %*% 4",
      "list(i == 3 %*% 4, N == 100)",
      "list(sigma[alpha]^2 == 1, sigma[phi]^2 == 27.5)",
      "list(sigma[alpha]^2 == 4, sigma[phi]^2 == 110)"
    )
  )) |>
  mutate(scenario_extra_tex = factor(
    scenario,
    levels = c(9, 25:28),
    labels = c(
      "Reference scenario",
      "$i = 3 \\times 4$",
      "$i = 3 \\times 4, N = 100$",
      "$\\sigma_{\\alpha}^2 = 1, \\sigma_{\\phi}^2 = 27.5$",
      "$\\sigma_{\\alpha}^2 = 4, \\sigma_{\\phi}^2 = 110$"
    )
  ))
scenarios <- scenarios_wide |>
  pivot_longer(cols = -scenario) |>
  rename(true = value)

# Combine
simres <- left_join(simres, scenarios, by = c("scenario", "name"))

# Study non-convergence further
simres_convergence <- simres %>%
  distinct(trt, scenario, i, model, converged) %>%
  group_by(trt, scenario, model) %>%
  summarise(nc = sum(converged), pc = mean(converged), pnc = 1 - pc) %>%
  ungroup() |>
  left_join(scenarios_wide, by = "scenario") |>
  mutate(omega = omega1) %>%
  mutate(omega = round(omega, 3), nu = round(nu, 3)) |>
  mutate(omega = factor(omega, levels = c(-0.693, -0.105, 0.000, 0.693), labels = c("list(omega[1], omega[2]) == log(0.5)", "list(omega[1], omega[2]) == log(0.9)", "list(omega[1], omega[2]) == log(1.0)", "list(omega[1], omega[2]) == log(2.0)"))) |>
  mutate(delta = factor(delta, levels = c(0, 5, 25), labels = c("delta == 0.0", "delta == 5.0", "delta == 25"))) %>%
  mutate(nu = factor(nu, levels = c(-0.2, 0.0), labels = c("nu == -0.2", "nu == 0.0"))) %>%
  mutate(model = factor(model, levels = c(2, 1), labels = c("Linear Mixed Model", "Joint Model"))) |>
  mutate(lbl = glue("{formattable::percent(pnc, 1)}")) |>
  left_join(scenarios_labels)
# Plots:
# CI:
p_conv_ci <- ggplot(filter(simres_convergence, trt == "ci"), aes(x = fct_rev(scenario_lbl), y = pnc, color = model, fill = model)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  geom_col(position = position_dodge(width = 4 / 5), width = 0.5) +
  geom_text(aes(label = lbl, family = "Atkinson Hyperlegible"), position = position_dodge(width = 4 / 5), color = "grey25", hjust = -0.1, size = 2.5) +
  scale_color_rda +
  scale_fill_rda +
  theme_rda +
  theme(legend.position = "top", axis.text.y = element_text(hjust = 0)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip(ylim = c(0, 1)) +
  labs(x = "", y = "Proportion of Non-Converged Iterations", color = "", fill = "")
p_conv_ci
ggsave(p_conv_ci, filename = here("figures/p_conv_ci.pdf"), device = cairo_pdf, width = 9, height = 7)
# GI:
p_conv_gi <- ggplot(filter(simres_convergence, trt == "gi"), aes(x = fct_rev(scenario_lbl), y = pnc, color = model, fill = model)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  geom_col(position = position_dodge(width = 4 / 5), width = 0.5) +
  geom_text(aes(label = lbl, family = "Atkinson Hyperlegible"), position = position_dodge(width = 4 / 5), color = "grey25", hjust = -0.1, size = 2.5) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  scale_fill_rda +
  theme_rda +
  theme(legend.position = "top", axis.text.y = element_text(hjust = 0)) +
  coord_flip(ylim = c(0, 1)) +
  labs(x = "", y = "Proportion of Non-Converged Iterations", color = "", fill = "")
p_conv_gi
ggsave(p_conv_gi, filename = here("figures/p_conv_gi.pdf"), device = cairo_pdf, width = 9, height = 7)

# Cleanup
rm(p_conv_ci, p_conv_gi, simres_convergence)
gc()

# Then, run simsum for each model and parameter
simres_summdf <- map_dfr(.x = unique(simres$name), .f = function(p) {
  map_dfr(.x = unique(simres$model), .f = function(m) {
    dd <- filter(simres, error == 0 & converged == 1 & name == p & model == m)
    s <- simsum(data = dd, estvarname = "b", se = "stderr", by = c("trt", "scenario"), true = "true")
    s <- tidy(summary(s))
    s[["model"]] <- m
    s[["name"]] <- p
    return(s)
  })
}, .progress = TRUE) |>
  mutate(scenario = as.numeric(as.character(scenario))) |>
  left_join(scenarios_wide, by = "scenario")

# Create labels, e.g., for plots and tables
simres_summdf <- simres_summdf |>
  mutate(omega = omega1) %>%
  mutate(omega = round(omega, 3)) |>
  mutate(omega_tex = factor(omega, levels = c(-0.693, -0.105, 0.000, 0.693), labels = c("$\\omega_1 = \\omega_2 = \\log(0.5)$", "$\\omega_1 = \\omega_2 = \\log(0.9)$", "$\\omega_1 = \\omega_2 = \\log(1.0)$", "$\\omega_1 = \\omega_2 = \\log(2.0)$"))) %>%
  mutate(omega = factor(omega, levels = c(-0.693, -0.105, 0.000, 0.693), labels = c("list(omega[1], omega[2]) == log(0.5)", "list(omega[1], omega[2]) == log(0.9)", "list(omega[1], omega[2]) == log(1.0)", "list(omega[1], omega[2]) == log(2.0)"))) |>
  mutate(delta_tex = factor(delta, levels = c(0, 5, 25), labels = c("$\\delta = 0.0$", "$\\delta = 5.0$", "$\\delta = 25$"))) %>%
  mutate(delta_gi = factor(delta, levels = c(0, 5, 25), labels = c("delta[0] == delta[1] == delta[2] == delta[3] == 0.00", "list(delta[0] == 0.00, delta[1] == 2.50, delta[2] == 5.00, delta[3] == 6.25)", "list(delta[0] == 0.00, delta[1] == 12.50, delta[2] == 25.00, delta[3] == 31.25)"))) |>
  mutate(delta_gi_tex = factor(delta, levels = c(0, 5, 25), labels = c("$\\delta_0 = \\delta_1 = \\delta_2 = \\delta_3 = 0.00$", "$\\delta_0 = 0.00, \\delta_1 = 2.50, \\delta_2 = 5.00, \\delta_3 = 6.25$", "$\\delta_0 = 0.00, \\delta_1 = 12.50, \\delta_2 = 25.00, \\delta_3 = 31.25$"))) |>
  mutate(delta = factor(delta, levels = c(0, 5, 25), labels = c("delta == 0.0", "delta == 5.0", "delta == 25"))) %>%
  mutate(nu_tex = factor(nu, levels = c(-0.2, 0.0), labels = c("$\\nu = -0.2$", "$\\nu = 0.0$"))) %>%
  mutate(nu = factor(nu, levels = c(-0.2, 0.0), labels = c("nu == -0.2", "nu == 0.0"))) %>%
  mutate(modelshort = factor(model, levels = c(2, 1), labels = c("LMM", "JM"))) |>
  mutate(model = factor(model, levels = c(2, 1), labels = c("Linear Mixed Model", "Joint Model"))) |>
  mutate(name_label = case_when(
    name == "beta1" ~ "beta[1]",
    name == "beta2" ~ "beta[2]",
    name == "beta3" ~ "beta[3]",
    name == "beta4" ~ "beta[4]",
    name == "beta5" ~ "beta[5]",
    name == "omega1" ~ "omega[1]",
    name == "omega3" ~ "omega[3]",
    name == "ln_lambda" ~ "log(lambda)",
    name == "ln_p" ~ "log(p)",
    name == "sigma2_alpha" ~ "sigma[alpha]^2",
    name == "sigma2_phi" ~ "sigma[phi]^2",
    name == "sigma2_epsilon" ~ "sigma[epsilon]^2",
    name == "icca" ~ "rho[a]",
    name == "iccw" ~ "rho[d]",
    name == "delta0" ~ "delta[0]",
    name == "delta1" ~ "delta[1]",
    name == "delta2" ~ "delta[2]",
    name == "delta3" ~ "delta[3]",
    TRUE ~ name
  )) |>
  mutate(name_tex = case_when(
    name == "beta1" ~ "$\\beta_1$",
    name == "beta2" ~ "$\\beta_2$",
    name == "beta3" ~ "$\\beta_3$",
    name == "beta4" ~ "$\\beta_4$",
    name == "beta5" ~ "$\\beta_5$",
    name == "omega1" ~ "$\\omega_1$",
    name == "omega3" ~ "$\\omega_3$",
    name == "ln_lambda" ~ "$\\log(\\lambda)$",
    name == "ln_p" ~ "$\\log(p)$",
    name == "sigma2_alpha" ~ "$\\sigma_{\\alpha}^2$",
    name == "sigma2_phi" ~ "$\\sigma_{\\phi}^2$",
    name == "sigma2_epsilon" ~ "$\\sigma_{\\epsilon}^2$",
    name == "icca" ~ "$\\rho_a$",
    name == "iccw" ~ "$\\rho_d$",
    name == "delta" ~ "$\\delta$",
    name == "delta0" ~ "$\\delta_0$",
    name == "delta1" ~ "$\\delta_1$",
    name == "delta2" ~ "$\\delta_2$",
    name == "delta3" ~ "$\\delta_3$",
    name == "nu" ~ "$\\nu$"
  )) |>
  left_join(scenarios_labels)
# Export simres_summdf, could be useful
saveRDS(object = simres_summdf, file = here("simulation/data/10-simres_summdf.RDS"))

# Plot bias of longitudinal treatment effect(s)
p_bias_ci_delta <- simres_summdf %>%
  filter(stat == "bias" & name == "delta" & trt == "ci" & scenario <= 24) |>
  mutate(omega = fct_rev(factor(omega))) |>
  ggplot(aes(x = omega, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_rda +
  theme_rda +
  facet_grid(delta ~ nu, labeller = label_parsed) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Bias (95% C.I.)", color = "")
p_bias_ci_delta
ggsave(p_bias_ci_delta, filename = here("figures/p_bias_ci_delta.pdf"), device = cairo_pdf, width = 7, height = 4)
#
p_bias_ci_delta_extra <- simres_summdf %>%
  filter(stat == "bias" & name == "delta" & trt == "ci" & (scenario > 24 | scenario == 9)) |>
  ggplot(aes(x = fct_rev(scenario_extra_lbl), y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_rda +
  theme_rda +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Bias (95% C.I.)", color = "")
p_bias_ci_delta_extra
ggsave(p_bias_ci_delta_extra, filename = here("figures/p_bias_ci_delta_extra.pdf"), device = cairo_pdf, width = 5, height = 3)
#
p_bias_gi_delta <- simres_summdf %>%
  filter(stat == "bias" & grepl("delta", name) & trt == "gi" & scenario <= 24) |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model, shape = name_label)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 2 / 3)) +
  geom_point(position = position_dodge(width = 2 / 3)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_rda +
  theme_rda +
  facet_grid(delta + nu ~ omega, labeller = label_parsed) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Bias (95% C.I.)", color = "") +
  guides(shape = "none")
p_bias_gi_delta
ggsave(p_bias_gi_delta, filename = here("figures/p_bias_gi_delta.pdf"), device = cairo_pdf, width = 7, height = 7)
#
p_bias_gi_delta_extra <- simres_summdf %>%
  filter(stat == "bias" & grepl("delta", name) & trt == "gi" & (scenario > 24 | scenario == 9)) |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model, shape = name_label)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 2 / 3)) +
  geom_point(position = position_dodge(width = 2 / 3)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_extra_lbl, labeller = label_parsed, nrow = 2) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Bias (95% C.I.)", color = "") +
  guides(shape = "none")
p_bias_gi_delta_extra
ggsave(p_bias_gi_delta_extra, filename = here("figures/p_bias_gi_delta_extra.pdf"), device = cairo_pdf, width = 6, height = 4)

# Plot relative bias of longitudinal treatment effect(s)
p_rbias_ci_delta <- simres_summdf %>%
  filter(stat == "rbias" & name == "delta" & trt == "ci" & scenario <= 24) |>
  mutate(omega = fct_rev(factor(omega))) |>
  ggplot(aes(x = omega, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_grid(delta ~ nu, labeller = label_parsed) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Relative Bias (95% C.I.)", color = "")
p_rbias_ci_delta
ggsave(p_rbias_ci_delta, filename = here("figures/p_rbias_ci_delta.pdf"), device = cairo_pdf, width = 7, height = 4)
#
p_rbias_ci_delta_extra <- simres_summdf %>%
  filter(stat == "rbias" & name == "delta" & trt == "ci" & (scenario > 24 | scenario == 9)) |>
  ggplot(aes(x = fct_rev(scenario_extra_lbl), y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Relative Bias (95% C.I.)", color = "")
p_rbias_ci_delta_extra
ggsave(p_rbias_ci_delta_extra, filename = here("figures/p_rbias_ci_delta_extra.pdf"), device = cairo_pdf, width = 5, height = 3)
#
p_rbias_gi_delta <- simres_summdf %>%
  filter(stat == "rbias" & grepl("delta", name) & trt == "gi" & scenario <= 24) |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model, shape = name_label)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 2 / 3)) +
  geom_point(position = position_dodge(width = 2 / 3)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_grid(delta + nu ~ omega, labeller = label_parsed) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Relative Bias (95% C.I.)", color = "") +
  guides(shape = "none")
p_rbias_gi_delta
ggsave(p_rbias_gi_delta, filename = here("figures/p_rbias_gi_delta.pdf"), device = cairo_pdf, width = 7, height = 7)
#
p_rbias_gi_delta_extra <- simres_summdf %>%
  filter(stat == "rbias" & grepl("delta", name) & trt == "gi" & (scenario > 24 | scenario == 9)) |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model, shape = name_label)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 2 / 3)) +
  geom_point(position = position_dodge(width = 2 / 3)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_extra_lbl, labeller = label_parsed, nrow = 2) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Relative Bias (95% C.I.)", color = "") +
  guides(shape = "none")
p_rbias_gi_delta_extra
ggsave(p_rbias_gi_delta_extra, filename = here("figures/p_rbias_gi_delta_extra.pdf"), device = cairo_pdf, width = 6, height = 4)

# Table for the same quantities
#
sink(file = glue("tables/tbias_ci_delta.tex"))
simres_summdf %>%
  filter(stat == "bias" & name == "delta" & trt == "ci") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(Scenario = glue("{omega_tex}, {delta_tex}, {nu_tex}")) |>
  mutate(Scenario = ifelse(scenario > 24, glue("{scenario_extra_tex}"), Scenario)) |>
  select(Scenario, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    linesep = "",
    caption = "Bias of treatment effect on the longitudinal outcome for the constant treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. \\label{tab:tbias_ci_delta}",
    escape = FALSE
  ) |>
  kable_styling(latex_options = c("repeat_header")) |>
  pack_rows("Main scenarios:", 1, 24) %>%
  pack_rows("Additional scenarios:", 25, 28) |>
  print()
sink()
#
sink(file = glue("tables/trbias_ci_delta.tex"))
simres_summdf %>%
  filter(stat == "rbias" & name == "delta" & trt == "ci") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(y = ifelse(!is.finite(est), "---", y)) |>
  mutate(Scenario = glue("{omega_tex}, {delta_tex}, {nu_tex}")) |>
  mutate(Scenario = ifelse(scenario > 24, glue("{scenario_extra_tex}"), Scenario)) |>
  select(Scenario, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    linesep = "",
    caption = "Relative bias of treatment effect on the longitudinal outcome for the constant treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. \\label{tab:trbias_ci_delta}",
    escape = FALSE
  ) |>
  kable_styling(latex_options = c("repeat_header")) |>
  pack_rows("Main scenarios:", 1, 24) %>%
  pack_rows("Additional scenarios:", 25, 28) |>
  print()
sink()
#
sink(file = glue("tables/tbias_gi_delta.tex"))
simres_summdf %>%
  filter(stat == "bias" & grepl("delta", name) & trt == "gi") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(Scenario = glue("{omega_tex}, {delta_gi_tex}, {nu_tex}")) |>
  mutate(Scenario = ifelse(scenario > 24, glue("{scenario_extra_tex}"), Scenario)) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    linesep = "",
    caption = "Bias of treatment effect on the longitudinal outcome for the general time on treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. \\label{tab:tbias_gi_delta}",
    escape = FALSE
  ) |>
  kable_styling(latex_options = c("repeat_header")) |>
  pack_rows("Main scenarios:", 1, 24 * 4) %>%
  pack_rows("Additional scenarios:", 24 * 4 + 1, 28 * 4) |>
  landscape() |>
  print()
sink()
#
sink(file = glue("tables/trbias_gi_delta.tex"))
simres_summdf %>%
  filter(stat == "rbias" & grepl("delta", name) & trt == "gi") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(y = ifelse(!is.finite(est), "---", y)) |>
  mutate(Scenario = glue("{omega_tex}, {delta_gi_tex}, {nu_tex}")) |>
  mutate(Scenario = ifelse(scenario > 24, glue("{scenario_extra_tex}"), Scenario)) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    linesep = "",
    caption = "Relative bias of treatment effect on the longitudinal outcome for the general time on treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. \\label{tab:trbias_gi_delta}",
    escape = FALSE
  ) |>
  kable_styling(latex_options = c("repeat_header")) |>
  pack_rows("Main scenarios:", 1, 24 * 4) %>%
  pack_rows("Additional scenarios:", 24 * 4 + 1, 28 * 4) |>
  landscape() |>
  print()
sink()

# Plot bias for the period effects
p_bias_ci_beta <- simres_summdf %>%
  filter(stat == "bias" & grepl("beta", name) & trt == "ci" & scenario <= 24) |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  geom_point() +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_rda +
  theme_rda +
  facet_grid(delta + nu ~ omega, labeller = label_parsed) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Bias (95% C.I.)", color = "")
p_bias_ci_beta
ggsave(p_bias_ci_beta, filename = here("figures/p_bias_ci_beta.pdf"), device = cairo_pdf, width = 7, height = 7)
#
p_bias_ci_beta_extra <- simres_summdf %>%
  filter(stat == "bias" & grepl("beta", name) & trt == "ci" & (scenario > 24 | scenario == 9)) |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  geom_point() +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_extra_lbl, labeller = label_parsed, nrow = 2) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Bias (95% C.I.)", color = "")
p_bias_ci_beta_extra
ggsave(p_bias_ci_beta_extra, filename = here("figures/p_bias_ci_beta_extra.pdf"), device = cairo_pdf, width = 6, height = 4)
#
p_bias_gi_beta <- simres_summdf %>%
  filter(stat == "bias" & grepl("beta", name) & trt == "gi" & scenario <= 24) |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  geom_point() +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_rda +
  theme_rda +
  facet_grid(delta + nu ~ omega, labeller = label_parsed) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Bias (95% C.I.)", color = "")
p_bias_gi_beta
ggsave(p_bias_gi_beta, filename = here("figures/p_bias_gi_beta.pdf"), device = cairo_pdf, width = 7, height = 7)
#
p_bias_gi_beta_extra <- simres_summdf %>%
  filter(stat == "bias" & grepl("beta", name) & trt == "gi" & (scenario > 24 | scenario == 9)) |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  geom_point() +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_extra_lbl, labeller = label_parsed, nrow = 2) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Bias (95% C.I.)", color = "")
p_bias_gi_beta_extra
ggsave(p_bias_gi_beta_extra, filename = here("figures/p_bias_gi_beta_extra.pdf"), device = cairo_pdf, width = 6, height = 4)

# Plot relative bias for the period effects
p_rbias_ci_beta <- simres_summdf %>%
  filter(stat == "rbias" & grepl("beta", name) & trt == "ci" & scenario <= 24) |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  geom_point() +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_grid(delta + nu ~ omega, labeller = label_parsed) +
  theme(legend.position = "bottom") +
  coord_flip(ylim = c(-0.3, 0.3)) +
  labs(x = "", y = "Relative Bias (95% C.I.)", color = "")
p_rbias_ci_beta
ggsave(p_rbias_ci_beta, filename = here("figures/p_rbias_ci_beta.pdf"), device = cairo_pdf, width = 7, height = 7)
#
p_rbias_ci_beta_extra <- simres_summdf %>%
  filter(stat == "rbias" & grepl("beta", name) & trt == "ci" & (scenario > 24 | scenario == 9)) |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  geom_point() +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_extra_lbl, labeller = label_parsed, nrow = 2) +
  theme(legend.position = "bottom") +
  coord_flip(ylim = c(-0.3, 0.3)) +
  labs(x = "", y = "Relative Bias (95% C.I.)", color = "")
p_rbias_ci_beta_extra
ggsave(p_rbias_ci_beta_extra, filename = here("figures/p_rbias_ci_beta_extra.pdf"), device = cairo_pdf, width = 6, height = 4)
#
p_rbias_gi_beta <- simres_summdf %>%
  filter(stat == "rbias" & grepl("beta", name) & trt == "gi" & scenario <= 24) |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  geom_point() +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_grid(delta + nu ~ omega, labeller = label_parsed) +
  theme(legend.position = "bottom") +
  coord_flip(ylim = c(-0.3, 0.3)) +
  labs(x = "", y = "Bias (95% C.I.)", color = "")
p_rbias_gi_beta
ggsave(p_rbias_gi_beta, filename = here("figures/p_rbias_gi_beta.pdf"), device = cairo_pdf, width = 7, height = 7)
#
p_rbias_gi_beta_extra <- simres_summdf %>%
  filter(stat == "rbias" & grepl("beta", name) & trt == "gi" & (scenario > 24 | scenario == 9)) |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  geom_point() +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_extra_lbl, labeller = label_parsed, nrow = 2) +
  theme(legend.position = "bottom") +
  coord_flip(ylim = c(-0.3, 0.3)) +
  labs(x = "", y = "Bias (95% C.I.)", color = "")
p_rbias_gi_beta_extra
ggsave(p_rbias_gi_beta_extra, filename = here("figures/p_rbias_gi_beta_extra.pdf"), device = cairo_pdf, width = 6, height = 4)

# Table for the same quantities
#
sink(file = glue("tables/tbias_ci_beta.tex"))
simres_summdf %>%
  filter(stat == "bias" & grepl("beta", name) & trt == "ci") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(Scenario = glue("{omega_tex}, {delta_tex}, {nu_tex}")) |>
  mutate(Scenario = ifelse(scenario > 24, glue("{scenario_extra_tex}"), Scenario)) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    linesep = "",
    caption = "Bias of period effects on the longitudinal outcome for the constant treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. \\label{tab:tbias_ci_beta}",
    escape = FALSE
  ) |>
  kable_styling(latex_options = c("repeat_header")) |>
  pack_rows("Main scenarios:", 1, 24 * 5) %>%
  pack_rows("Additional scenarios:", 24 * 5 + 1, 28 * 5) |>
  print()
sink()
#
sink(file = glue("tables/tbias_gi_beta.tex"))
simres_summdf %>%
  filter(stat == "bias" & grepl("beta", name) & trt == "gi") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(Scenario = glue("{omega_tex}, {delta_tex}, {nu_tex}")) |>
  mutate(Scenario = ifelse(scenario > 24, glue("{scenario_extra_tex}"), Scenario)) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    linesep = "",
    caption = "Bias of period effects on the longitudinal outcome for the general time on treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. \\label{tab:tbias_gi_beta}",
    escape = FALSE
  ) |>
  kable_styling(latex_options = c("repeat_header")) |>
  pack_rows("Main scenarios:", 1, 24 * 5) %>%
  pack_rows("Additional scenarios:", 24 * 5 + 1, 28 * 5) |>
  landscape() |>
  print()
sink()
#
sink(file = glue("tables/trbias_ci_beta.tex"))
simres_summdf %>%
  filter(stat == "rbias" & grepl("beta", name) & trt == "ci") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(Scenario = glue("{omega_tex}, {delta_tex}, {nu_tex}")) |>
  mutate(Scenario = ifelse(scenario > 24, glue("{scenario_extra_tex}"), Scenario)) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    linesep = "",
    caption = "Relative bias of period effects on the longitudinal outcome for the constant treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. \\label{tab:trbias_ci_beta}",
    escape = FALSE
  ) |>
  kable_styling(latex_options = c("repeat_header")) |>
  pack_rows("Main scenarios:", 1, 24 * 5) %>%
  pack_rows("Additional scenarios:", 24 * 5 + 1, 28 * 5) |>
  print()
sink()
#
sink(file = glue("tables/trbias_gi_beta.tex"))
simres_summdf %>%
  filter(stat == "rbias" & grepl("beta", name) & trt == "gi") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(Scenario = glue("{omega_tex}, {delta_tex}, {nu_tex}")) |>
  mutate(Scenario = ifelse(scenario > 24, glue("{scenario_extra_tex}"), Scenario)) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    linesep = "",
    caption = "Relative bias of period effects on the longitudinal outcome for the general time on treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. \\label{tab:trbias_gi_beta}",
    escape = FALSE
  ) |>
  kable_styling(latex_options = c("repeat_header")) |>
  pack_rows("Main scenarios:", 1, 24 * 5) %>%
  pack_rows("Additional scenarios:", 24 * 5 + 1, 28 * 5) |>
  landscape() |>
  print()
sink()

# Plot coverage probability of longitudinal treatment effect(s)
p_cover_ci_delta <- simres_summdf %>%
  filter(stat == "cover" & name == "delta" & trt == "ci" & scenario <= 24) |>
  mutate(omega = fct_rev(factor(omega))) |>
  ggplot(aes(x = omega, y = est, color = model)) +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_grid(delta ~ nu, labeller = label_parsed) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Coverage Probability (95% C.I.)", color = "")
p_cover_ci_delta
ggsave(p_cover_ci_delta, filename = here("figures/p_cover_ci_delta.pdf"), device = cairo_pdf, width = 7, height = 4)
#
p_cover_ci_delta_extra <- simres_summdf %>%
  filter(stat == "cover" & name == "delta" & trt == "ci" & (scenario > 24 | scenario == 9)) |>
  ggplot(aes(x = fct_rev(scenario_extra_lbl), y = est, color = model)) +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Coverage Probability (95% C.I.)", color = "")
p_cover_ci_delta_extra
ggsave(p_cover_ci_delta_extra, filename = here("figures/p_cover_ci_delta_extra.pdf"), device = cairo_pdf, width = 5, height = 3)
#
p_cover_gi_delta <- simres_summdf %>%
  filter(stat == "cover" & grepl("delta", name) & trt == "gi" & scenario <= 24) |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model, shape = name_label)) +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 2 / 3)) +
  geom_point(position = position_dodge(width = 2 / 3)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_grid(delta + nu ~ omega, labeller = label_parsed) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Coverage Probability (95% C.I.)", color = "") +
  guides(shape = "none")
p_cover_gi_delta
ggsave(p_cover_gi_delta, filename = here("figures/p_cover_gi_delta.pdf"), device = cairo_pdf, width = 7, height = 7)
#
p_cover_gi_delta_extra <- simres_summdf %>%
  filter(stat == "cover" & grepl("delta", name) & trt == "gi" & (scenario > 24 | scenario == 9)) |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model, shape = name_label)) +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 2 / 3)) +
  geom_point(position = position_dodge(width = 2 / 3)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_extra_lbl, labeller = label_parsed, nrow = 2) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Coverage Probability (95% C.I.)", color = "") +
  guides(shape = "none")
p_cover_gi_delta_extra
ggsave(p_cover_gi_delta_extra, filename = here("figures/p_cover_gi_delta_extra.pdf"), device = cairo_pdf, width = 6, height = 4)

# Table for the same quantities
#
sink(file = glue("tables/tcover_ci_delta.tex"))
simres_summdf %>%
  filter(stat == "cover" & name == "delta" & trt == "ci") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(Scenario = glue("{omega_tex}, {delta_tex}, {nu_tex}")) |>
  mutate(Scenario = ifelse(scenario > 24, glue("{scenario_extra_tex}"), Scenario)) |>
  select(Scenario, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    linesep = "",
    caption = "Coverage probability of treatment effect on the longitudinal outcome for the constant treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. \\label{tab:tcover_ci_delta}",
    escape = FALSE
  ) |>
  kable_styling(latex_options = c("repeat_header")) |>
  pack_rows("Main scenarios:", 1, 24) %>%
  pack_rows("Additional scenarios:", 25, 28) |>
  print()
sink()
#
sink(file = glue("tables/tcover_gi_delta.tex"))
simres_summdf %>%
  filter(stat == "cover" & grepl("delta", name) & trt == "gi") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(Scenario = glue("{omega_tex}, {delta_gi_tex}, {nu_tex}")) |>
  mutate(Scenario = ifelse(scenario > 24, glue("{scenario_extra_tex}"), Scenario)) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    linesep = "",
    caption = "Coverage probability of treatment effect on the longitudinal outcome for the general time on treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. \\label{tab:tcover_gi_delta}",
    escape = FALSE
  ) |>
  kable_styling(latex_options = c("repeat_header")) |>
  pack_rows("Main scenarios:", 1, 24 * 4) %>%
  pack_rows("Additional scenarios:", 24 * 4 + 1, 28 * 4) |>
  landscape() |>
  print()
sink()
#

# Plot coverage probability of period effects(s)
p_cover_ci_beta <- simres_summdf %>%
  filter(stat == "cover" & grepl("beta", name) & trt == "ci" & scenario <= 24) |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_grid(delta + nu ~ omega, labeller = label_parsed) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Coverage Probability (95% C.I.)", color = "")
p_cover_ci_beta
ggsave(p_cover_ci_beta, filename = here("figures/p_cover_ci_beta.pdf"), device = cairo_pdf, width = 7, height = 7)
#
p_cover_ci_beta_extra <- simres_summdf %>%
  filter(stat == "cover" & grepl("beta", name) & trt == "ci" & (scenario > 24 | scenario == 9)) |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_extra_lbl, labeller = label_parsed, nrow = 2) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Coverage Probability (95% C.I.)", color = "")
p_cover_ci_beta_extra
ggsave(p_cover_ci_beta_extra, filename = here("figures/p_cover_ci_beta_extra.pdf"), device = cairo_pdf, width = 6, height = 4)
#
p_cover_gi_beta <- simres_summdf %>%
  filter(stat == "cover" & grepl("beta", name) & trt == "gi" & scenario <= 24) |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_grid(delta + nu ~ omega, labeller = label_parsed) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Coverage Probability (95% C.I.)", color = "")
p_cover_gi_beta
ggsave(p_cover_gi_beta, filename = here("figures/p_cover_gi_beta.pdf"), device = cairo_pdf, width = 7, height = 7)
#
p_cover_gi_beta_extra <- simres_summdf %>%
  filter(stat == "cover" & grepl("beta", name) & trt == "gi" & (scenario > 24 | scenario == 9)) |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_extra_lbl, labeller = label_parsed, nrow = 2) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Coverage Probability (95% C.I.)", color = "")
p_cover_gi_beta_extra
ggsave(p_cover_gi_beta_extra, filename = here("figures/p_cover_gi_beta_extra.pdf"), device = cairo_pdf, width = 6, height = 4)

# Table for the same quantities
#
sink(file = glue("tables/tcover_ci_beta.tex"))
simres_summdf %>%
  filter(stat == "cover" & grepl("beta", name) & trt == "ci") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(Scenario = glue("{omega_tex}, {delta_tex}, {nu_tex}")) |>
  mutate(Scenario = ifelse(scenario > 24, glue("{scenario_extra_tex}"), Scenario)) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    linesep = "",
    caption = "Coverage probability of period effects on the longitudinal outcome for the constant treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. \\label{tab:tcover_ci_beta}",
    escape = FALSE
  ) |>
  kable_styling(latex_options = c("repeat_header")) |>
  pack_rows("Main scenarios:", 1, 24 * 5) %>%
  pack_rows("Additional scenarios:", 24 * 5 + 1, 28 * 5) |>
  print()
sink()
#
sink(file = glue("tables/tcover_gi_beta.tex"))
simres_summdf %>%
  filter(stat == "cover" & grepl("beta", name) & trt == "gi") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(Scenario = glue("{omega_tex}, {delta_gi_tex}, {nu_tex}")) |>
  mutate(Scenario = ifelse(scenario > 24, glue("{scenario_extra_tex}"), Scenario)) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    linesep = "",
    caption = "Coverage probability of period effects on the longitudinal outcome for the general time on treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. \\label{tab:tcover_gi_beta}",
    escape = FALSE
  ) |>
  kable_styling(latex_options = c("repeat_header")) |>
  pack_rows("Main scenarios:", 1, 24 * 5) %>%
  pack_rows("Additional scenarios:", 24 * 5, 28 * 5) |>
  landscape() |>
  print()
sink()

# Bias of variance components
p_bias_ci_var <- simres_summdf %>%
  filter(stat == "bias" & grepl("sigma", name) & trt == "ci" & scenario <= 24) |>
  mutate(name_label = fct_rev(factor(name_label, levels = c("sigma[alpha]^2", "sigma[phi]^2", "sigma[epsilon]^2")))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_rda +
  theme_rda +
  facet_grid(delta + nu ~ omega, labeller = label_parsed) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Bias (95% C.I.)", color = "")
p_bias_ci_var
ggsave(p_bias_ci_var, filename = here("figures/p_bias_ci_var.pdf"), device = cairo_pdf, width = 7, height = 7)
#
p_bias_ci_var_extra <- simres_summdf %>%
  filter(stat == "bias" & grepl("sigma", name) & trt == "ci" & (scenario > 24 | scenario == 9)) |>
  mutate(name_label = fct_rev(factor(name_label, levels = c("sigma[alpha]^2", "sigma[phi]^2", "sigma[epsilon]^2")))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_extra_lbl, labeller = label_parsed, nrow = 2) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Bias (95% C.I.)", color = "")
p_bias_ci_var_extra
ggsave(p_bias_ci_var_extra, filename = here("figures/p_bias_ci_var_extra.pdf"), device = cairo_pdf, width = 6, height = 4)
#
p_bias_gi_var <- simres_summdf %>%
  filter(stat == "bias" & grepl("sigma", name) & trt == "gi" & scenario <= 24) |>
  mutate(name_label = fct_rev(factor(name_label, levels = c("sigma[alpha]^2", "sigma[phi]^2", "sigma[epsilon]^2")))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_rda +
  theme_rda +
  facet_grid(delta + nu ~ omega, labeller = label_parsed) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Bias (95% C.I.)", color = "")
p_bias_gi_var
ggsave(p_bias_gi_var, filename = here("figures/p_bias_gi_var.pdf"), device = cairo_pdf, width = 7, height = 7)
#
p_bias_gi_var_extra <- simres_summdf %>%
  filter(stat == "bias" & grepl("sigma", name) & trt == "gi" & (scenario > 24 | scenario == 9)) |>
  mutate(name_label = fct_rev(factor(name_label, levels = c("sigma[alpha]^2", "sigma[phi]^2", "sigma[epsilon]^2")))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_extra_lbl, labeller = label_parsed, nrow = 2) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Bias (95% C.I.)", color = "")
p_bias_gi_var_extra
ggsave(p_bias_gi_var_extra, filename = here("figures/p_bias_gi_var_extra.pdf"), device = cairo_pdf, width = 6, height = 4)
#
p_rbias_ci_var <- simres_summdf %>%
  filter(stat == "rbias" & grepl("sigma", name) & trt == "ci" & scenario <= 24) |>
  mutate(name_label = fct_rev(factor(name_label, levels = c("sigma[alpha]^2", "sigma[phi]^2", "sigma[epsilon]^2")))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_grid(delta + nu ~ omega, labeller = label_parsed) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Relative Bias (95% C.I.)", color = "")
p_rbias_ci_var
ggsave(p_rbias_ci_var, filename = here("figures/p_rbias_ci_var.pdf"), device = cairo_pdf, width = 7, height = 7)
#
p_rbias_ci_var_extra <- simres_summdf %>%
  filter(stat == "rbias" & grepl("sigma", name) & trt == "ci" & (scenario > 24 | scenario == 9)) |>
  mutate(name_label = fct_rev(factor(name_label, levels = c("sigma[alpha]^2", "sigma[phi]^2", "sigma[epsilon]^2")))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_extra_lbl, labeller = label_parsed, nrow = 2) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Relative Bias (95% C.I.)", color = "")
p_rbias_ci_var_extra
ggsave(p_rbias_ci_var_extra, filename = here("figures/p_rbias_ci_var_extra.pdf"), device = cairo_pdf, width = 6, height = 4)
#
p_rbias_gi_var <- simres_summdf %>%
  filter(stat == "rbias" & grepl("sigma", name) & trt == "gi" & scenario <= 24) |>
  mutate(name_label = fct_rev(factor(name_label, levels = c("sigma[alpha]^2", "sigma[phi]^2", "sigma[epsilon]^2")))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_grid(delta + nu ~ omega, labeller = label_parsed) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Relative Bias (95% C.I.)", color = "")
p_rbias_gi_var
ggsave(p_rbias_gi_var, filename = here("figures/p_rbias_gi_var.pdf"), device = cairo_pdf, width = 7, height = 7)
#
p_rbias_gi_var_extra <- simres_summdf %>%
  filter(stat == "rbias" & grepl("sigma", name) & trt == "gi" & (scenario > 24 | scenario == 9)) |>
  mutate(name_label = fct_rev(factor(name_label, levels = c("sigma[alpha]^2", "sigma[phi]^2", "sigma[epsilon]^2")))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_extra_lbl, labeller = label_parsed, nrow = 2) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Relative Bias (95% C.I.)", color = "")
p_rbias_gi_var_extra
ggsave(p_rbias_gi_var_extra, filename = here("figures/p_rbias_gi_var_extra.pdf"), device = cairo_pdf, width = 6, height = 4)
#

# Table for the same quantities
#
sink(file = glue("tables/tbias_ci_var.tex"))
simres_summdf %>%
  filter(stat == "bias" & grepl("sigma", name) & trt == "ci") |>
  mutate(name_tex = factor(name_tex, levels = c("$\\sigma_{\\alpha}^2$", "$\\sigma_{\\phi}^2$", "$\\sigma_{\\epsilon}^2$"))) |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(Scenario = glue("{omega_tex}, {delta_tex}, {nu_tex}")) |>
  mutate(Scenario = ifelse(scenario > 24, glue("{scenario_extra_tex}"), Scenario)) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    linesep = "",
    caption = "Bias of variance components for the constant treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. \\label{tab:tbias_ci_var}",
    escape = FALSE
  ) |>
  kable_styling(latex_options = c("repeat_header")) |>
  pack_rows("Main scenarios:", 1, 24 * 3) %>%
  pack_rows("Additional scenarios:", 24 * 3 + 1, 28 * 3) |>
  print()
sink()
#
sink(file = glue("tables/tbias_gi_var.tex"))
simres_summdf %>%
  filter(stat == "bias" & grepl("sigma", name) & trt == "gi") |>
  mutate(name_tex = factor(name_tex, levels = c("$\\sigma_{\\alpha}^2$", "$\\sigma_{\\phi}^2$", "$\\sigma_{\\epsilon}^2$"))) |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(Scenario = glue("{omega_tex}, {delta_tex}, {nu_tex}")) |>
  mutate(Scenario = ifelse(scenario > 24, glue("{scenario_extra_tex}"), Scenario)) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    linesep = "",
    caption = "Bias of variance components for the general time on treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. \\label{tab:tbias_gi_var}",
    escape = FALSE
  ) |>
  kable_styling(latex_options = c("repeat_header")) |>
  pack_rows("Main scenarios:", 1, 24 * 3) %>%
  pack_rows("Additional scenarios:", 24 * 3 + 1, 28) |>
  print()
sink()
#
sink(file = glue("tables/trbias_ci_var.tex"))
simres_summdf %>%
  filter(stat == "rbias" & grepl("sigma", name) & trt == "ci") |>
  mutate(name_tex = factor(name_tex, levels = c("$\\sigma_{\\alpha}^2$", "$\\sigma_{\\phi}^2$", "$\\sigma_{\\epsilon}^2$"))) |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(Scenario = glue("{omega_tex}, {delta_tex}, {nu_tex}")) |>
  mutate(Scenario = ifelse(scenario > 24, glue("{scenario_extra_tex}"), Scenario)) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    linesep = "",
    caption = "Relative bias of variance components for the constant treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. \\label{tab:trbias_ci_var}",
    escape = FALSE
  ) |>
  kable_styling(latex_options = c("repeat_header")) |>
  pack_rows("Main scenarios:", 1, 24 * 3) %>%
  pack_rows("Additional scenarios:", 24 * 3 + 1, 28 * 3) |>
  print()
sink()
#
sink(file = glue("tables/trbias_gi_var.tex"))
simres_summdf %>%
  filter(stat == "rbias" & grepl("sigma", name) & trt == "gi") |>
  mutate(name_tex = factor(name_tex, levels = c("$\\sigma_{\\alpha}^2$", "$\\sigma_{\\phi}^2$", "$\\sigma_{\\epsilon}^2$"))) |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(Scenario = glue("{omega_tex}, {delta_tex}, {nu_tex}")) |>
  mutate(Scenario = ifelse(scenario > 24, glue("{scenario_extra_tex}"), Scenario)) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    linesep = "",
    caption = "Relative bias of variance components for the general time on treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. \\label{tab:trbias_gi_var}",
    escape = FALSE
  ) |>
  kable_styling(latex_options = c("repeat_header")) |>
  pack_rows("Main scenarios:", 1, 24 * 3) %>%
  pack_rows("Additional scenarios:", 24 * 3 + 1, 28 * 3) |>
  print()
sink()

# Coverage of variance components
p_cover_ci_var <- simres_summdf %>%
  filter(stat == "cover" & grepl("sigma", name) & trt == "ci" & scenario <= 24) |>
  mutate(name_label = fct_rev(factor(name_label, levels = c("sigma[alpha]^2", "sigma[phi]^2", "sigma[epsilon]^2")))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_grid(delta + nu ~ omega, labeller = label_parsed) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Coverage Probability (95% C.I.)", color = "")
p_cover_ci_var
ggsave(p_cover_ci_var, filename = here("figures/p_cover_ci_var.pdf"), device = cairo_pdf, width = 7, height = 7)
#
p_cover_ci_var_extra <- simres_summdf %>%
  filter(stat == "cover" & grepl("sigma", name) & trt == "ci" & (scenario > 24 | scenario == 9)) |>
  mutate(name_label = fct_rev(factor(name_label, levels = c("sigma[alpha]^2", "sigma[phi]^2", "sigma[epsilon]^2")))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_extra_lbl, labeller = label_parsed, nrow = 2) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Coverage Probability (95% C.I.)", color = "")
p_cover_ci_var_extra
ggsave(p_cover_ci_var_extra, filename = here("figures/p_cover_ci_var_extra.pdf"), device = cairo_pdf, width = 6, height = 4)
#
p_cover_gi_var <- simres_summdf %>%
  filter(stat == "cover" & grepl("sigma", name) & trt == "gi" & scenario <= 24) |>
  mutate(name_label = fct_rev(factor(name_label, levels = c("sigma[alpha]^2", "sigma[phi]^2", "sigma[epsilon]^2")))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_grid(delta + nu ~ omega, labeller = label_parsed) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Coverage Probability (95% C.I.)", color = "")
p_cover_gi_var
ggsave(p_cover_gi_var, filename = here("figures/p_cover_gi_var.pdf"), device = cairo_pdf, width = 7, height = 7)
#
p_cover_gi_var_extra <- simres_summdf %>%
  filter(stat == "cover" & grepl("sigma", name) & trt == "gi" & (scenario > 24 | scenario == 9)) |>
  mutate(name_label = fct_rev(factor(name_label, levels = c("sigma[alpha]^2", "sigma[phi]^2", "sigma[epsilon]^2")))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_extra_lbl, labeller = label_parsed, nrow = 2) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Coverage Probability (95% C.I.)", color = "")
p_cover_gi_var_extra
ggsave(p_cover_gi_var_extra, filename = here("figures/p_cover_gi_var_extra.pdf"), device = cairo_pdf, width = 6, height = 4)

# Table for the same quantities
#
sink(file = glue("tables/tcover_ci_var.tex"))
simres_summdf %>%
  filter(stat == "cover" & grepl("sigma", name) & trt == "ci") |>
  mutate(name_tex = factor(name_tex, levels = c("$\\sigma_{\\alpha}^2$", "$\\sigma_{\\phi}^2$", "$\\sigma_{\\epsilon}^2$"))) |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(Scenario = glue("{omega_tex}, {delta_tex}, {nu_tex}")) |>
  mutate(Scenario = ifelse(scenario > 24, glue("{scenario_extra_tex}"), Scenario)) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    linesep = "",
    caption = "Coverage probability of variance components for the constant treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. \\label{tab:tcover_ci_var}",
    escape = FALSE
  ) |>
  kable_styling(latex_options = c("repeat_header")) |>
  pack_rows("Main scenarios:", 1, 24 * 3) %>%
  pack_rows("Additional scenarios:", 24 * 3 + 1, 28 * 3) |>
  print()
sink()
#
sink(file = glue("tables/tcover_gi_var.tex"))
simres_summdf %>%
  filter(stat == "cover" & grepl("sigma", name) & trt == "gi") |>
  mutate(name_tex = factor(name_tex, levels = c("$\\sigma_{\\alpha}^2$", "$\\sigma_{\\phi}^2$", "$\\sigma_{\\epsilon}^2$"))) |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(Scenario = glue("{omega_tex}, {delta_tex}, {nu_tex}")) |>
  mutate(Scenario = ifelse(scenario > 24, glue("{scenario_extra_tex}"), Scenario)) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    linesep = "",
    caption = "Coverage probability of variance components for the general time on treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. \\label{tab:tcover_gi_var}",
    escape = FALSE
  ) |>
  kable_styling(latex_options = c("repeat_header")) |>
  pack_rows("Main scenarios:", 1, 24 * 3) %>%
  pack_rows("Additional scenarios:", 24 * 3 + 1, 28 * 3) |>
  print()
sink()

# Bias of ICCs
p_bias_ci_icc <- simres_summdf %>%
  filter(stat == "bias" & grepl("icc", name) & trt == "ci" & scenario <= 24) |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::comma, breaks = c(-0.3, -0.2, -0.1, 0.0)) +
  scale_color_rda +
  theme_rda +
  facet_grid(nu + delta ~ omega, labeller = label_parsed) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Bias (95% C.I.)", color = "")
p_bias_ci_icc
ggsave(p_bias_ci_icc, filename = here("figures/p_bias_ci_icc.pdf"), device = cairo_pdf, width = 7, height = 6)
#
p_bias_ci_icc_extra <- simres_summdf %>%
  filter(stat == "bias" & grepl("icc", name) & trt == "ci" & (scenario > 24 | scenario == 9)) |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::comma, breaks = c(-0.3, -0.2, -0.1, 0.0)) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_extra_lbl, labeller = label_parsed, nrow = 2) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Bias (95% C.I.)", color = "")
p_bias_ci_icc_extra
ggsave(p_bias_ci_icc_extra, filename = here("figures/p_bias_ci_icc_extra.pdf"), device = cairo_pdf, width = 6, height = 4)
#
p_bias_gi_icc <- simres_summdf %>%
  filter(stat == "bias" & grepl("icc", name) & trt == "gi" & scenario <= 24) |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::comma, breaks = c(-0.3, -0.2, -0.1, 0.0)) +
  scale_color_rda +
  theme_rda +
  facet_grid(nu + delta ~ omega, labeller = label_parsed) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Bias (95% C.I.)", color = "")
p_bias_gi_icc
ggsave(p_bias_gi_icc, filename = here("figures/p_bias_gi_icc.pdf"), device = cairo_pdf, width = 7, height = 6)
#
p_bias_gi_icc_extra <- simres_summdf %>%
  filter(stat == "bias" & grepl("icc", name) & trt == "gi" & (scenario > 24 | scenario == 9)) |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::comma, breaks = c(-0.3, -0.2, -0.1, 0.0)) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_extra_lbl, labeller = label_parsed, nrow = 2) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Bias (95% C.I.)", color = "")
p_bias_gi_icc_extra
ggsave(p_bias_gi_icc_extra, filename = here("figures/p_bias_gi_icc_extra.pdf"), device = cairo_pdf, width = 6, height = 4)
#
p_rbias_ci_icc <- simres_summdf %>%
  filter(stat == "rbias" & grepl("icc", name) & trt == "ci" & scenario <= 24) |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_grid(nu + delta ~ omega, labeller = label_parsed) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Relative Bias (95% C.I.)", color = "")
p_rbias_ci_icc
ggsave(p_rbias_ci_icc, filename = here("figures/p_rbias_ci_icc.pdf"), device = cairo_pdf, width = 7, height = 6)
#
p_rbias_ci_icc_extra <- simres_summdf %>%
  filter(stat == "rbias" & grepl("icc", name) & trt == "ci" & (scenario > 24 | scenario == 9)) |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_extra_lbl, labeller = label_parsed, nrow = 2) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Relative Bias (95% C.I.)", color = "")
p_rbias_ci_icc_extra
ggsave(p_rbias_ci_icc_extra, filename = here("figures/p_rbias_ci_icc_extra.pdf"), device = cairo_pdf, width = 6, height = 4)
#
p_rbias_gi_icc <- simres_summdf %>%
  filter(stat == "rbias" & grepl("icc", name) & trt == "gi" & scenario <= 24) |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_grid(nu + delta ~ omega, labeller = label_parsed) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Relative Bias (95% C.I.)", color = "")
p_rbias_gi_icc
ggsave(p_rbias_gi_icc, filename = here("figures/p_rbias_gi_icc.pdf"), device = cairo_pdf, width = 7, height = 6)
#
p_rbias_gi_icc_extra <- simres_summdf %>%
  filter(stat == "rbias" & grepl("icc", name) & trt == "gi" & (scenario > 24 | scenario == 9)) |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_extra_lbl, labeller = label_parsed, nrow = 2) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Relative Bias (95% C.I.)", color = "")
p_rbias_gi_icc_extra
ggsave(p_rbias_gi_icc_extra, filename = here("figures/p_rbias_gi_icc_extra.pdf"), device = cairo_pdf, width = 6, height = 4)
#

# Table for the same quantities
#
sink(file = glue("tables/tbias_ci_icc.tex"))
simres_summdf %>%
  filter(stat == "bias" & grepl("icc", name) & trt == "ci") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(Scenario = glue("{omega_tex}, {delta_tex}, {nu_tex}")) |>
  mutate(Scenario = ifelse(scenario > 24, glue("{scenario_extra_tex}"), Scenario)) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    linesep = "",
    caption = "Bias of ICCs for the constant treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. \\label{tab:tbias_ci_icc}",
    escape = FALSE
  ) |>
  kable_styling(latex_options = c("repeat_header")) |>
  pack_rows("Main scenarios:", 1, 24 * 2) %>%
  pack_rows("Additional scenarios:", 24 * 2 + 1, 28 * 2) |>
  print()
sink()
#
sink(file = glue("tables/tbias_gi_icc.tex"))
simres_summdf %>%
  filter(stat == "bias" & grepl("icc", name) & trt == "gi") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(Scenario = glue("{omega_tex}, {delta_tex}, {nu_tex}")) |>
  mutate(Scenario = ifelse(scenario > 24, glue("{scenario_extra_tex}"), Scenario)) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    linesep = "",
    caption = "Bias of ICCs for the general time on treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. \\label{tab:tbias_gi_icc}",
    escape = FALSE
  ) |>
  kable_styling(latex_options = c("repeat_header")) |>
  pack_rows("Main scenarios:", 1, 24 * 2) %>%
  pack_rows("Additional scenarios:", 24 * 2 + 1, 28 * 2) |>
  print()
sink()
#
sink(file = glue("tables/trbias_ci_icc.tex"))
simres_summdf %>%
  filter(stat == "rbias" & grepl("icc", name) & trt == "ci") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(Scenario = glue("{omega_tex}, {delta_tex}, {nu_tex}")) |>
  mutate(Scenario = ifelse(scenario > 24, glue("{scenario_extra_tex}"), Scenario)) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    linesep = "",
    caption = "Relative bias of ICCs for the constant treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. \\label{tab:trbias_ci_icc}",
    escape = FALSE
  ) |>
  kable_styling(latex_options = c("repeat_header")) |>
  pack_rows("Main scenarios:", 1, 24 * 2) %>%
  pack_rows("Additional scenarios:", 24 * 2 + 1, 28 * 2) |>
  print()
sink()
#
sink(file = glue("tables/trbias_gi_icc.tex"))
simres_summdf %>%
  filter(stat == "rbias" & grepl("icc", name) & trt == "gi") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(Scenario = glue("{omega_tex}, {delta_tex}, {nu_tex}")) |>
  mutate(Scenario = ifelse(scenario > 24, glue("{scenario_extra_tex}"), Scenario)) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    linesep = "",
    caption = "Relative bias of ICCs for the general time on treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. \\label{tab:trbias_gi_icc}",
    escape = FALSE
  ) |>
  kable_styling(latex_options = c("repeat_header")) |>
  pack_rows("Main scenarios:", 1, 24 * 2) %>%
  pack_rows("Additional scenarios:", 24 * 2 + 1, 28 * 2) |>
  print()
sink()

# Coverage of ICCs
p_cover_ci_icc <- simres_summdf %>%
  filter(stat == "cover" & grepl("icc", name) & trt == "ci" & scenario <= 24) |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_grid(nu + delta ~ omega, labeller = label_parsed) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Coverage Probability (95% C.I.)", color = "")
p_cover_ci_icc
ggsave(p_cover_ci_icc, filename = here("figures/p_cover_ci_icc.pdf"), device = cairo_pdf, width = 7, height = 6)
#
p_cover_ci_icc_extra <- simres_summdf %>%
  filter(stat == "cover" & grepl("icc", name) & trt == "ci" & (scenario > 24 | scenario == 9)) |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_extra_lbl, labeller = label_parsed, nrow = 2) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Coverage Probability (95% C.I.)", color = "")
p_cover_ci_icc_extra
ggsave(p_cover_ci_icc_extra, filename = here("figures/p_cover_ci_icc_extra.pdf"), device = cairo_pdf, width = 7, height = 6)
#
p_cover_gi_icc <- simres_summdf %>%
  filter(stat == "cover" & grepl("icc", name) & trt == "gi" & scenario <= 24) |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_grid(nu + delta ~ omega, labeller = label_parsed) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Coverage Probability (95% C.I.)", color = "")
p_cover_gi_icc
ggsave(p_cover_gi_icc, filename = here("figures/p_cover_gi_icc.pdf"), device = cairo_pdf, width = 7, height = 7)
#
p_cover_gi_icc_extra <- simres_summdf %>%
  filter(stat == "cover" & grepl("icc", name) & trt == "gi" & (scenario > 24 | scenario == 9)) |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_extra_lbl, labeller = label_parsed, nrow = 2) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Coverage Probability (95% C.I.)", color = "")
p_cover_gi_icc_extra
ggsave(p_cover_gi_icc_extra, filename = here("figures/p_cover_gi_icc_extra.pdf"), device = cairo_pdf, width = 7, height = 7)

# Table for the same quantities
#
sink(file = glue("tables/tcover_ci_icc.tex"))
simres_summdf %>%
  filter(stat == "cover" & grepl("icc", name) & trt == "ci") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(Scenario = glue("{omega_tex}, {delta_tex}, {nu_tex}")) |>
  mutate(Scenario = ifelse(scenario > 24, glue("{scenario_extra_tex}"), Scenario)) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    linesep = "",
    caption = "Coverage probability of ICCs for the constant treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. \\label{tab:tcover_ci_icc}",
    escape = FALSE
  ) |>
  kable_styling(latex_options = c("repeat_header")) |>
  pack_rows("Main scenarios:", 1, 24 * 2) %>%
  pack_rows("Additional scenarios:", 24 * 2 + 1, 28 * 2) |>
  print()
sink()
#
sink(file = glue("tables/tcover_gi_icc.tex"))
simres_summdf %>%
  filter(stat == "cover" & grepl("icc", name) & trt == "gi") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(Scenario = glue("{omega_tex}, {delta_tex}, {nu_tex}")) |>
  mutate(Scenario = ifelse(scenario > 24, glue("{scenario_extra_tex}"), Scenario)) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    linesep = "",
    caption = "Coverage probability of ICCs for the general time on treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. \\label{tab:tcover_gi_icc}",
    escape = FALSE
  ) |>
  kable_styling(latex_options = c("repeat_header")) |>
  pack_rows("Main scenarios:", 1, 24 * 2) %>%
  pack_rows("Additional scenarios:", 24 * 2 + 1, 28 * 2) |>
  print()
sink()

# Combined bias for treatment effects (ci, gi models)
p_bias_cigi_delta_extra <- simres_summdf %>%
  filter(stat == "bias" & grepl("^delta", name) & (scenario > 24 | scenario == 9)) |>
  mutate(trt = factor(trt, levels = c("ci", "gi"), labels = c("Constant ~ Intervention", "General ~ Time ~ on ~ Treatment"))) |>
  ggplot(aes(x = fct_rev(name_label), y = est, color = model, shape = name_label)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 2 / 3)) +
  geom_point(position = position_dodge(width = 2 / 3)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_rda +
  theme_rda +
  facet_grid(scenario_extra_lbl ~ trt, , labeller = label_parsed, scales = "free_x") +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Bias (95% C.I.)", color = "", shape = "") +
  guides(shape = "none")
p_bias_cigi_delta_extra
ggsave(p_bias_cigi_delta_extra, filename = "figures/p_bias_cigi_delta_extra.pdf", device = cairo_pdf, width = 7, height = 7.7)

# Largest MCSE for bias, across scenarios, by delta, beta_j, ICC:
df_mcse <- simres_summdf |>
  filter(stat %in% c("bias", "nsim")) |>
  select(trt, scenario, model, stat, est, mcse, name_label) |>
  pivot_wider(names_from = "stat", values_from = c("est", "mcse")) |>
  mutate(trt = factor(trt, levels = c("ci", "gi"), labels = c("Constant Treatment Effect", "General Time on Treatment")))
p_mcse_bias <- df_mcse |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = mcse_bias)) +
  geom_boxplot() +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_log10(labels = scales::comma) +
  facet_wrap(~trt) +
  theme_rda +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Monte Carlo Standard Error for Bias")
p_mcse_bias
ggsave(p_mcse_bias, filename = "figures/p_mcse_bias.pdf", device = cairo_pdf, width = 6, height = 6)
#
p_mcse_bias_conv <- df_mcse |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = mcse_bias, color = est_nsim)) +
  geom_jitter(height = 0, width = 0.2) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_log10(labels = scales::comma) +
  scale_color_viridis_c(option = "A", direction = -1, end = 4 / 5) +
  facet_wrap(~trt) +
  theme_rda +
  theme(legend.position = "top") +
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(4, "in"))) +
  coord_flip() +
  labs(x = "", y = "Monte Carlo Standard Error for Bias", color = "Converged Repetitions")
p_mcse_bias_conv
ggsave(p_mcse_bias_conv, filename = "figures/p_mcse_bias_conv.pdf", device = cairo_pdf, width = 6, height = 6)

# Checking MCSE of rbias for delta, delta[j] parameters no more than 2%
df_mcse_rbias <- simres_summdf |>
  filter(stat %in% c("rbias", "nsim") & grepl("delta", name_label)) |>
  select(trt, scenario, model, stat, est, mcse, name_label) |>
  pivot_wider(names_from = "stat", values_from = c("est", "mcse")) |>
  mutate(trt = factor(trt, levels = c("ci", "gi"), labels = c("Constant Treatment Effect", "General Time on Treatment")))
p_mcse_rbias <- df_mcse_rbias |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = mcse_rbias)) +
  geom_boxplot() +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_log10(labels = scales::comma) +
  facet_wrap(~trt) +
  theme_rda +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Monte Carlo Standard Error for Relative Bias")
p_mcse_rbias
ggsave(p_mcse_rbias, filename = "figures/p_mcse_rbias.pdf", device = cairo_pdf, width = 6, height = 3)
#
p_mcse_rbias_conv <- df_mcse_rbias |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = mcse_rbias, color = est_nsim)) +
  geom_jitter(height = 0, width = 0.2) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_log10(labels = scales::comma) +
  scale_color_viridis_c(option = "A", direction = -1, end = 4 / 5) +
  facet_wrap(~trt) +
  theme_rda +
  theme(legend.position = "top") +
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(4, "in"))) +
  coord_flip() +
  labs(x = "", y = "Monte Carlo Standard Error for Relative Bias", color = "Converged Repetitions")
p_mcse_rbias_conv
ggsave(p_mcse_rbias_conv, filename = "figures/p_mcse_rbias_conv.pdf", device = cairo_pdf, width = 6, height = 3)

# Dropout rates by simulation scenario
lf <- list.files(path = here("simulation/data/simdata/"), pattern = "01-|04-")
# lf <- sample(x = lf, size = 1000)
drdt <- map_dfr(.x = lf, .f = function(xx) {
  filedf <- data.frame(tmp = xx) |>
    mutate(tmp2 = str_sub(tmp, end = -5)) |>
    separate(col = "tmp2", into = c("drop1", "analysis", "trt", "scenario", "index")) |>
    select(-drop1, -analysis, -index)
  this <- read_dta(file = here(glue("simulation/data/simdata/{xx}"))) |>
    zap_formats() |>
    zap_label() |>
    zap_labels() |>
    zap_missing() |>
    zap_widths() |>
    mutate(tmp = xx) |>
    select(-cumx, -y, -yobs, -eventtime, -actual_eventtime) |>
    left_join(filedf, by = "tmp") |>
    mutate(scenario = as.numeric(scenario)) |>
    select(-tmp) |>
    group_by(trt, scenario, j, x) |>
    summarise(mean_d = mean(d, na.rm = TRUE), .groups = "drop") |>
    ungroup()
  return(this)
}, .progress = TRUE)
drdt_for_plot <- drdt |>
  group_by(trt, scenario, j, x) |>
  summarise(mean_d = mean(mean_d, na.rm = TRUE)) |>
  ungroup() |>
  left_join(rename(scenarios_wide, i_design = i, j_design = j), by = "scenario") |>
  mutate(facet = glue("list(list(omega[1], omega[2]) == log({exp(omega1)}), delta == {delta}, nu == {nu})")) |>
  mutate(facet = factor(facet)) |>
  mutate(facet = fct_reorder(.f = facet, .x = scenario)) |>
  mutate(x = factor(x, levels = c(0, 1), labels = c("No Treatment", "Treatment"))) |>
  mutate(mean_d_lab = percent(mean_d, 2)) |>
  left_join(scenarios_labels)
# CI:
p_dr_ci <- ggplot(filter(drdt_for_plot, trt == "ci"), aes(x = factor(j), y = mean_d, fill = x)) +
  geom_col(position = position_dodge(width = 0.75), width = 2 / 3) +
  geom_text(aes(label = mean_d_lab), position = position_dodge(width = 0.75), angle = 90, hjust = -0.1, family = "Atkinson Hyperlegible", size = 3) +
  scale_y_continuous(labels = function(x) percent(x, 0)) +
  facet_wrap(~scenario_lbl, labeller = label_parsed, ncol = 3) +
  theme_rda +
  scale_fill_grey() +
  theme(legend.position = "top") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Period", y = "Dropout Proportion", fill = "")
p_dr_ci
ggsave(p_dr_ci, filename = "figures/p_dr_ci.pdf", device = cairo_pdf, width = 17, height = 14)
# GI:
p_dr_gi <- ggplot(filter(drdt_for_plot, trt == "gi"), aes(x = factor(j), y = mean_d, fill = x)) +
  geom_col(position = position_dodge(width = 0.75), width = 2 / 3) +
  geom_text(aes(label = mean_d_lab), position = position_dodge(width = 0.75), angle = 90, hjust = -0.1, family = "Atkinson Hyperlegible", size = 3) +
  scale_y_continuous(labels = function(x) percent(x, 0)) +
  facet_wrap(~scenario_lbl, labeller = label_parsed, ncol = 3) +
  theme_rda +
  scale_fill_grey() +
  theme(legend.position = "top") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Period", y = "Dropout Proportion", fill = "")
p_dr_gi
ggsave(p_dr_gi, filename = "figures/p_dr_gi.pdf", device = cairo_pdf, width = 17, height = 14)
