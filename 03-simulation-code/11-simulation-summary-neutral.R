# With this script, we summarise neutral scenario

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
simres_ci <- read_dta(file = here("simulation/data/05-full-results-ci.dta")) |>
  mutate(trt = "ci") |>
  filter(scenario == 9)
simres_gi <- read_dta(file = "simulation/data/05-full-results-gi.dta") |>
  mutate(trt = "gi") |>
  filter(scenario == 9)
simres_ci_neutral <- read_dta(file = here("simulation/data/07-full-results-ci.dta")) |>
  mutate(trt = "ci", scenario = 99)
simres_gi_neutral <- read_dta(file = here("simulation/data/07-full-results-gi.dta")) |>
  mutate(trt = "gi", scenario = 99)

# Combine
simres <- bind_rows(simres_ci, simres_ci_neutral, simres_gi, simres_gi_neutral) |>
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
scenarios_wide <- bind_rows(
  readRDS(file = here("simulation/data/01-scenarios.RDS")) |> filter(scenario == 9),
  readRDS(file = here("simulation/data/06-scenarios-neutral.RDS")) |> mutate(scenario = 99)
)
scenarios_labels <- scenarios_wide |>
  select(scenario) |>
  mutate(scenario_lbl = factor(
    scenario,
    levels = c(9, 99),
    labels = c("Reference", "Neutral")
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
  left_join(scenarios_labels) |>
  mutate(trt = factor(trt, levels = c("ci", "gi"), labels = c("Constant Treatment Effect", "General Time on Treatment")))
# Plot:
p_conv_neutral <- ggplot(simres_convergence, aes(x = fct_rev(scenario_lbl), y = pnc, color = model, fill = model)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  geom_col(position = position_dodge(width = 4 / 5), width = 0.5) +
  geom_text(aes(label = lbl, family = "Atkinson Hyperlegible"), position = position_dodge(width = 4 / 5), color = "grey25", hjust = -0.1, size = 2.5) +
  scale_color_rda +
  scale_fill_rda +
  theme_rda +
  theme(legend.position = "top") +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~trt, ncol = 1) +
  coord_flip(ylim = c(0, 1)) +
  labs(x = "", y = "Proportion of Non-Converged Iterations", color = "", fill = "")
p_conv_neutral
ggsave(p_conv_neutral, filename = here("figures/p_conv_neutral.pdf"), device = cairo_pdf, width = 5, height = 4)

# Cleanup
rm(p_conv_neutral, simres_convergence)
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
saveRDS(object = simres_summdf, file = here("simulation/data/11-simres_summdf-neutral.RDS"))

# Plot bias of longitudinal treatment effect(s)
p_bias_ci_delta <- simres_summdf %>%
  filter(stat == "bias" & name == "delta" & trt == "ci") |>
  ggplot(aes(x = fct_rev(scenario_lbl), y = est, color = model)) +
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
p_bias_ci_delta
ggsave(p_bias_ci_delta, filename = here("figures/p_bias_ci_delta_neutral.pdf"), device = cairo_pdf, width = 5, height = 3)
#
p_bias_gi_delta <- simres_summdf %>%
  filter(stat == "bias" & grepl("delta", name) & trt == "gi") |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model, shape = name_label)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 2 / 3)) +
  geom_point(position = position_dodge(width = 2 / 3)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_lbl, ncol = 1) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Bias (95% C.I.)", color = "") +
  guides(shape = "none")
p_bias_gi_delta
ggsave(p_bias_gi_delta, filename = here("figures/p_bias_gi_delta_neutral.pdf"), device = cairo_pdf, width = 5, height = 4)

# Plot relative bias of longitudinal treatment effect(s)
p_rbias_ci_delta <- simres_summdf %>%
  filter(stat == "rbias" & name == "delta" & trt == "ci") |>
  ggplot(aes(x = fct_rev(scenario_lbl), y = est, color = model)) +
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
p_rbias_ci_delta
ggsave(p_rbias_ci_delta, filename = here("figures/p_rbias_ci_delta_neutral.pdf"), device = cairo_pdf, width = 5, height = 3)
#
p_rbias_gi_delta <- simres_summdf %>%
  filter(stat == "rbias" & grepl("delta", name) & trt == "gi") |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model, shape = name_label)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 2 / 3)) +
  geom_point(position = position_dodge(width = 2 / 3)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_lbl, ncol = 1) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Relative Bias (95% C.I.)", color = "") +
  guides(shape = "none")
p_rbias_gi_delta
ggsave(p_rbias_gi_delta, filename = here("figures/p_rbias_gi_delta_neutral.pdf"), device = cairo_pdf, width = 5, height = 3)

# Table for the same quantities
#
sink(file = glue("tables/tbias_ci_delta_neutral.tex"))
simres_summdf %>%
  filter(stat == "bias" & name == "delta" & trt == "ci") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(Scenario = scenario_lbl) |>
  select(Scenario, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    linesep = "",
    caption = "Bias of treatment effect on the longitudinal outcome for the constant treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. Comparison of scenarios based on a joint model and on a neutral data-generating mechanism. \\label{tab:tbias_ci_delta_neutral}",
    escape = FALSE
  ) |>
  print()
sink()
#
sink(file = glue("tables/trbias_ci_delta_neutral.tex"))
simres_summdf %>%
  filter(stat == "rbias" & name == "delta" & trt == "ci") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(y = ifelse(!is.finite(est), "---", y)) |>
  mutate(Scenario = scenario_lbl) |>
  select(Scenario, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    linesep = "",
    caption = "Relative bias of treatment effect on the longitudinal outcome for the constant treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. Comparison of scenarios based on a joint model and on a neutral data-generating mechanism. \\label{tab:trbias_ci_delta_neutral}",
    escape = FALSE
  ) |>
  print()
sink()
#
sink(file = glue("tables/tbias_gi_delta_neutral.tex"))
simres_summdf %>%
  filter(stat == "bias" & grepl("delta", name) & trt == "gi") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(Scenario = scenario_lbl) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    linesep = "",
    caption = "Bias of treatment effect on the longitudinal outcome for the general time on treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. Comparison of scenarios based on a joint model and on a neutral data-generating mechanism. \\label{tab:tbias_gi_delta_neutral}",
    escape = FALSE
  ) |>
  print()
sink()
#
sink(file = glue("tables/trbias_gi_delta_neutral.tex"))
simres_summdf %>%
  filter(stat == "rbias" & grepl("delta", name) & trt == "gi") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(y = ifelse(!is.finite(est), "---", y)) |>
  mutate(Scenario = scenario_lbl) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    linesep = "",
    caption = "Relative bias of treatment effect on the longitudinal outcome for the general time on treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. Comparison of scenarios based on a joint model and on a neutral data-generating mechanism. \\label{tab:trbias_gi_delta_neutral}",
    escape = FALSE
  ) |>
  print()
sink()

# Plot bias for the period effects
p_bias_ci_beta <- simres_summdf %>%
  filter(stat == "bias" & grepl("beta", name) & trt == "ci") |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  geom_point() +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_lbl, ncol = 1) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Bias (95% C.I.)", color = "")
p_bias_ci_beta
ggsave(p_bias_ci_beta, filename = here("figures/p_bias_ci_beta_neutral.pdf"), device = cairo_pdf, width = 5, height = 4)
#
p_bias_gi_beta <- simres_summdf %>%
  filter(stat == "bias" & grepl("beta", name) & trt == "gi") |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  geom_point() +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_lbl, ncol = 1) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Bias (95% C.I.)", color = "")
p_bias_gi_beta
ggsave(p_bias_gi_beta, filename = here("figures/p_bias_gi_beta_neutral.pdf"), device = cairo_pdf, width = 5, height = 4)

# Plot relative bias for the period effects
p_rbias_ci_beta <- simres_summdf %>%
  filter(stat == "rbias" & grepl("beta", name) & trt == "ci") |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  geom_point() +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_lbl, ncol = 1) +
  theme(legend.position = "bottom") +
  coord_flip(ylim = c(-0.3, 0.3)) +
  labs(x = "", y = "Relative Bias (95% C.I.)", color = "")
p_rbias_ci_beta
ggsave(p_rbias_ci_beta, filename = here("figures/p_rbias_ci_beta_neutral.pdf"), device = cairo_pdf, width = 5, height = 4)
#
p_rbias_gi_beta <- simres_summdf %>%
  filter(stat == "rbias" & grepl("beta", name) & trt == "gi") |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  geom_point() +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_lbl, ncol = 1) +
  theme(legend.position = "bottom") +
  coord_flip(ylim = c(-0.3, 0.3)) +
  labs(x = "", y = "Relative Bias (95% C.I.)", color = "")
p_rbias_gi_beta
ggsave(p_rbias_gi_beta, filename = here("figures/p_rbias_gi_beta_neutral.pdf"), device = cairo_pdf, width = 5, height = 4)

# Table for the same quantities
#
sink(file = glue("tables/tbias_ci_beta_neutral.tex"))
simres_summdf %>%
  filter(stat == "bias" & grepl("beta", name) & trt == "ci") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(Scenario = scenario_lbl) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    linesep = "",
    caption = "Bias of period effects on the longitudinal outcome for the constant treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. Comparison of scenarios based on a joint model and on a neutral data-generating mechanism. \\label{tab:tbias_ci_beta_neutral}",
    escape = FALSE
  ) |>
  print()
sink()
#
sink(file = glue("tables/tbias_gi_beta_neutral.tex"))
simres_summdf %>%
  filter(stat == "bias" & grepl("beta", name) & trt == "gi") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(Scenario = scenario_lbl) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    linesep = "",
    caption = "Bias of period effects on the longitudinal outcome for the general time on treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. Comparison of scenarios based on a joint model and on a neutral data-generating mechanism. \\label{tab:tbias_gi_beta_neutral}",
    escape = FALSE
  ) |>
  print()
sink()
#
sink(file = glue("tables/trbias_ci_beta_neutral.tex"))
simres_summdf %>%
  filter(stat == "rbias" & grepl("beta", name) & trt == "ci") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(Scenario = scenario_lbl) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    linesep = "",
    caption = "Relative bias of period effects on the longitudinal outcome for the constant treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. Comparison of scenarios based on a joint model and on a neutral data-generating mechanism. \\label{tab:trbias_ci_beta_neutral}",
    escape = FALSE
  ) |>
  print()
sink()
#
sink(file = glue("tables/trbias_gi_beta_neutral.tex"))
simres_summdf %>%
  filter(stat == "rbias" & grepl("beta", name) & trt == "gi") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(Scenario = scenario_lbl) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    linesep = "",
    caption = "Relative bias of period effects on the longitudinal outcome for the general time on treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. Comparison of scenarios based on a joint model and on a neutral data-generating mechanism. \\label{tab:trbias_gi_beta_neutral}",
    escape = FALSE
  ) |>
  print()
sink()

# Plot coverage probability of longitudinal treatment effect(s)
p_cover_ci_delta <- simres_summdf %>%
  filter(stat == "cover" & name == "delta" & trt == "ci") |>
  ggplot(aes(x = fct_rev(scenario_lbl), y = est, color = model)) +
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
p_cover_ci_delta
ggsave(p_cover_ci_delta, filename = here("figures/p_cover_ci_delta_neutral.pdf"), device = cairo_pdf, width = 5, height = 3)
#
p_cover_gi_delta <- simres_summdf %>%
  filter(stat == "cover" & grepl("delta", name) & trt == "gi") |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model, shape = name_label)) +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 2 / 3)) +
  geom_point(position = position_dodge(width = 2 / 3)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_lbl, ncol = 1) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Coverage Probability (95% C.I.)", color = "") +
  guides(shape = "none")
p_cover_gi_delta
ggsave(p_cover_gi_delta, filename = here("figures/p_cover_gi_delta_neutral.pdf"), device = cairo_pdf, width = 5, height = 4)

# Table for the same quantities
#
sink(file = glue("tables/tcover_ci_delta_neutral.tex"))
simres_summdf %>%
  filter(stat == "cover" & name == "delta" & trt == "ci") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(Scenario = scenario_lbl) |>
  select(Scenario, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    linesep = "",
    caption = "Coverage probability of treatment effect on the longitudinal outcome for the constant treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Comparison of scenarios based on a joint model and on a neutral data-generating mechanism. \\label{tab:tcover_ci_delta_neutral}",
    escape = FALSE
  ) |>
  print()
sink()
#
sink(file = glue("tables/tcover_gi_delta_neutral.tex"))
simres_summdf %>%
  filter(stat == "cover" & grepl("delta", name) & trt == "gi") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(Scenario = scenario_lbl) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    linesep = "",
    caption = "Coverage probability of treatment effect on the longitudinal outcome for the general time on treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Comparison of scenarios based on a joint model and on a neutral data-generating mechanism. \\label{tab:tcover_gi_delta_neutral}",
    escape = FALSE
  ) |>
  print()
sink()
#

# Plot coverage probability of period effects(s)
p_cover_ci_beta <- simres_summdf %>%
  filter(stat == "cover" & grepl("beta", name) & trt == "ci") |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_lbl, ncol = 1) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Coverage Probability (95% C.I.)", color = "")
p_cover_ci_beta
ggsave(p_cover_ci_beta, filename = here("figures/p_cover_ci_beta_neutral.pdf"), device = cairo_pdf, width = 5, height = 4)
#
p_cover_gi_beta <- simres_summdf %>%
  filter(stat == "cover" & grepl("beta", name) & trt == "gi") |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_lbl, ncol = 1) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Coverage Probability (95% C.I.)", color = "")
p_cover_gi_beta
ggsave(p_cover_gi_beta, filename = here("figures/p_cover_gi_beta_neutral.pdf"), device = cairo_pdf, width = 5, height = 4)

# Table for the same quantities
#
sink(file = glue("tables/tcover_ci_beta_neutral.tex"))
simres_summdf %>%
  filter(stat == "cover" & grepl("beta", name) & trt == "ci") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(Scenario = scenario_lbl) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    linesep = "",
    caption = "Coverage probability of period effects on the longitudinal outcome for the constant treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Comparison of scenarios based on a joint model and on a neutral data-generating mechanism. \\label{tab:tcover_ci_beta_neutral}",
    escape = FALSE
  ) |>
  print()
sink()
#
sink(file = glue("tables/tcover_gi_beta_neutral.tex"))
simres_summdf %>%
  filter(stat == "cover" & grepl("beta", name) & trt == "gi") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(Scenario = scenario_lbl) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    linesep = "",
    caption = "Coverage probability of period effects on the longitudinal outcome for the general time on treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Comparison of scenarios based on a joint model and on a neutral data-generating mechanism. \\label{tab:tcover_gi_beta_neutral}",
    escape = FALSE
  ) |>
  print()
sink()

# Bias of variance components
p_bias_ci_var <- simres_summdf %>%
  filter(stat == "bias" & grepl("sigma", name) & trt == "ci") |>
  mutate(name_label = fct_rev(factor(name_label, levels = c("sigma[alpha]^2", "sigma[phi]^2", "sigma[epsilon]^2")))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_lbl, ncol = 1) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Bias (95% C.I.)", color = "")
p_bias_ci_var
ggsave(p_bias_ci_var, filename = here("figures/p_bias_ci_var_neutral.pdf"), device = cairo_pdf, width = 5, height = 4)
#
p_bias_gi_var <- simres_summdf %>%
  filter(stat == "bias" & grepl("sigma", name) & trt == "gi") |>
  mutate(name_label = fct_rev(factor(name_label, levels = c("sigma[alpha]^2", "sigma[phi]^2", "sigma[epsilon]^2")))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_lbl, ncol = 1) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Bias (95% C.I.)", color = "")
p_bias_gi_var
ggsave(p_bias_gi_var, filename = here("figures/p_bias_gi_var_neutral.pdf"), device = cairo_pdf, width = 5, height = 4)
#
p_rbias_ci_var <- simres_summdf %>%
  filter(stat == "rbias" & grepl("sigma", name) & trt == "ci") |>
  mutate(name_label = fct_rev(factor(name_label, levels = c("sigma[alpha]^2", "sigma[phi]^2", "sigma[epsilon]^2")))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_lbl, ncol = 1) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Relative Bias (95% C.I.)", color = "")
p_rbias_ci_var
ggsave(p_rbias_ci_var, filename = here("figures/p_rbias_ci_var_neutral.pdf"), device = cairo_pdf, width = 5, height = 4)
#
p_rbias_gi_var <- simres_summdf %>%
  filter(stat == "rbias" & grepl("sigma", name) & trt == "gi") |>
  mutate(name_label = fct_rev(factor(name_label, levels = c("sigma[alpha]^2", "sigma[phi]^2", "sigma[epsilon]^2")))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_lbl, ncol = 1) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Relative Bias (95% C.I.)", color = "")
p_rbias_gi_var
ggsave(p_rbias_gi_var, filename = here("figures/p_rbias_gi_var_neutral.pdf"), device = cairo_pdf, width = 5, height = 4)

# Table for the same quantities
#
sink(file = glue("tables/tbias_ci_var_neutral.tex"))
simres_summdf %>%
  filter(stat == "bias" & grepl("sigma", name) & trt == "ci") |>
  mutate(name_tex = factor(name_tex, levels = c("$\\sigma_{\\alpha}^2$", "$\\sigma_{\\phi}^2$", "$\\sigma_{\\epsilon}^2$"))) |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(Scenario = scenario_lbl) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    linesep = "",
    caption = "Bias of variance components for the constant treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. Comparison of scenarios based on a joint model and on a neutral data-generating mechanism. \\label{tab:tbias_ci_var_neutral}",
    escape = FALSE
  ) |>
  print()
sink()
#
sink(file = glue("tables/tbias_gi_var_neutral.tex"))
simres_summdf %>%
  filter(stat == "bias" & grepl("sigma", name) & trt == "gi") |>
  mutate(name_tex = factor(name_tex, levels = c("$\\sigma_{\\alpha}^2$", "$\\sigma_{\\phi}^2$", "$\\sigma_{\\epsilon}^2$"))) |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(Scenario = scenario_lbl) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    linesep = "",
    caption = "Bias of variance components for the general time on treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. Comparison of scenarios based on a joint model and on a neutral data-generating mechanism. \\label{tab:tbias_gi_var_neutral}",
    escape = FALSE
  ) |>
  print()
sink()
#
sink(file = glue("tables/trbias_ci_var_neutral.tex"))
simres_summdf %>%
  filter(stat == "rbias" & grepl("sigma", name) & trt == "ci") |>
  mutate(name_tex = factor(name_tex, levels = c("$\\sigma_{\\alpha}^2$", "$\\sigma_{\\phi}^2$", "$\\sigma_{\\epsilon}^2$"))) |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(Scenario = scenario_lbl) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    linesep = "",
    caption = "Relative bias of variance components for the constant treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. Comparison of scenarios based on a joint model and on a neutral data-generating mechanism. \\label{tab:trbias_ci_var_neutral}",
    escape = FALSE
  ) |>
  print()
sink()
#
sink(file = glue("tables/trbias_gi_var_neutral.tex"))
simres_summdf %>%
  filter(stat == "rbias" & grepl("sigma", name) & trt == "gi") |>
  mutate(name_tex = factor(name_tex, levels = c("$\\sigma_{\\alpha}^2$", "$\\sigma_{\\phi}^2$", "$\\sigma_{\\epsilon}^2$"))) |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(Scenario = scenario_lbl) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    linesep = "",
    caption = "Relative bias of variance components for the general time on treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. Comparison of scenarios based on a joint model and on a neutral data-generating mechanism. \\label{tab:trbias_gi_var_neutral}",
    escape = FALSE
  ) |>
  print()
sink()

# Coverage of variance components
p_cover_ci_var <- simres_summdf %>%
  filter(stat == "cover" & grepl("sigma", name) & trt == "ci") |>
  mutate(name_label = fct_rev(factor(name_label, levels = c("sigma[alpha]^2", "sigma[phi]^2", "sigma[epsilon]^2")))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_lbl, ncol = 1) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Coverage Probability (95% C.I.)", color = "")
p_cover_ci_var
ggsave(p_cover_ci_var, filename = here("figures/p_cover_ci_var_neutral.pdf"), device = cairo_pdf, width = 5, height = 4)
#
p_cover_gi_var <- simres_summdf %>%
  filter(stat == "cover" & grepl("sigma", name) & trt == "gi") |>
  mutate(name_label = fct_rev(factor(name_label, levels = c("sigma[alpha]^2", "sigma[phi]^2", "sigma[epsilon]^2")))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_lbl, ncol = 1) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Coverage Probability (95% C.I.)", color = "")
p_cover_gi_var
ggsave(p_cover_gi_var, filename = here("figures/p_cover_gi_var_neutral.pdf"), device = cairo_pdf, width = 5, height = 4)

# Table for the same quantities
#
sink(file = glue("tables/tcover_ci_var_neutral.tex"))
simres_summdf %>%
  filter(stat == "cover" & grepl("sigma", name) & trt == "ci") |>
  mutate(name_tex = factor(name_tex, levels = c("$\\sigma_{\\alpha}^2$", "$\\sigma_{\\phi}^2$", "$\\sigma_{\\epsilon}^2$"))) |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(Scenario = scenario_lbl) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    linesep = "",
    caption = "Coverage probability of variance components for the constant treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Comparison of scenarios based on a joint model and on a neutral data-generating mechanism. \\label{tab:tcover_ci_var_neutral}",
    escape = FALSE
  ) |>
  print()
sink()
#
sink(file = glue("tables/tcover_gi_var_neutral.tex"))
simres_summdf %>%
  filter(stat == "cover" & grepl("sigma", name) & trt == "gi") |>
  mutate(name_tex = factor(name_tex, levels = c("$\\sigma_{\\alpha}^2$", "$\\sigma_{\\phi}^2$", "$\\sigma_{\\epsilon}^2$"))) |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(Scenario = scenario_lbl) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    linesep = "",
    caption = "Coverage probability of variance components for the general time on treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Comparison of scenarios based on a joint model and on a neutral data-generating mechanism. \\label{tab:tcover_gi_var_neutral}",
    escape = FALSE
  ) |>
  print()
sink()

# Bias of ICCs
p_bias_ci_icc <- simres_summdf %>%
  filter(stat == "bias" & grepl("icc", name) & trt == "ci") |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::comma, breaks = c(-0.3, -0.2, -0.1, 0.0)) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_lbl, ncol = 1) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Bias (95% C.I.)", color = "")
p_bias_ci_icc
ggsave(p_bias_ci_icc, filename = here("figures/p_bias_ci_icc_neutral.pdf"), device = cairo_pdf, width = 5, height = 4)
#
p_bias_gi_icc <- simres_summdf %>%
  filter(stat == "bias" & grepl("icc", name) & trt == "gi") |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::comma, breaks = c(-0.3, -0.2, -0.1, 0.0)) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_lbl, ncol = 1) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Bias (95% C.I.)", color = "")
p_bias_gi_icc
ggsave(p_bias_gi_icc, filename = here("figures/p_bias_gi_icc_neutral.pdf"), device = cairo_pdf, width = 5, height = 4)
#
p_rbias_ci_icc <- simres_summdf %>%
  filter(stat == "rbias" & grepl("icc", name) & trt == "ci") |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_lbl, ncol = 1) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Relative Bias (95% C.I.)", color = "")
p_rbias_ci_icc
ggsave(p_rbias_ci_icc, filename = here("figures/p_rbias_ci_icc_neutral.pdf"), device = cairo_pdf, width = 5, height = 4)
#
p_rbias_gi_icc <- simres_summdf %>%
  filter(stat == "rbias" & grepl("icc", name) & trt == "gi") |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_lbl, ncol = 1) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Relative Bias (95% C.I.)", color = "")
p_rbias_gi_icc
ggsave(p_rbias_gi_icc, filename = here("figures/p_rbias_gi_icc_neutral.pdf"), device = cairo_pdf, width = 5, height = 4)

# Table for the same quantities
#
sink(file = glue("tables/tbias_ci_icc_neutral.tex"))
simres_summdf %>%
  filter(stat == "bias" & grepl("icc", name) & trt == "ci") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(Scenario = scenario_lbl) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    linesep = "",
    caption = "Bias of ICCs for the constant treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. Comparison of scenarios based on a joint model and on a neutral data-generating mechanism. \\label{tab:tbias_ci_icc_neutral}",
    escape = FALSE
  ) |>
  print()
sink()
#
sink(file = glue("tables/tbias_gi_icc_neutral.tex"))
simres_summdf %>%
  filter(stat == "bias" & grepl("icc", name) & trt == "gi") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(Scenario = scenario_lbl) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    linesep = "",
    caption = "Bias of ICCs for the general time on treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. Comparison of scenarios based on a joint model and on a neutral data-generating mechanism. \\label{tab:tbias_gi_icc_neutral}",
    escape = FALSE
  ) |>
  print()
sink()
#
sink(file = glue("tables/trbias_ci_icc_neutral.tex"))
simres_summdf %>%
  filter(stat == "rbias" & grepl("icc", name) & trt == "ci") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(Scenario = scenario_lbl) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    linesep = "",
    caption = "Relative bias of ICCs for the constant treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. Comparison of scenarios based on a joint model and on a neutral data-generating mechanism. \\label{tab:trbias_ci_icc_neutral}",
    escape = FALSE
  ) |>
  print()
sink()
#
sink(file = glue("tables/trbias_gi_icc_neutral.tex"))
simres_summdf %>%
  filter(stat == "rbias" & grepl("icc", name) & trt == "gi") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(y = ifelse(lower > 0 | upper < 0, glue("\\textbf{[y]}", .open = "[", .close = "]"), y)) |>
  mutate(Scenario = scenario_lbl) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    linesep = "",
    caption = "Relative bias of ICCs for the general time on treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Statistically significant biases are highlighted in bold. Comparison of scenarios based on a joint model and on a neutral data-generating mechanism. \\label{tab:trbias_gi_icc_neutral}",
    escape = FALSE
  ) |>
  print()
sink()

# Coverage of ICCs
p_cover_ci_icc <- simres_summdf %>%
  filter(stat == "cover" & grepl("icc", name) & trt == "ci") |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_lbl, ncol = 1) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Coverage Probability (95% C.I.)", color = "")
p_cover_ci_icc
ggsave(p_cover_ci_icc, filename = here("figures/p_cover_ci_icc_neutral.pdf"), device = cairo_pdf, width = 5, height = 4)
#
p_cover_gi_icc <- simres_summdf %>%
  filter(stat == "cover" & grepl("icc", name) & trt == "gi") |>
  mutate(name_label = fct_rev(factor(name_label))) |>
  ggplot(aes(x = name_label, y = est, color = model)) +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 1 / 2)) +
  geom_point(position = position_dodge(width = 1 / 2)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_rda +
  theme_rda +
  facet_wrap(~scenario_lbl, ncol = 1) +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = "", y = "Coverage Probability (95% C.I.)", color = "")
p_cover_gi_icc
ggsave(p_cover_gi_icc, filename = here("figures/p_cover_gi_icc_neutral.pdf"), device = cairo_pdf, width = 5, height = 4)

# Table for the same quantities
#
sink(file = glue("tables/tcover_ci_icc_neutral.tex"))
simres_summdf %>%
  filter(stat == "cover" & grepl("icc", name) & trt == "ci") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(Scenario = scenario_lbl) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    linesep = "",
    caption = "Coverage probability of ICCs for the constant treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Comparison of scenarios based on a joint model and on a neutral data-generating mechanism. \\label{tab:tcover_ci_icc_neutral}",
    escape = FALSE
  ) |>
  print()
sink()
#
sink(file = glue("tables/tcover_gi_icc_neutral.tex"))
simres_summdf %>%
  filter(stat == "cover" & grepl("icc", name) & trt == "gi") |>
  mutate(y = glue("{comma(est, 3)} ({comma(lower, 3)}, {comma(upper, 3)})")) |>
  mutate(Scenario = scenario_lbl) |>
  select(scenario, Scenario, name_tex, modelshort, y) |>
  pivot_wider(names_from = "modelshort", values_from = "y") |>
  arrange(scenario, name_tex) |>
  select(-scenario) |>
  rename(Parameter = name_tex) |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    linesep = "",
    caption = "Coverage probability of ICCs for the general time on treatment effect parametrisation, with 95\\% confidence intervals based on Monte Carlo errors. LMM denotes the linear mixed model, while JM denotes the joint model. Comparison of scenarios based on a joint model and on a neutral data-generating mechanism. \\label{tab:tcover_gi_icc_neutral}",
    escape = FALSE
  ) |>
  print()
sink()
