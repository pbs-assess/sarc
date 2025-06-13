library(tidyverse)
library(stringr)
library(scales)
library(flextable)

library(brms)
library(bayesplot)
library(tidybayes)

library(patchwork)

theme_set(ggsidekick::theme_sleek(base_size = 12))

fit_dir <- here::here("data-generated", "models")
options(brms.file_refit = "on_change") # re-fit cached models if changes
dir.create("cache", showWarnings = FALSE)

e_table <- readRDS(here::here("data-generated", "encounter-spp-table-systematic-years.rds"))
main_spp <- e_table |>
  filter(`0` > 10 & `1` > 10)

# Maturity ~ Sarc Presence
# ------------------------------------------------------------------------------
ld <- readRDS("data-generated/length-dat.rds") |>
  filter(sex %in% c("female", "male"))
options(mc.cores = parallel::detectCores() - 2)
fit_1f <- brm(
  mature ~ 0 + Intercept + length_std * sarc_presence +
    (1 + length_std * sarc_presence | species),
  family = bernoulli(),
  data = ld |> filter(sex == "female"),
  iter = 2000L,
  warmup = 1000L,
  chains = 4L,
  cores = 4L,
  backend = "cmdstanr",
  file = "cache/length-fit1f",
  seed = 293829,
  prior =
    prior(normal(0, 2), class = b) +
    prior(student_t(3, 0, 2), class = sd) +
    prior(normal(0, 10), class = b, coef = Intercept),
  control = list(max_treedepth = 12, adapt_delta = 0.8)
)

# options(mc.cores = parallel::detectCores() - 2)
fit_1m <- update(fit_1f, newdata = ld |> filter(sex == "male"), file = "cache/length-fit1m")

# Maturity ~ presence brms
# ------------------------
# prior_summary(fit_1f) # show priors
# get_variables(fit_1f) # get variable names

# Fixed effects
post_1f <- fit_1f |>
  spread_draws(b_Intercept, b_length_std, b_sarc_presence, `b_length_std:sarc_presence`) |>
  mutate(sex = "female")
post_1m <- fit_1m |>
  spread_draws(b_Intercept, b_length_std, b_sarc_presence, `b_length_std:sarc_presence`) |>
  mutate(sex = "male")
post_fe <- bind_rows(post_1f, post_1m) |>
  pivot_longer(cols = b_Intercept:`b_length_std:sarc_presence`, values_to = "fe_coef") |>
  mutate(term = gsub("^b_", "", name)) |>
  mutate(term = gsub("length_std", "Length", term),
         term = gsub("sarc_presence", "Infection", term)) |>
  mutate(term = factor(term, levels = c("Intercept", "Length", "Infection", "Length:Infection")))

# pd_length_fe <- post_fe |>
#   ggplot(aes(x = fe_coef)) +
#   facet_grid(term ~ sex) +
#   geom_vline(xintercept = 0) +
#   ggsidekick::theme_sleek(base_size = 10) +
#   geom_density(fill = "grey90") +
#   # facet_wrap(~ term, ncol = 1) +
#   coord_cartesian(xlim = c(-3, 6), ylim = c(0, 1.7), expand = FALSE) +
#   xlab("Coefficient estimate") + ylab("Posterior density") +
#   theme(strip.clip = "off")
# pd_length_fe
# ggsave(here::here("figures", "maturity-length-posterior-densities.png"), width = 4.4, height = 3.8)

# Species effects
post_re_1f <- fit_1f |> spread_draws(r_species[species, term]) |>
  mutate(sex = "female")
post_re_1m <- fit_1m |> spread_draws(r_species[species, term]) |>
  mutate(sex = "male")
post <- bind_rows(post_re_1f, post_re_1m) |>
  mutate(term = gsub("^b_", "", term)) |>
  mutate(term = gsub("length_std", "Length", term),
         term = gsub("sarc_presence", "Infection", term)) |>
  mutate(term = factor(term, levels = c("Intercept", "Length", "Infection", "Length:Infection")))

# ggplot(post, aes(x = r_species, y = term)) +
#   facet_wrap(~ sex) +
#   geom_vline(xintercept = 0) +
#   stat_halfeye()

post_spp <- left_join(post_fe, post) |>
  mutate(combined = r_species + fe_coef) |>
  mutate(species = gsub("\\.", " ", species))

# Supplemental figure - species-level coefficients
post_spp |>
  filter(term != "Intercept") |>
  mutate(species = forcats::fct_rev(species),
         species = forcats::fct_relabel(species, str_to_title),
         sex = forcats::fct_relabel(sex, str_to_title)) |>
  ggplot() +
  aes(x = combined, y = species, colour = sex, fill = sex) +
  facet_wrap(~ term, scale = "free_x") +
  geom_vline(xintercept = 0, colour = "grey50", linetype = "dotted") +
  # ggdist::stat_halfeye(.width = c(.5, .95), trim = TRUE, alpha = 0.5, size = 2,
    # linewidth = 1, position = ggstance::position_dodgev(height = 0.4)) +
  ggdist::stat_pointinterval(.width = c(0.95),
    size = 2, linewidth = 1,
    position = ggstance::position_dodgev(height = -0.5)) +
  scale_color_manual(values = c("Female" = "black", "Male" = "grey60")) +
  scale_fill_manual(values = c("Female" = "black", "Male" = "grey60")) +
  guides(color = guide_legend(title = "Sex"), fill = guide_legend(title = "Sex")) +
  ggsidekick::theme_sleek(base_size = 12) +
  xlab("Coefficient estimate") +
  theme(axis.title.y.left = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 11))
ggsave(here::here("figures", "maturity-length-species-coef.png"), width = 7.1, height = 5.5)
ggsave(here::here("figures", "maturity-length-species-coef.pdf"), width = 7.1, height = 5.5)

# ------------------------------------------------------------
# Plot length ogives
# ------------------
# Make new data for predictions
nd <- expand.grid(species = unique(ld$species),
                  fork_length = seq(min(ld$fork_length), max(ld$fork_length), length.out = 200L),
                  sarc_presence = c(0, 1))
nd$length_std <- (nd$fork_length - mean(ld$fork_length)) / sd(ld$fork_length)
nd_f <- nd |> mutate(sex = "female")
nd_m <- nd |> mutate(sex = "male")

# Function to get posterior predictions and confidence intervals
get_pp_summary <- function(fit, newdata, re_formula = NULL, probs = c(0.05, 0.5, 0.95)) {
  # Get posterior linear predictions
  p <- brms::posterior_linpred(fit, newdata = newdata, re_formula = re_formula)

  out <- newdata
  out$lwr <- apply(p, 2, quantile, probs = probs[1])
  out$est <- apply(p, 2, quantile, probs = probs[2])
  out$upr <- apply(p, 2, quantile, probs = probs[3])

  return(out)
}

# Population-level relationship
# -----------------------------
p_f_mean <- brms::posterior_linpred(fit_1f, newdata = nd, re_formula = NA)
nd_f_mean <- get_pp_summary(fit_1f, newdata = nd_f, re_formula = NA)
nd_m_mean <- get_pp_summary(fit_1m, newdata = nd_m, re_formula = NA)

out <- bind_rows(nd_f_mean, nd_m_mean) |> as_tibble()

out |>
  mutate(sarc_pres_label = factor(sarc_presence, levels = c(0, 1), labels = c("No", "Yes"))) |>
ggplot(aes(x = fork_length, y = plogis(est),
           colour = sarc_pres_label, fill = sarc_pres_label)) +
  facet_wrap(~ sex, labeller = labeller(sex = str_to_title)) +
  geom_ribbon(aes(ymin = plogis(lwr), ymax = plogis(upr)), colour = NA, alpha = 0.3) +
  geom_line() +
  coord_cartesian(expand = FALSE, xlim = c(50, 750), ylim = c(-0.08, 1.08)) +
  geom_segment(data = ld |> filter(mature == 1),
    mapping = aes(x = fork_length, y = 1.01 + 0.03 * sarc_presence, yend = 1.04 + 0.03 * sarc_presence),
    alpha = 0.6, position = position_dodge2(width = 0.5)) +
  geom_segment(data = ld |> filter(mature == 0),
    mapping = aes(x = fork_length, y = -0.04 - 0.03 * sarc_presence, yend = -0.01 - 0.03 * sarc_presence),
    alpha = 0.6, position = position_dodge2(width = 0.5)) +
  labs(x = "Fork length (mm)", y = "Probability of maturity",
       colour = "*Sarcotaces* sp. present", fill = "*Sarcotaces* sp. present") +
  theme(legend.title = ggtext::element_markdown()) +
  scale_fill_manual(values = c("No" = "grey50", "Yes" = "red")) +
  scale_colour_manual(values = c("No" = "grey50", "Yes" = "red")) +
  theme(legend.position = "top")
ggsave(here::here("figures", "maturity-length-ogive.png"), width = 7.1, height = 4.1)
ggsave(here::here("figures", "maturity-length-ogive.pdf"), width = 7.1, height = 4.1)

# ----
# Species-level
nd_f_spp <- get_pp_summary(fit_1f, newdata = nd_f, re_formula = NULL)
nd_m_spp <- get_pp_summary(fit_1m, newdata = nd_m, re_formula = NULL)

out_spp <- bind_rows(nd_f_spp, nd_m_spp) |> as_tibble()

make_fig <- function(dat) {
  .ld <- filter(ld, species %in% dat$species)
  dat |>
    mutate(sarc_pres_label = factor(sarc_presence, levels = c(0, 1), labels = c("No", "Yes"))) |>
    ggplot(aes(x = fork_length, y = plogis(est),
      colour = sarc_pres_label, fill = sarc_pres_label)) +
    facet_grid(stringr::str_to_title(species)~sex, labeller = labeller(species = function(x) stringr::str_to_title(x))) +
    geom_ribbon(aes(ymin = plogis(lwr), ymax = plogis(upr)), colour = NA, alpha = 0.3) +
    geom_line() +
    coord_cartesian(expand = FALSE, xlim = c(50, 750), ylim = c(-0.08, 1.08)) +
    geom_segment(data = .ld |> filter(mature == 1),
      mapping = aes(x = fork_length, y = 1.01 + 0.03 * sarc_presence, yend = 1.04 + 0.03 * sarc_presence),
      alpha = 0.6, position = position_dodge2(width = 0.5)) +
    geom_segment(data = .ld |> filter(mature == 0),
      mapping = aes(x = fork_length, y = -0.04 - 0.03 * sarc_presence, yend = -0.01 - 0.03 * sarc_presence),
      alpha = 0.6, position = position_dodge2(width = 0.5)) +
    labs(x = "Fork length (mm)", y = "Probability of maturity",
      colour = "*Sarcotaces* sp. present", fill = "*Sarcotaces* sp. present") +
    theme(legend.title = ggtext::element_markdown()) +
    scale_fill_manual(values = c("No" = "grey50", "Yes" = "red")) +
    scale_colour_manual(values = c("No" = "grey50", "Yes" = "red")) +
    theme(legend.position = "top",
          strip.clip = "off")
}

sp1 <- unique(out_spp$species)[1:5]
sp2 <- unique(out_spp$species)[6:10]
g1 <- out_spp |> filter(species %in% sp1) |> make_fig()
g2 <- out_spp |> filter(species %in% sp2) |> make_fig()

g1 + g2 + plot_layout(axes = "collect", guides = "collect") & theme(legend.position = "top")

# ggsave(here::here("figures", "maturity-length-ogive-by-species.png"), width =3.5, height = 10)
ggsave(here::here("figures", "maturity-length-ogive-by-species.pdf"), width = 8.5, height = 8.5)

# ------------------------------------------------------------
# Compare expected length when probability of maturity > 0.5
# ------------------------------------------------------------
get_p50 <- function(.int, .slope, .sd = sd(ld$fork_length), .mean = mean(ld$fork_length), p = 0.5) {
  xx <- -(log((1/p) - 1) + .int) / .slope
  (xx * .sd) + .mean
}

p50df <- post_spp |>
  bind_rows(post_fe |> mutate(species = "population", combined = fe_coef)) |>
  group_by(sex, species) |>
  group_split() |>
  purrr::map_dfr(\(x) {
    intercept <- filter(x, term == "Intercept") |> pull(combined)
    length_slope <- filter(x, term == "Length") |> pull(combined)
    sarc_adj <- filter(x, term == "Infection") |> pull(combined)
    sarc_interaction <- filter(x, term == "Length:Infection") |> pull(combined)
    data.frame(
      p50_sarc_0 = get_p50(intercept, length_slope),
      p50_sarc_1 = get_p50(intercept + sarc_adj, length_slope + sarc_interaction),
      species = x$species[1],
      sex = x$sex[1]
    )
  }) |> as_tibble()

p50_diff <- p50df |>
  mutate(diff = p50_sarc_1 - p50_sarc_0)

p50_diff_summary <- p50_diff |>
  group_by(species, sex) |>
  summarise(mid = median(diff),
            lwr = quantile(diff, probs = 0.05),
            upr = quantile(diff, probs = 0.95))
saveRDS(p50_diff_summary, here::here("data-generated", "p50-length-diff-summary.rds"))

p50_diff_spp_levels <- p50_diff_summary |>
  filter(sex == "female") |>
  filter(species != "population") |>
  arrange(-mid) |>
  pull("species") %>%
  c("population", .)

# Create a reusable plotting function
plot_p50_diff <- function(data, x_limits = c(-100, 100)) {
  dat <- data |>
    mutate(species = factor(species, levels = p50_diff_spp_levels),
           species = forcats::fct_relabel(species, str_to_title),
           sex = forcats::fct_relabel(sex, str_to_title))

  p <- ggplot()
  if ("population" %in% tolower(dat$species)) {
    p <- p + geom_rect(aes(ymin = -Inf, ymax = 1.5, xmin = -Inf, xmax = Inf),
                      fill = "grey92", color = NA, inherit.aes = FALSE)
  }

  p <- p +
    geom_vline(xintercept = 0, colour = "grey50", linetype = "dotted") +
    ggdist::stat_pointinterval(data = dat, aes(x = diff, y = species, colour = sex, fill = sex),
      .width = c(0.5, 0.95),
      point_size = 2,
      position = ggstance::position_dodgev(height = -0.5)) +
    theme(axis.title.y.left = element_blank()) +
    ggsidekick::theme_sleek(base_size = 12) +
    xlab("Difference in length at 50% maturity<br>if infected with *Sarcotaces* sp.<br>(mm)") +
    theme(axis.title = ggtext::element_markdown()) +
    ylab("") +
    scale_color_manual(values = c("Female" = "black", "Male" = "grey60")) +
    scale_fill_manual(values = c("Female" = "black", "Male" = "grey60")) +
    guides(color = guide_legend(title = "Sex"), fill = guide_legend(title = "Sex")) +
    scale_x_continuous(limits = x_limits, oob = scales::squish) +
    theme(legend.position = "top",
          legend.text = element_text(size = 10))
  p
}

# Plot all species
plot_p50_diff(p50_diff)
ggsave(here::here("figures", "maturity-length-p50-by-species.png"), width = 4.7, height = 6.7)
ggsave(here::here("figures", "maturity-length-p50-by-species.pdf"), width = 4.7, height = 6.7)

# Plot main species with highlighting
plot_p50_diff(p50_diff |> filter(species %in% c(main_spp$species, "population")),
  x_limits = c(-60, 60))
ggsave(here::here("figures", "maturity-length-p50-by-main-species.png"), width = 4.4, height = 5.2)
ggsave(here::here("figures", "maturity-length-p50-by-main-species.pdf"), width = 4.4, height = 5.2)

# Reference table for values in text
p50_diff_tab <- p50_diff_summary |>
  mutate(
    Median = round(mid, 1),
    `95% CI` = paste0(round(lwr, 1), " to ", round(upr, 1))
  ) |>
  select(sex, species, Median, `95% CI`) |>
  arrange(sex, -Median) |>
  as_grouped_data(groups = c("sex")) |>
  flextable() |>
  theme_booktabs() |>
  autofit() |>
  set_caption("Differece in length at 50% maturity by sex (infected - uninfected)")
saveRDS(p50_diff_tab, here::here("data-generated", "p50-length-diff-table.rds"))

# ------------------------------------------------------------------------------
# Prior and posterior predictive checks
# -----
brms::pp_check(fit_1f, ndraws = 50)

# Observations
y_obs_f <- ld |> filter(sex == "female") |> pull(mature)
# Posterior predictions
yrep_f <- posterior_predict(fit_1f)
ppc_stat(y = y_obs_f, yrep = yrep_f, stat = mean)
ppc_dens_overlay(y_obs_f, yrep_f)

# Species means
ppc_stat_grouped(y = y_obs_f, yrep = yrep_f, group = fit_1f$data$species, stat = mean)

# Compare simulated mean expectations with observed maturity proportions
plot_sim_mat <- function(fit, ld, sex = c("female", "male"), prior_only = FALSE,
  n_draws = 10) {

  sex <- match.arg(sex)

  # pp_f <- posterior_predict(fit_1f, ndraws = 20, re_formula = NULL)
  # pp_m <- posterior_predict(fit_1m, ndraws = 20, re_formula = NULL)
  breaks <- seq(min(ld$fork_length, na.rm = TRUE), max(ld$fork_length, na.rm = TRUE), length.out = 100)

  # Observed
  l_bins <- ld |>
    filter(sex == !!sex) |>
    mutate(fork_length_bin = cut(fork_length, breaks = breaks, include.lowest = TRUE)) |>
    group_by(species, fork_length_bin, sarc_presence) |>
    summarise(
      prop_mature = sum(mature) / n(),
      mean_fork_length = mean(fork_length, na.rm = TRUE),
      .groups = "drop"
    )

  # Predicted
  pp_pred <- fit |>
    add_predicted_draws(
      newdata = ld |>
        filter(sex == !!sex) |>
        select(sex, species, fork_length, sarc_presence, length_std) |>
        mutate(fork_length_bin = cut(fork_length, breaks = breaks, include.lowest = TRUE)),
      re_formula = NULL,
      ndraws = 10
    ) |>
    mutate(fork_length_bin = cut(fork_length, breaks = breaks, include.lowest = TRUE)) |>
    group_by(species, .draw, fork_length_bin, sarc_presence) |>
    summarise(
      prop_mature = sum(.prediction) / n(),
      mean_fork_length = mean(fork_length, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(.draw = as.character(.draw))

  # Combine and plot
  plot_title <- glue::glue(
    ifelse(prior_only, "Prior predictive check", "Posterior predictive check"),
    " - {sex}"
  )

  bind_rows(
    l_bins |> mutate(source = "observed", .draw = "observed"),
    pp_pred |> mutate(source = "predicted")
  ) |>
    mutate(.draw = factor(.draw, levels = c("observed", as.character(1:n_draws)))) |>
    ggplot(aes(x = mean_fork_length, y = prop_mature, colour = species)) +
    geom_point(data = ~filter(.x, source == "observed")) +
    geom_point(data = ~filter(.x, source == "predicted"), aes(group = interaction(species, .draw)),
      alpha = 1) +
    facet_grid(.draw ~ sarc_presence, scales = "free_x") +
    scale_color_brewer(palette = "Paired") +
    guides(colour = "none") +
    labs(y = "Proportion Mature", x = "Fork Length", title = plot_title)
}

# Posterior
plot_sim_mat(fit = fit_1f, sex = "female", ld = ld)
plot_sim_mat(fit = fit_1m, sex = "male", ld = ld)

# Prior
priors <- get_prior(fit_1f)
fit_1f_p <- update(fit_1f, sample_prior = "only", prior = priors)
fit_1m_p <- update(fit_1m, sample_prior = "only", prior = priors)

plot_sim_mat(fit = fit_1f_p, sex = "female", ld = ld, prior_only = TRUE)
plot_sim_mat(fit = fit_1m_p, sex = "male", ld = ld, prior_only = TRUE)

# Prior checking
# ------------------------
# priors <- (
#   # prior(normal(0, 2), class = b) +
#   prior(student_t(3, 1, 1), class = "b") +
#   prior(student_t(3, 0, 2), class = sd) +
#   prior(normal(0, 10), class = b, coef = Intercept)
# )
# fit_1f_p <- update(fit_1f, sample_prior = "only", prior = priors)
# fit_1m_p <- update(fit_1m, sample_prior = "only")

# mcmc_areas(as_draws_df(fit_1f_p), regex_pars = c("^b_"))
# mcmc_areas(as_draws_df(fit_1f_p), regex_pars = c("^sd_")) +
#   xlim(c(0, 25))

# # Get prior draws
# f_comb <- bind_rows(
#   as_draws_df(fit_1f_p) |> mutate(source = "prior"),
#   as_draws_df(fit_1f) |> mutate(source = "posterior")
# )
# f_comb |>
#   pivot_longer(cols = starts_with("b_"), names_to = "parameter") |>
#   ggplot(aes(x = value, fill = source)) +
#   geom_density(alpha = 0.5) +
#   facet_wrap(~parameter, scales = "free")
# f_comb |>
#   pivot_longer(cols = starts_with("sd_"), names_to = "parameter") |>
#   ggplot(aes(x = value, fill = source)) +
#   geom_density(alpha = 0.5) +
#   facet_wrap(~parameter, scales = "free") +
#   xlim(c(0, 20))

# m_comb <- bind_rows(
#   as_draws_df(fit_1m_p) |> mutate(source = "prior"),
#   as_draws_df(fit_1m) |> mutate(source = "posterior")
# )
# m_comb |>
#   pivot_longer(cols = starts_with("b_"), names_to = "parameter") |>
#   ggplot(aes(x = value, fill = source)) +
#   geom_density(alpha = 0.5) +
#   facet_wrap(~parameter, scales = "free")
# m_comb |>
#   pivot_longer(cols = starts_with("sd_"), names_to = "parameter") |>
#   ggplot(aes(x = value, fill = source)) +
#   geom_density(alpha = 0.5) +
#   facet_wrap(~parameter, scales = "free") +
#   xlim(c(0, 20))