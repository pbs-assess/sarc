library(tidyverse)
library(scales)

library(brms)
library(bayesplot)
library(tidybayes)

library(patchwork)

theme_set(ggsidekick::theme_sleek())

options(brms.file_refit = "on_change") # re-fit cached models if changes
dir.create("cache", showWarnings = FALSE)

fit_dir <- here::here("data-generated", "models")

# Maturity ~ Sarc Presence
# ------------------------------------------------------------------------------
ad <- readRDS(here::here("data-generated", "age-dat.rds")) |>
  filter(sex %in% c("male", "female")) |>
  mutate(sex = factor(sex))

# Note: the original MS Figure 4 had sample sizes much larger than this, but
# that is because it was including years where sampling does not appear to be
# systematic.
# test0 <- readRDS(here::here("data-generated", "clean-data-all-years.rds")) |>
#   filter(sex %in% c("male", "female")) |>
#   filter(species %in% c("yelloweye rockfish", "pacific ocean perch", "quillback rockfish", "yellowmouth rockfish"))
# test <- test0 |>
#   drop_na(specimen_age, sex, mature, sarc_presence)

# table(test$species, test$sarc_presence, test$sex)
# table(ad$species, ad$sarc_presence, ad$sex)

nd <- expand_grid(
    species = unique(ad$species),
    sex = unique(ad$sex),
    specimen_age = seq(min(ad$specimen_age), max(ad$specimen_age), length.out = 100L),
    sarc_presence = c(0, 1)
  ) |>
  droplevels() |>
  mutate(age_std = (specimen_age - mean(ad$specimen_age)) / sd(ad$specimen_age))

a_bins <- ad |>
  group_by(species, specimen_age) |>
  summarise(prop_mature = sum(mature == 1) / n(), n = n(), .groups = "drop")

# ggplot(ad, aes(specimen_age, mature, colour = as.factor(sarc_presence))) +
#   geom_jitter(width = 0, height = 0.1, alpha = 0.9) +
#   facet_grid(species ~ sarc_presence) +
#   geom_smooth(method = "glm", method.args = list(family = binomial()), se = FALSE) +
#   theme(legend.position = "bottom")

# Age at Maturity ~ presence brms
# -------------------------------
# Constrained the slope because there are so few data when considering age
table(ad$species, ad$sarc_presence, ad$sex)
table(ad$sarc_presence) |> sum()


options(mc.cores = parallel::detectCores() - 2)
# Below won't converge
# fit <- brm(
#   mature ~ 0 + Intercept + (age_std + sarc_presence) * sex +
#   (1 + (age_std + sarc_presence) * sex | species),
#   family = bernoulli(),
#   data = ad |> filter(species == "pacific ocean perch"),
#   iter = 4000L,
#   warmup = 1000L,
#   chains = 4L,
#   cores = 4L,
#   backend = "cmdstanr",
#   prior =
#     prior(normal(0, 5), class = b) +
#     prior(student_t(3, 0, 2), class = sd) +
#     prior(normal(0, 10), class = b, coef = Intercept),
#   control = list(max_treedepth = 12, adapt_delta = 0.9)
# )
# beepr::beep()
# # dir.create("data-generated", showWarnings = FALSE)
# saveRDS(fit, file.path(fit_dir, "maturity-age-stan-model.rds"))
# fit <- readRDS(file.path(fit_dir, "maturity-age-stan-model.rds"))


# color_scheme_set("viridis")
options(mc.cores = parallel::detectCores() - 2)
fit_1f <- brm(
  mature ~ 0 + Intercept + age_std + sarc_presence,
  family = bernoulli(),
  data = ad |> filter(sex == "female"),
  iter = 2000L,
  warmup = 1000L,
  seed = 24821,
  file = "cache/age-fit1f",
  chains = 4L,
  cores = 4L,
  backend = "cmdstanr",
  prior =
    prior(normal(0, 5), class = b) +
    prior(normal(0, 10), class = b, coef = Intercept),
  control = list(max_treedepth = 12, adapt_delta = 0.97)
)
# beepr::beep()

fit_1m <- update(fit_1f, newdata = ad |> filter(sex == "male"), file = "cache/age-fit1m")

# Including species level effects - doesn't converge anymore after cleaning up data
# fit_2f <- brm(
#   mature ~ 0 + Intercept + (age_std + sarc_presence) +
#     (1 + age_std + sarc_presence | species),
#   family = bernoulli(),
#   data = ad |> filter(sex == "female"),
#   iter = 4000L,
#   warmup = 1000L,
#   chains = 4L,
#   cores = 4L,
#   backend = "cmdstanr",
#   prior =
#     prior(normal(0, 5), class = b) +
#     prior(student_t(3, 0, 2), class = sd) +
#     prior(normal(0, 10), class = b, coef = Intercept),
#   control = list(max_treedepth = 12, adapt_delta = 0.97)
# )
# beepr::beep()
# saveRDS(fit_2f, file.path(fit_dir, "maturity-age-stan-model-female-by-species.rds"))
# fit_2f <- readRDS(file.path(fit_dir, "maturity-age-stan-model-female-by-species.rds"))

# fit_2m <- update(fit_2f, newdata = ad |> filter(sex == "male"))
# saveRDS(fit_2m, file.path(fit_dir, "maturity-age-stan-model-male-by-species.rds"))
# beepr::beep()
# fit_2m <- readRDS(file.path(fit_dir, "maturity-age-stan-model-male-by-species.rds"))

# ------------------------
# Prior checking
# ------------------------
priors <- (
  prior(normal(0, 5), class = b) +
  prior(student_t(3, 0, 2), class = sd) +
  prior(normal(0, 10), class = b, coef = Intercept)
)
fit_1f_p <- update(fit_1f, sample_prior = "only")
fit_1m_p <- update(fit_1m, sample_prior = "only")

mcmc_areas(as_draws_df(fit_1f_p), regex_pars = c("^b_"))
# mcmc_areas(as_draws_df(fit_2f_p), regex_pars = c("^sd_")) +
#   xlim(c(0, 25))

# Get prior draws
f_comb <- bind_rows(
  as_draws_df(fit_1f_p) |> mutate(source = "prior"),
  as_draws_df(fit_1f) |> mutate(source = "posterior")
)
f_comb |>
  pivot_longer(cols = starts_with("b_"), names_to = "parameter") |>
  ggplot(aes(x = value, fill = source)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~parameter, scales = "free")
# f_comb |>
#   pivot_longer(cols = starts_with("sd_"), names_to = "parameter") |>
#   ggplot(aes(x = value, fill = source)) +
#   geom_density(alpha = 0.5) +
#   facet_wrap(~parameter, scales = "free") +
#   xlim(c(0, 20))

m_comb <- bind_rows(
  as_draws_df(fit_1m_p) |> mutate(source = "prior"),
  as_draws_df(fit_1m) |> mutate(source = "posterior")
)
m_comb |>
  pivot_longer(cols = starts_with("b_"), names_to = "parameter") |>
  ggplot(aes(x = value, fill = source)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~parameter, scales = "free")
# m_comb |>
#   pivot_longer(cols = starts_with("sd_"), names_to = "parameter") |>
#   ggplot(aes(x = value, fill = source)) +
#   geom_density(alpha = 0.5) +
#   facet_wrap(~parameter, scales = "free") +
#   xlim(c(0, 20))

# Maturity ~ presence brms
# ------------------------
prior_summary(fit_1f)
get_variables(fit_1f)

# These plots aren't really informative right?
# # Fixed effects
# post_2f <- fit_2f |>
#   spread_draws(b_Intercept, b_age_std, b_sarc_presence) |>
#   mutate(sex = "female")
# post_2m <- fit_2m |>
#   spread_draws(b_Intercept, b_age_std, b_sarc_presence) |>
#   mutate(sex = "male")
# post_fe <- bind_rows(post_2f, post_2m) |>
#   pivot_longer(cols = b_Intercept:`b_sarc_presence`, values_to = "fe_coef") |>
#   mutate(term = gsub("^b_", "", name)) |>
#   mutate(term = gsub("age_std", "Age", term),
#          term = gsub("sarc_presence", "Infection", term)) |>
#   mutate(term = factor(term, levels = c("Intercept", "Age", "Infection")))

# # pd_age_fe <- post_fe |>
# #   ggplot(aes(x = fe_coef)) +
# #   facet_grid(term ~ sex) +
# #   geom_vline(xintercept = 0) +
# #   ggsidekick::theme_sleek(base_size = 10) +
# #   geom_density(fill = "grey90") +
# #   # facet_wrap(~ term, ncol = 1) +
# #   # coord_cartesian(xlim = c(-3, 6), ylim = c(0, 1.7), expand = FALSE) +
# #   xlab("Coefficient estimate") + ylab("Posterior density") +
# #   theme(strip.clip = "off")
# # pd_age_fe

# # ggsave(here::here("figures", "maturity-age-posterior-densities.png"), width = 4.4, height = 3.8)

# # Species effects
# post_re_2f <- fit_2f |> spread_draws(r_species[species, term]) |>
#   mutate(sex = "female")
# post_re_2m <- fit_2m |> spread_draws(r_species[species, term]) |>
#   mutate(sex = "male")
# post <- bind_rows(post_re_2f, post_re_2m) |>
#   mutate(term = gsub("^b_", "", term)) |>
#   mutate(term = gsub("age_std", "Age", term),
#          term = gsub("sarc_presence", "Infection", term)) |>
#   mutate(term = factor(term, levels = c("Intercept", "Age", "Infection")))

# ggplot(post, aes(x = r_species, y = term)) +
#   facet_wrap(~ sex) +
#   geom_vline(xintercept = 0) +
#   stat_halfeye()

# post_spp <- left_join(post_fe, post) |>
#   mutate(combined = r_species + fe_coef) |>
#   mutate(species = gsub("\\.", " ", species))

# post_spp |>
#   filter(term != "Intercept") |>
#   group_by(term, sex, species) |>
#   summarise(lwr = quantile(combined, probs = 0.0275), upr = quantile(combined, probs = 0.975), mid = median(combined)) |>
#   ggplot(aes(x = mid, y = species, xmin = lwr, xmax = upr)) +
#   facet_wrap(sex ~ term) +
#   geom_vline(xintercept = 0, lty = 2) +
#   geom_point(pch = 21) +
#   geom_linerange() + ggsidekick::theme_sleek() +
#   xlab("Coefficient estimate") + theme(axis.title.y.left = element_blank())
# ggsave(here::here("figures", "maturity-age-species-coef.png"))

# Overall ogives
# --------------
nd <- expand.grid(sex = unique(ad$sex),
                  species = unique(ad$species),
                  specimen_age = seq(min(ad$specimen_age), max(ad$specimen_age), length.out = 200L),
                  sarc_presence = c(0, 1))
nd$age_std <- (nd$specimen_age - mean(ad$specimen_age)) / sd(ad$specimen_age)
nd_f <- nd |> mutate(sex = "female")
nd_m <- nd |> mutate(sex = "male")

# Mean relationship
p_f_mean <- brms::posterior_linpred(fit_1f, newdata = nd_f, re_formula = NA)
nd_f_mean <- nd_f
nd_f_mean$lwr <- apply(p_f_mean, 2, quantile, probs = 0.05)
nd_f_mean$est <- apply(p_f_mean, 2, median)
nd_f_mean$upr <- apply(p_f_mean, 2, quantile, probs = 0.95)

p_m_mean <- brms::posterior_linpred(fit_1m, newdata = nd_m, re_formula = NA)
nd_m_mean <- nd_m
nd_m_mean$lwr <- apply(p_m_mean, 2, quantile, probs = 0.05)
nd_m_mean$est <- apply(p_m_mean, 2, median)
nd_m_mean$upr <- apply(p_m_mean, 2, quantile, probs = 0.95)

age_ogive <-
bind_rows(nd_f_mean, nd_m_mean) |>
  mutate(sarc_pres_label = factor(sarc_presence, levels = c(0, 1), labels = c("No", "Yes"))) |>
ggplot(aes(x = specimen_age, y = plogis(est),
           colour = sarc_pres_label, fill = sarc_pres_label)) +
  facet_wrap(~ sex, labeller = labeller(sex = str_to_title)) +
  geom_ribbon(aes(ymin = plogis(lwr), ymax = plogis(upr)), colour = NA, alpha = 0.3) +
  geom_line() +
  coord_cartesian(expand = FALSE, xlim = c(0, 75), ylim = c(-0.08, 1.08)) +
  scale_x_continuous(breaks = seq(0, 75, 10)) +
  geom_segment(data = ad |> filter(mature == 1),
    mapping = aes(x = specimen_age, y = 1.01 + 0.03 * sarc_presence, yend = 1.04 + 0.03 * sarc_presence),
    alpha = 0.6, position = position_dodge2(width = 0.5)) +
  geom_segment(data = ad |> filter(mature == 0),
    mapping = aes(x = specimen_age, y = -0.04 - 0.03 * sarc_presence, yend = -0.01 - 0.03 * sarc_presence),
    alpha = 0.6, position = position_dodge2(width = 0.5)) +
  labs(x = "Age (years)", y = "Probability of maturity",
       colour = "*Sarcotaces* sp. present", fill = "*Sarcotaces* sp. present") +
  theme(legend.title = ggtext::element_markdown()) +
  scale_fill_manual(values = c("No" = "grey50", "Yes" = "red")) +
  scale_colour_manual(values = c("No" = "grey50", "Yes" = "red")) +
  theme(legend.position = "top",
        legend.margin = margin(t = 0, r = 0, b = -5, l = 0),
        strip.text = element_text(size = 10))
age_ogive
ggsave(here::here("figures", "maturity-age-overall-ogive.pdf"), width = 7.1, height = 4.1)
# ggsave(here::here("figures", "maturity-age-overall-ogive.png"), width = 7.1, height = 4.1)
saveRDS(age_ogive, here::here("data-generated", "age-ogive-ggplot.rds"))

length_ogive <- readRDS(here::here("data-generated", "length-ogive-ggplot.rds"))

wrap_plots(list(length_ogive, age_ogive), ncol = 1,
  # axis_titles = "collect_y",
  guides = "collect") +
  plot_layout(widths = c(1, 20), guides = "collect") +
  plot_annotation(tag_levels = 'a', tag_suffix = ")") &
  theme(
    plot.tag.position = c(0.05, 0.96),
    axis.title.y = element_text(vjust = 6),
    plot.margin = margin(t = 2, r = 2, b = 2, l = 10),
    legend.position = "top",
    legend.margin = margin(t = 0, r = 0, b = -5, l = 0)
  )

ggsave(here::here("figures", "length-age-ogives.pdf"), width = 7.1, height = 6.1)
# ggsave(here::here("figures", "length-age-ogives.png"), width = 7.1, height = 6.1)

# --------
# Species-level
# nd_f <- nd |> filter(sex == "female")
# p <- brms::posterior_linpred(fit_2f, newdata = nd_f, re_formula = NULL)
# nd_f$lwr <- apply(p, 2, quantile, probs = 0.05)
# nd_f$est <- apply(p, 2, median)
# nd_f$upr <- apply(p, 2, quantile, probs = 0.95)

# nd_m <- nd |> filter(sex == "male")
# p <- brms::posterior_linpred(fit_2m, newdata = nd_m, re_formula = NULL)
# nd_m$lwr <- apply(p, 2, quantile, probs = 0.05)
# nd_m$est <- apply(p, 2, median)
# nd_m$upr <- apply(p, 2, quantile, probs = 0.95)


# ad$jit <- runif(nrow(ad), 0, 0.03)
# ad_plot <- ad |>
#   mutate(jit = ifelse(sarc_presence == 1, jit + 0.04, jit)) |>
#   mutate(mature_jit = ifelse(mature == 1, mature - jit, mature + jit))

# bind_rows(nd_f, nd_m) |>
#   mutate(sarc_pres_label = factor(sarc_presence, levels = c(0, 1), labels = c("No", "Yes"))) |>
# ggplot(aes(x = specimen_age, y = plogis(est),
#                colour = sarc_pres_label, fill = sarc_pres_label)) +
#   facet_grid(sex ~ species, scales = "free_x", labeller = labeller(sex = str_to_title, species = str_to_title)) +
#   geom_ribbon(aes(ymin = plogis(lwr), ymax = plogis(upr)), colour = NA, alpha = 0.3) +
#   geom_line() +
#   coord_cartesian(expand = FALSE, ylim = c(-0.08, 1.08)) +
#   # geom_point(data = ad_plot |> filter(sarc_presence == 0),
#   #   mapping = aes(specimen_age, mature_jit), inherit.aes = FALSE,
#   #   alpha = 0.3, pch = 21, size = 1.8) +
#   # geom_point(data = ad_plot |> filter(sarc_presence == 1),
#   #   mapping = aes(specimen_age, mature_jit), inherit.aes = FALSE,
#   #   alpha = 0.8, pch = 19, size = 2, colour = "red") +
#   geom_segment(data = ad_plot |> filter(mature == 1),
#     mapping = aes(x = specimen_age, y = 1.01 + 0.03 * sarc_presence, yend = 1.04 + 0.03 * sarc_presence),
#     alpha = 0.8, position = position_jitter(width = 0.2, height = 0)) +
#   geom_segment(data = ad_plot |> filter(mature == 0),
#     mapping = aes(x = specimen_age, y = -0.04 - 0.03 * sarc_presence, yend = -0.01 - 0.03 * sarc_presence),
#     alpha = 0.8, position = position_jitter(width = 0.2, height = 0)) +
#   labs(x = "Age (years)", y = "Probability of maturity",
#        colour = "*Sarcotaces* sp. present", fill = "*Sarcotaces* sp. present") +
#   theme(legend.title = ggtext::element_markdown()) +
#   scale_fill_manual(values = c("No" = "grey50", "Yes" = "red")) +
#   scale_colour_manual(values = c("No" = "grey50", "Yes" = "red")) +
#   theme(legend.position = "top")
# ggsave(here::here("figures", "maturity-age-ogives-by-species.png"), width = 7, height = 3.6)



# # Compare expected length when probability of maturity > 0.5
# get_p50 <- function(.int, .slope, .sd = sd(ad$specimen_age), .mean = mean(ad$specimen_age), p = 0.5) {
#   xx <- -(log((1/p) - 1) + .int) / .slope
#   (xx * .sd) + .mean
# }

# p50df <- group_by(post_spp, sex, species) |>
#   group_split() |>
#   purrr::map_dfr(\(x) {
#     intercept <- filter(x, term == "Intercept") |> pull(combined)
#     length_slope <- filter(x, term == "Length") |> pull(combined)
#     sarc_adj <- filter(x, term == "Infection") |> pull(combined)
#     sarc_interaction <- filter(x, term == "Length:Infection") |> pull(combined)
#     data.frame(
#       p50_sarc_0 = get_p50(intercept, length_slope),
#       p50_sarc_1 = get_p50(intercept + sarc_adj, length_slope + sarc_interaction),
#       species = x$species[1],
#       sex = x$sex[1]
#     )
#   }) |> as_tibble()

# p50_diff <- p50df |>
#   mutate(diff = p50_sarc_1 - p50_sarc_0) |>
#   group_by(species, sex) |>
#   summarise(
#     lwr = quantile(diff, probs = 0.05),
#     lwr2 = quantile(diff, probs = 0.25),
#     upr = quantile(diff, probs = 0.95),
#     upr2 = quantile(diff, probs = 0.75),
#     mid = median(diff)
#   )

# p50_diff |>
#   mutate(species = factor(species, levels = p50_diff |> filter(sex == "female") |> arrange(-mid) |> pull("species"))) |>
#   ggplot(aes(mid, species, xmin = lwr, xmax = upr)) +
#   geom_linerange(lwd = 0.4) +
#   geom_linerange(aes(xmin = lwr2, xmax = upr2), lwd = .7) +
#   geom_point(pch = 19) +
#   geom_vline(xintercept = 0, lty = 2)+
#   theme(axis.title.y.left = element_blank()) +
#   ggsidekick::theme_sleek() +
#   xlab("Difference in length at 50% maturity<br>if infected with *Sarcotaces* sp.<br>(mm)") +
#   theme(axis.title = ggtext::element_markdown()) +
#   ylab("") +
#   facet_grid(. ~ sex, labeller = labeller(sex = str_to_title)) +
#   scale_y_discrete(label = str_to_title) +
#   scale_x_continuous(limits = c(-100, 100), oob = scales::squish) # NOTE: copper rockfish cutoff
# ggsave(here::here("figures", "maturity-length-p50-by-species.png"), width = 8, height = 7)


# ------------------------------------------------------------------------------
# Posterior predictive checks
# -----
color_scheme_set("mix-brightblue-gray")
brms::pp_check(fit_2f, ndraws = 50)
mcmc_trace(fit_2f, regex_pars = "^b")
mcmc_trace(fit_2f, regex_pars = "^r")
mcmc_trace(fit_2f)
plot(fit_2f)

# Observations
y_obs_f <- ad |> filter(sex == "female") |> pull(mature)
# Posterior predictions
yrep_f <- posterior_predict(fit_2f)
ppc_stat(y = y_obs_f, yrep = yrep_f, stat = mean)
ppc_dens_overlay(y_obs_f, yrep_f)

# https://discourse.mc-stan.org/t/posterior-predictive-checks-kurtosis-and-skew/11136/3
skew <- function(x) {
  xdev <- x - mean(x)
  n <- length(x)
  r <- sum(xdev^3) / sum(xdev^2)^1.5
  return(r * sqrt(n) * (1 - 1/n)^1.5)
}

color_scheme_set("blue")
ppc_stat(y_obs_f, yrep_f, stat = "skew")
ppc_stat_grouped(y_obs_f, yrep_f, stat = "skew", group = fit_2f$data$species)

# Species means
ppc_stat_grouped(y = y_obs_f, yrep = yrep_f, group = fit_2f$data$species, stat = mean)

# Compare simulated mean expectations with observed maturity proportions
breaks <- seq(min(ad$specimen_age, na.rm = TRUE), max(ad$specimen_age, na.rm = TRUE), length.out = 100)

pp_f <- posterior_predict(fit_2f, ndraws = 20, re_formula = NULL)
pp_m <- posterior_predict(fit_2m, ndraws = 20, re_formula = NULL)

# Observed
a_bins <- ad |>
  mutate(age_bin = cut(specimen_age, breaks = breaks, include.lowest = TRUE)) |>
  group_by(species, sex, age_bin, sarc_presence) |>
  summarise(
    prop_mature = sum(mature) / n(),
    mean_age = mean(specimen_age, na.rm = TRUE),
    .groups = "drop"
  )

# a_bins |>
#   filter(sarc_presence == 0) |>
# ggplot(data = _, aes(x = mean_age, y = prop_mature)) +
#   geom_point() +
#   guides(colour = "none") +
#   facet_grid(sex ~ ., scales = "free_x")

# Simulated
pp_f_species <- fit_2f |>
  # add_epred_draws(
  add_predicted_draws(
    newdata = ad |>
      filter(sex == "female") |>
      select(sex, species, specimen_age, sarc_presence, age_std) |>
      mutate(age_bin = cut(specimen_age, breaks = breaks, include.lowest = TRUE)),
  re_formula = NULL, ndraws = 10)

pp_m_species <- fit_2m |>
  add_predicted_draws(
    newdata = ad |>
      filter(sex == "male") |>
      select(sex, species, specimen_age, sarc_presence, age_std) |>
      mutate(age_bin = cut(specimen_age, breaks = breaks, include.lowest = TRUE)),
  re_formula = NULL, ndraws = 10)

pred <- bind_rows(pp_f_species, pp_m_species) |>
  mutate(age_bin = cut(specimen_age, breaks = breaks, include.lowest = TRUE)) |>
  group_by(species, sex, .draw, age_bin, sarc_presence) |>
  summarise(
    # prop_mature = sum(.epred) / n(),
    prop_mature = sum(.prediction) / n(),
    mean_age = mean(specimen_age, na.rm = TRUE),
    .groups = "drop"
  )

bind_rows(a_bins, pred) |>
  filter(sex == "female") |>
ggplot(data = _, aes(x = mean_age, y = prop_mature)) +
  geom_point(aes(colour = species)) +
  guides(colour = "none") +
  scale_color_brewer(palette = "Paired") +
  facet_grid(.draw ~ sarc_presence, scales = "free_x")

bind_rows(a_bins, pred) |>
  filter(sex == "male") |>
ggplot(data = _, aes(x = mean_age, y = prop_mature)) +
  geom_point(aes(colour = species)) +
  guides(colour = "none") +
  scale_color_brewer(palette = "Paired") +
  facet_grid(.draw ~ sarc_presence, scales = "free_x")
