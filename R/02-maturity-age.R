library(tidyverse)
library(scales)

library(glmmTMB)
library(broom.mixed)

library(brms)
library(bayesplot)
library(tidybayes)

library(patchwork)

theme_set(ggsidekick::theme_sleek())
# Maturity ~ Sarc Presence
# ------------------------------------------------------------------------------
ad <- readRDS(here::here("data-generated", "age-dat.rds"))

# TODO: check age_std calculation
nd <- expand_grid(
    species = unique(ad$species),
    specimen_age = seq(min(ad$specimen_age), max(ad$specimen_age), length.out = 100L),
    sarc_presence = c(0, 1)
  ) |>
  mutate(age_std = (specimen_age - mean(ad$specimen_age)) / sd(ad$specimen_age))

a_bins <- ad |>
  group_by(species, specimen_age) |>
  summarise(prop_mature = sum(mature == 1) / n(), n = n(), .groups = "drop")

ggplot(ad, aes(specimen_age, mature, colour = as.factor(sarc_presence))) +
  geom_jitter(width = 0, height = 0.1, alpha = 0.9) +
  facet_grid(species ~ sarc_presence) +
  geom_smooth(method = "glm", method.args = list(family = binomial()), se = FALSE) +
  theme(legend.position = "bottom")

# Maturity: glmmTMB
# -------------------
am1 <- glmmTMB(
  mature ~ age_std * sarc_presence + (1 + age_std | species), # more complex structures don't converge
  # mature ~ age_std * sarc_presence + (1 + age_std | species), # more complex structures don't converge
  family = binomial(),
  data = ad
)
summary(am1)

pred <- predict(am1, newdata = nd, se.fit = TRUE)

nd1 <- nd |>
  mutate(
    fit = plogis(pred$fit),
    lower = plogis(pred$fit - 1.96 * pred$se.fit),
    upper = plogis(pred$fit + 1.96 * pred$se.fit)
  )

ggplot(nd1, aes(x = specimen_age, y = fit)) +
  facet_wrap(~ species) +
  geom_line(aes(colour = factor(sarc_presence))) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(sarc_presence)), alpha = 0.3) +
  geom_hline(yintercept = 0.5) +
  geom_point(data = a_bins, aes(x = specimen_age, y = prop_mature))


# Age at Maturity ~ presence brms
# -------------------------------
# Constrained the slope because there are so few data when considering age
table(ad$species, ad$sarc_presence)

# fit1 <- brm(
#   mature ~ 0 + Intercept + age_std + sarc_presence +
#     (1 + age_std | species),
#   family = bernoulli(),
#   data = ad,
#   iter = 1000L,
#   warmup = 500L,
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
# dir.create("data-generated", showWarnings = FALSE)
# saveRDS(fit1, "data-generated/maturity-age-stan-model.rds")
fit1 <- readRDS("data-generated/maturity-age-stan-model.rds")

# fit_prior <- update(
#   fit1,
#   prior =
#     prior(normal(0, 5), class = b) +
#     prior(student_t(3, 0, 2), class = sd) +
#     prior(normal(0, 10), class = b, coef = Intercept),
#   sample_prior = "only",
#   recompile = TRUE  # Ensures the model is recompiled with the new settings
# )

# fit2 <- update(
#   fit1,
#   formula = mature ~ 0 + Intercept + age_std,
#   prior =
#     prior(normal(0, 5), class = b) +
#     prior(student_t(3, 0, 2), class = sd) +
#     prior(normal(0, 10), class = b, coef = Intercept),
#   recompile = TRUE  # Ensures the model is recompiled with the new settings
# )
# beepr::beep()

# Generate simulated data based on prior samples of main effect,
p <- brms::posterior_linpred(fit1, newdata = nd, re_formula = NA)
nd$lwr <- apply(p, 2, quantile, probs = 0.05)
nd$est <- apply(p, 2, median)
nd$upr <- apply(p, 2, quantile, probs = 0.95)

ad$jit <- runif(nrow(ad), 0, 0.03)
ad_plot <- ad |>
  mutate(jit = ifelse(sarc_presence == 1, jit + 0.04, jit)) |>
  mutate(mature_jit = ifelse(mature == 1, mature - jit, mature + jit))

nd |>
  mutate(sarc_pres_label = factor(sarc_presence, levels = c(0, 1), labels = c("No", "Yes"))) |>
ggplot(aes(x = specimen_age, y = plogis(est),
               colour = sarc_pres_label, fill = sarc_pres_label)) +
  # facet_wrap(~stringr::str_to_title(species)) +
  geom_ribbon(aes(ymin = plogis(lwr), ymax = plogis(upr)), colour = NA, alpha = 0.3) +
  geom_line() +
  # coord_cartesian(expand = FALSE, ylim = c(-0.01, 1.01)) +
  # geom_point(data = ad_plot |> filter(sarc_presence == 0),
  #   mapping = aes(specimen_age, mature_jit), inherit.aes = FALSE,
  #   alpha = 0.3, pch = 21, size = 1.8) +
  # geom_point(data = ad_plot |> filter(sarc_presence == 1),
  #   mapping = aes(specimen_age, mature_jit), inherit.aes = FALSE,
  #   alpha = 0.8, pch = 19, size = 2, colour = "red") +
  geom_segment(data = ad_plot |> filter(mature == 1),
    mapping = aes(x = specimen_age, y = 1.01 + 0.03 * sarc_presence, yend = 1.04 + 0.03 * sarc_presence),
    alpha = 0.8, position = position_dodge2(width = 0.5)) +
  geom_segment(data = ad_plot |> filter(mature == 0),
    mapping = aes(x = specimen_age, y = -0.04 - 0.03 * sarc_presence, yend = -0.01 - 0.03 * sarc_presence),
    alpha = 0.8, position = position_dodge2(width = 0.5)) +
  labs(x = "Age (years)", y = "Probability of maturity",
       colour = "*Sarcotaces* sp. present", fill = "*Sarcotaces* sp. present") +
  theme(legend.title = ggtext::element_markdown()) +
  scale_fill_manual(values = c("No" = "grey50", "Yes" = "red")) +
  scale_colour_manual(values = c("No" = "grey50", "Yes" = "red")) +
  theme(legend.position = "top")

ggsave("maturity-age-overall.png")

# Convert to a data frame for plotting
p <- brms::posterior_linpred(fit1, newdata = nd)
nd$lwr <- apply(p, 2, quantile, probs = 0.05)
nd$est <- apply(p, 2, median)
nd$upr <- apply(p, 2, quantile, probs = 0.95)

ad$jit <- runif(nrow(ad), 0, 0.03)
ad_plot <- ad |>
  mutate(jit = ifelse(sarc_presence == 1, jit + 0.04, jit)) |>
  mutate(mature_jit = ifelse(mature == 1, mature - jit, mature + jit))

nd |>
  mutate(sarc_pres_label = factor(sarc_presence, levels = c(0, 1), labels = c("No", "Yes"))) |>
ggplot(aes(x = specimen_age, y = plogis(est),
               colour = sarc_pres_label, fill = sarc_pres_label)) +
  facet_wrap(~stringr::str_to_title(species)) +
  geom_ribbon(aes(ymin = plogis(lwr), ymax = plogis(upr)), colour = NA, alpha = 0.3) +
  geom_line() +
  # coord_cartesian(expand = FALSE, ylim = c(-0.01, 1.01)) +
  # geom_point(data = ad_plot |> filter(sarc_presence == 0),
  #   mapping = aes(specimen_age, mature_jit), inherit.aes = FALSE,
  #   alpha = 0.3, pch = 21, size = 1.8) +
  # geom_point(data = ad_plot |> filter(sarc_presence == 1),
  #   mapping = aes(specimen_age, mature_jit), inherit.aes = FALSE,
  #   alpha = 0.8, pch = 19, size = 2, colour = "red") +
  geom_segment(data = ad_plot |> filter(mature == 1),
    mapping = aes(x = specimen_age, y = 1.01 + 0.03 * sarc_presence, yend = 1.04 + 0.03 * sarc_presence),
    alpha = 0.6, position = position_dodge2(width = 0.5)) +
  geom_segment(data = ad_plot |> filter(mature == 0),
    mapping = aes(x = specimen_age, y = -0.04 - 0.03 * sarc_presence, yend = -0.01 - 0.03 * sarc_presence),
    alpha = 0.6, position = position_dodge2(width = 0.5)) +
  labs(x = "Age (years)", y = "Probability of maturity",
       colour = "*Sarcotaces* sp. present", fill = "*Sarcotaces* sp. present") +
  theme(legend.title = ggtext::element_markdown()) +
  scale_fill_manual(values = c("No" = "grey50", "Yes" = "red")) +
  scale_colour_manual(values = c("No" = "grey50", "Yes" = "red")) +
  theme(legend.position = "top") +
  geom_point(data = a_bins, inherit.aes = FALSE,
    aes(x = specimen_age, y = prop_mature), size = 1, colour = "blue")

ggsave("maturity-age-by-species.png")

prior_summary(fit1) # just shows the list of priors used and what can be assigned?
get_variables(fit1) # get variable names

summary(fit1)

# TODO posterior predictive checks
pp_check(fit1)
pp_check(fit1, type = "error_hist", ndraws = 11, binwidth = 0.1)
pp_check(fit1, type = "stat_2d")

# Fixed effects
# post_fe <- fit1 |> spread_draws(b_Intercept, b_age_std, b_sarc_presence, `b_age_std:sarc_presence`)
post_fe <- fit1 |> spread_draws(b_Intercept, b_age_std, b_sarc_presence)
post_fe <- post_fe |>
  pivot_longer(cols = b_Intercept:`b_sarc_presence`, values_to = "fe_coef") |>
  # pivot_longer(cols = b_Intercept:`b_age_std:sarc_presence`, values_to = "fe_coef") |>
  mutate(term = gsub("^b_", "", name))

post_fe |>
  # filter(term != "b_Intercept") |>
  ggplot(aes(x = fe_coef)) +
  geom_vline(xintercept = 0) +
  ggsidekick::theme_sleek() +
  geom_density(fill = "grey90") +
  facet_wrap(~ term, ncol = 1) +
  coord_cartesian(xlim = c(-3, 6), ylim = c(0, 1.7), expand = FALSE) +
  xlab("Coefficient estimate") + ylab("Density")

# Species effects
post_re <- fit1 |> spread_draws(r_species[species, term])
post <- left_join(post_re, post_fe)

ggplot(post, aes(x = r_species, y = term)) +
  geom_vline(xintercept = 0) +
  stat_halfeye()

# plot(fit1)


# Look at species-level effects
nd1 <- nd
p <- brms::posterior_linpred(fit1, newdata = nd)
nd1$lwr <- apply(p, 2, quantile, probs = 0.05)
nd1$est <- apply(p, 2, median)
nd1$upr <- apply(p, 2, quantile, probs = 0.95)

set.seed(3)
ad$jit <- runif(nrow(ad), 0, 0.06)
ad$mature_jit <- ifelse(ad$mature == 1, ad$mature - ad$jit, ad$mature + ad$jit)

a_bins <- ad |>
  group_by(species, specimen_age) |>
  summarise(prop_mature = sum(mature == 1) / n(), .groups = "drop")

nd1 |>
  mutate(sarc_pres_label = factor(sarc_presence, levels = c(0, 1), labels = c("No", "Yes"))) |>
ggplot(aes(x = specimen_age, y = plogis(est),  ymin = plogis(lwr), ymax = plogis(upr),
               colour = sarc_pres_label, fill = sarc_pres_label)) +
  geom_ribbon(colour = NA, alpha = 0.3) +
  geom_line() +
  facet_wrap(~stringr::str_to_title(species)) +
  coord_cartesian(expand = FALSE, ylim = c(-0.01, 1.01)) +
  geom_point(data = ad |> filter(sarc_presence == 0),
    mapping = aes(specimen_age, mature_jit), inherit.aes = FALSE,
    alpha = 0.3, pch = 21, size = 1.8) +
  geom_point(data = ad |> filter(sarc_presence == 1),
    mapping = aes(specimen_age, mature_jit), inherit.aes = FALSE,
    alpha = 0.8, pch = 19, size = 2, colour = "red") +
  labs(x = "Age (years)", y = "Probability of maturity",
       colour = "*Sarcotaces* sp. present", fill = "*Sarcotaces* sp. present") +
  theme(legend.title = ggtext::element_markdown()) +
  scale_fill_manual(values = c("No" = "grey50", "Yes" = "red")) +
  scale_colour_manual(values = c("No" = "grey50", "Yes" = "red")) +
  theme(legend.position = "top") +
  geom_point(data = a_bins, inherit.aes = FALSE,
    aes(x = specimen_age, y = prop_mature), size = 2, colour = "blue")

# Question: Why are these curves SO off??

# Population-level effect
nd2 <- nd
p <- brms::posterior_linpred(fit1, newdata = nd, re_formula = NA)
nd2$lwr <- apply(p, 2, quantile, probs = 0.05)
nd2$est <- apply(p, 2, median)
nd2$upr <- apply(p, 2, quantile, probs = 0.95)

set.seed(3)
ad$jit <- runif(nrow(ad), 0, 0.06)
ad$mature_jit <- ifelse(ad$mature == 1, ad$mature - ad$jit, ad$mature + ad$jit)

nd2 |>
  mutate(sarc_pres_label = factor(sarc_presence, levels = c(0, 1), labels = c("No", "Yes"))) |>
ggplot(aes(x = specimen_age, y = plogis(est),  ymin = plogis(lwr), ymax = plogis(upr),
               colour = sarc_pres_label, fill = sarc_pres_label)) +
  geom_ribbon(colour = NA, alpha = 0.3) +
  geom_line() +
  coord_cartesian(expand = FALSE, ylim = c(-0.01, 1.01)) +
  geom_point(data = ad |> filter(sarc_presence == 0),
    mapping = aes(specimen_age, mature_jit), inherit.aes = FALSE,
    alpha = 0.3, pch = 21, size = 1.8) +
  geom_point(data = ad |> filter(sarc_presence == 1),
    mapping = aes(specimen_age, mature_jit), inherit.aes = FALSE,
    alpha = 0.8, pch = 19, size = 2, colour = "red") +
  labs(x = "Age (years)", y = "Probability of maturity",
       colour = "*Sarcotaces* sp. present", fill = "*Sarcotaces* sp. present") +
  theme(legend.title = ggtext::element_markdown()) +
  scale_fill_manual(values = c("No" = "grey50", "Yes" = "red")) +
  scale_colour_manual(values = c("No" = "grey50", "Yes" = "red")) +
  theme(legend.position = "top") +
  geom_point(data = a_bins, inherit.aes = FALSE,
    aes(x = specimen_age, y = prop_mature, shape = species), size = 2, colour = "blue")



test <- glmmTMB(
  formula = prop_mature ~ specimen_age + (1 | species),
  family = binomial(),
  weights = n,
  data = a_bins
)

test_nd <- expand.grid(
  species = unique(a_bins$species),
  specimen_age = min(a_bins$specimen_age):max(a_bins$specimen_age),
  n = 100
)
pred <- predict(test, newdata = test_nd, se.fit = TRUE)

test_nd$fit = plogis(pred$fit)
test_nd$lower = plogis(pred$fit - 1.96 * pred$se.fit)
test_nd$upper = plogis(pred$fit + 1.96 * pred$se.fit)

ggplot(test_nd, aes(x = specimen_age, y = fit)) +
  facet_wrap(~ species) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
  geom_hline(yintercept = 0.5) +
  geom_point(data = a_bins, aes(y = prop_mature))






f <- mature ~ 0 + Intercept + length_std * sarc_presence +
  (1 + length_std * sarc_presence | species)

mod_prior1 <- brm(f,
  data = length_dat,
  prior =
    prior(normal(0, 5), class = b) +
    prior(student_t(3, 0, 2), class = sd) +
    prior(normal(0, 1.5), class = b, coef = Intercept),
    # prior(normal(0, 10), class = b, coef = Intercept),

  cores = 4L,
  sample_prior = "only"
)
mod_prior
mod_prior1 # with normal(0, 1.5) on Intercept

mod_posterior <- brm(f,
  data = length_dat,
  cores = 4L,
  prior = prior(normal(0, 5), class = b) +
  prior(student_t(3, 0, 2), class = sd) +
  prior(normal(0, 10), class = b, coef = Intercept))



mat_dat <- d |> filter(specimen_sex_desc %in% c("FEMALE", "MALE"))
# Start with difference of presence/absence of sarcs at maturity levels
fit <- brm(
  formula = sarc_presence ~ mo(maturity_code) + (mo(maturity_code) | species),
  family = bernoulli(),
  data = d,
  iter = 1000L,
  warmup = 500L,
  chains = 4L,
  cores = 4L,
  backend = "cmdstanr"
  )


# TODO:
#@ for the maturity thing (figure 6)
#@ negative binomial
#@ could answer the questions of
#@ - if there are more in males than females
#@ - are there more sarcs in different categories, namely immature compared to the mature categories

#@ simulate some data from the working model and compare against the real data to
#@ do the posterior checks

#@ ogives
#@ quang has suggested probit
#@ stock assessments do different things
#@ - bin the fork length or age, calculate the proportion mature, and then plot that out and
#@ compare it to the curve from the model
#@ this would be a reasonable thing to check

#@ for me to do
#@ - would be to have sarc number predictor on length/age at maturity
#@ - glm to
#@ - check proportion of ogives

 filter(d, sarc_) |> table(d$sarc_count, d$species)

# clean up the main analysis and output 'draft figures'
# look at sarc count as continuous
# maturity model
# write a revised part of results and methods


# # Proportion of maturity stages observed
# d |> filter(specimen_sex_desc %in% c("FEMALE", "MALE")) |>
# ggplot(aes(x = maturity_code)) +
#   geom_histogram() +
#   scale_x_continuous(breaks = 1:7) +
#   facet_wrap(~ specimen_sex_desc)