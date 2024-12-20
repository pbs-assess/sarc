library(tidyverse)
library(scales)

library(glmmTMB)
library(broom.mixed)

library(brms)
library(bayesplot)
library(tidybayes)

library(patchwork)

theme_set(ggsidekick::theme_sleek())

fit_dir <- here::here("data-generated", "models")

# Maturity ~ Sarc Presence
# ------------------------------------------------------------------------------
ld <- readRDS(here::here("data-generated", "length-dat.rds"))

# ggplot(ld, aes(fork_length, mature, colour = as.factor(sarc_presence))) +
#   geom_jitter(width = 0, height = 0.1, alpha = 0.9) +
#   facet_grid(species ~ sarc_presence) +
#   geom_smooth(method = "glm", method.args = list(family = binomial()), se = FALSE) +
#   theme(legend.position = "bottom")

# Maturity: glmmTMB
# -------------------
# m1 <- glmmTMB(
#   mature ~ length_std * sarc_presence + (1 + length_std + sarc_presence | species),
#   family = binomial(),
#   data = ld
# )
# summary(m1)

# m2 <- glmmTMB(
#   mature ~ length_std * sarc_count + (1 + length_std + sarc_count | species),
#   family = binomial(),
#   data = ld
# )
# summary(m2)

# nd <- expand_grid(
#     species = unique(ld$species),
#     fork_length = seq(min(ld$fork_length), max(ld$fork_length), length.out = 100L),
#     sarc_presence = c(0, 1),
#     sarc_count = min(ld$sarc_count, na.rm = T):max(ld$sarc_count, na.rm = T)
#   ) |>
#   mutate(length_std = (fork_length - mean(fork_length)) / sd(fork_length))

# nd2 <- nd

# pred <- predict(m2, newdata = nd2, re.form = NA,
#                 type = "link", se.fit = TRUE)
# nd2 <- nd2 |>
#   mutate(
#     fit = plogis(pred$fit),
#     lower = plogis(pred$fit - 1.96 * pred$se.fit),
#     upper = plogis(pred$fit + 1.96 * pred$se.fit)
#   )

# ggplot(nd2, aes(x = fork_length, y = fit, colour = factor(sarc_presence))) +
#   geom_line() +
#   geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(sarc_presence)), alpha = 0.3) +
#   facet_wrap(~ species) +
#   geom_vline(xintercept = mean(ld$fork_length))

# nd2 |>
#   filter(species %in% c("yelloweye rockfish", "rougheye/blackspotted",
#     "yellowmouth rockfish", "silvergray rockfish",
#     "pacific ocean perch")) |>
# ggplot(aes(x = fork_length, y = fit, colour = sarc_count, group = sarc_count)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = lower, ymax = upper, fill = sarc_count), alpha = 0.3) +
#   geom_vline(xintercept = mean(ld$fork_length)) +
#   facet_grid(sarc_count ~ species)

# Looks like there isn't enough information to look at the effect within species,
# but it does seem like the mean effect is as we might expect, with higher
# sarc infection intensity, the effect on maturity ogives leads to a reduction
# in the proportion of fish at mean length that are mature
# - whether that is because they die when they have high infection rates, only
# young fish are infected at high rates, can't exactly say

# TODO - add sex make a new model for each sex --> do for all models
# TODO: only 1999-2000; 2019:2022

# -------------------
# Maturity: brms
# -------------------
ld <- readRDS("data-generated/length-dat.rds")
fit_1f <- brm(
  mature ~ 0 + Intercept + length_std * sarc_presence +
    (1 + length_std * sarc_presence | species),
  family = bernoulli(),
  data = ld |> filter(sex == "female"),
  iter = 1000L,
  warmup = 500L,
  chains = 4L,
  cores = 4L,
  backend = "cmdstanr",
  prior =
    prior(normal(0, 5), class = b) +
    prior(student_t(3, 0, 2), class = sd) +
    prior(normal(0, 10), class = b, coef = Intercept),
  control = list(max_treedepth = 12, adapt_delta = 0.9)
)
dir.create("data-generated", showWarnings = FALSE)
saveRDS(fit_1f, file.path(fit_dir, "maturity-stan-model-female.rds"))
beepr::beep()
fit_1f <- readRDS(file.path(fit_dir, "maturity-stan-model-female.rds"))

fit_1m <- update(fit_1f, newdata = ld |> filter(sex == "female"))
saveRDS(fit_1m, file.path(fit_dir, "maturity-stan-model-male.rds"))
beepr::beep()
fit_1m <- readRDS(file.path(fit_dir, "maturity-stan-model-male.rds"))


# Maturity ~ presence brms
# ------------------------
prior_summary(fit_1f)
get_variables(fit_1f)

# Fixed effects
post_fe <- fit_1f |> spread_draws(b_Intercept, b_length_std, b_sarc_presence, `b_length_std:sarc_presence`)
post_fe <- post_fe |>
  pivot_longer(cols = b_Intercept:`b_length_std:sarc_presence`, values_to = "fe_coef") |>
  mutate(term = gsub("^b_", "", name)) |>
  mutate(term = gsub("length_std", "Length", term),
         term = gsub("sarc_presence", "Infection", term)) |>
  mutate(term = factor(term, levels = c("Intercept", "Length", "Infection", "Length:Infection")))

pd_length_fe <- post_fe |>
  ggplot(aes(x = fe_coef)) +
  geom_vline(xintercept = 0) +
  ggsidekick::theme_sleek() +
  geom_density(fill = "grey90") +
  facet_wrap(~ term, ncol = 1) +
  coord_cartesian(xlim = c(-3, 6), ylim = c(0, 1.7), expand = FALSE) +
  xlab("Coefficient estimate") + ylab("Posterior density") +
  ggtitle("Maturity ~ length: Females")
pd_length_fe
ggsave(here::here("figures", "maturity-length-posterior-densities.png"))

# Species effects
post_re <- fit1 |> spread_draws(r_species[species, term])
post <- left_join(post_re, post_fe)

ggplot(post, aes(x = r_species, y = term)) +
  geom_vline(xintercept = 0) +
  stat_halfeye()
# Question attempted interpretation:
# - wide spread of Intercept across r_species suggests that the most variation
#   in species effect is in the proportion at maturity at the mean fork length
#   in the dataset?
# - species effect next greatest on how much sarc presence affects maturity ogive

get_variables(fit1) # get variable names
prior_summary(fit1) # just shows the list of priors used and what can be assigned?
# plot(fit1)
tidy(fit1)

nd <- expand.grid(species = unique(ld$species),
                  fork_length = seq(min(ld$fork_length), max(ld$fork_length), length.out = 200L),
                  sarc_presence = c(0, 1))
nd$length_std <- (nd$fork_length - mean(ld$fork_length)) / sd(ld$fork_length)

# Mean relationship
p <- brms::posterior_linpred(fit1, newdata = nd, re_formula = NA)
nd$lwr <- apply(p, 2, quantile, probs = 0.05)
nd$est <- apply(p, 2, median)
nd$upr <- apply(p, 2, quantile, probs = 0.95)

nd |>
  mutate(sarc_pres_label = factor(sarc_presence, levels = c(0, 1), labels = c("No", "Yes"))) |>
ggplot(aes(x = fork_length, y = plogis(est),
           colour = sarc_pres_label, fill = sarc_pres_label)) +
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
ggsave(here::here("figures", "maturity-length-infection.png"), width = 4.2, height = 4.3)

# ----

p <- brms::posterior_linpred(fit1, newdata = nd)
nd$lwr <- apply(p, 2, quantile, probs = 0.05)
nd$est <- apply(p, 2, median)
nd$upr <- apply(p, 2, quantile, probs = 0.95)

# reminder to self, fork length is standardised

set.seed(3)

ld$jit <- runif(nrow(ld), 0, 0.03)
ld <- ld |>
  mutate(jit = ifelse(sarc_presence == 1, jit + 0.04, jit)) |>
  mutate(mature_jit = ifelse(mature == 1, mature - jit, mature + jit))

l_bins <- ld |>
  group_by(species, fork_length) |>
  summarise(prop_mature = sum(mature == 1) / n(),
            prop_inf = sum(sarc_presence) / n(),
            .groups = "drop")

nd |>
  mutate(sarc_pres_label = factor(sarc_presence, levels = c(0, 1), labels = c("No", "Yes"))) |>
ggplot(aes(x = fork_length, y = plogis(est),
           colour = sarc_pres_label, fill = sarc_pres_label)) +
  facet_wrap(~stringr::str_to_title(species), ncol = 3) +
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
  # geom_point(data = l_bins, inherit.aes = FALSE,
  #   aes(x = fork_length, y = prop_mature, size = sqrt(prop_inf)),
  #   colour = "blue", shape = 21,  alpha = 0.3) +
  # geom_point(data = l_bins, inherit.aes = FALSE, aes(x = fork_length, y = prop_mature), colour = "blue", alpha = 0.3) +
  theme(legend.position = "top")
ggsave(here::here("figures", "maturity-length-infection-by-species.png"), width = 7.5, height = 10.5)


# Maturity ~ count brms
# ------------------------
# Not really enough data to look at species effects of sarc count on maturity,
# but I think there is enough to make a general statement about what it does on
# average - that a higher intensity of infection does likely have a negative
# effect on maturity at length/age
# TODO: filter out n_encounters == 0
# TODO look
# TODO plot residuals to make sure that the raw sarc counts are better, might need to have log(sarc_count)
e_spp <- readRDS(here::here("data-generated", "clean-data-encounter-summary.rds"))
spp_filter <- e_spp |> filter(n_encounters > 0) |> pull(species)
cld <- ld |>
  filter(species %in% spp_filter) |>
  # filter(year >= 2019) |>
  drop_na(sarc_count)

fit_2f_lc <- brm(
  mature ~ 0 + Intercept + length_std * log(sarc_count + 0.1) +
    (1 + length_std * log(sarc_count + 0.1) | species),
  family = bernoulli(),
  data =
    cld |>
      filter(sex == "female"),
  iter = 1000L,
  warmup = 500L,
  chains = 4L,
  cores = 4L,
  backend = "cmdstanr",
  prior =
    prior(normal(0, 5), class = b) +
    prior(student_t(3, 0, 2), class = sd) +
    prior(normal(0, 10), class = b, coef = Intercept),
  control = list(max_treedepth = 12, adapt_delta = 0.95)
)
saveRDS(fit_2f_lc, file.path(fit_dir, "maturity-count-stan-model-female-lc.rds"))
beepr::beep()
fit_2f_lc <- readRDS(file.path(fit_dir, "maturity-count-stan-model-female-lc.rds"))

library(bayesplot)
# fit <- fit_2f_lc
fit <- fit_2f

test <- fit$data
fitted_values <- posterior_epred(fit)
predicted_values <- posterior_predict(fit)
observed_values <- fit$data$mature
mean_fitted <- colMeans(fitted_values)
mean_observed <- mean(observed_values)

test$fitted_prob <- colMeans(fitted_values)
ggplot(test, aes(x = fitted_prob, y = mature)) +
  geom_jitter(height = 0.05, width = 0) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE) +
  labs(x = "Fitted Probabilities", y = "Observed Outcomes",
       title = "Fitted vs Observed") +
  theme_minimal()

test$fitted_prob <- colMeans(fitted_values)
test$bin <- cut(test$fitted_prob, breaks = seq(0, 1, by = 0.1))

calibration <- test %>%
  group_by(bin) %>%
  summarize(mean_pred = mean(fitted_prob),
            mean_obs = mean(mature))

ggplot(calibration, aes(x = mean_pred, y = mean_obs)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Mean Predicted Probability", y = "Mean Observed Outcome",
       title = "Calibration Plot") +
  theme_minimal()


test <- cld |>
  filter(year >= 2019) |>
  filter(sex == "female") |>
  add_residual_draws(fit, ndraws = 100)

# LOO
loo_test <- loo(fit, save_psis = TRUE)

   |>
  ggplot(aes(x = .row, y = .residual)) +
  stat_pointinterval()

resid <- residuals(fit)
plot(resid)

pp_check(fit)
pp_check(fit, type = "hist", nsamples = 100)


test <- ld |>
  add_predicted_draws(fit_2f, ndraws = 100) |>
  summarise(
    p_lower = mean(.prediction < y_lower),
    p_upper = mean(.prediction < y_upper),
    p_residual = runif(1, p_lower, p_upper),
    z_residual = qnorm(p_residual),
    .groups = "drop_last"
  ) |>
  ggplot(aes(x = .row, y = z_residual)) +
  geom_point()

#
table(ld$sarc_count) |>
enframe() |>
pivot_wider(names_from = name, values_from = value) |>
flextable::flextable()

ld |>
  filter(sarc_presence == 1) |>
  ggplot(aes(x = species, y = sarc_count)) +
  geom_jitter(width = 0.2, height = 0) +
  scale_y_continuous(breaks = 1:10) +
  theme_minimal() +
  coord_flip() +
  labs(x = "Species", y = "Sarc Count")

# Fixed effects
post_fe2 <- fit2 |> spread_draws(b_Intercept, b_length_std, b_sarc_count, `b_length_std:sarc_count`)
post_fe2 <- post_fe2 |>
  pivot_longer(cols = b_Intercept:`b_length_std:sarc_count`, values_to = "fe_coef") |>
  mutate(term = gsub("^b_", "", name)) |>
  mutate(term = gsub("length_std", "Length", term),
         term = gsub("sarc_count", "Sarcotaces count", term)) |>
  mutate(term = factor(term, levels = c("Intercept", "Length", "Sarcotaces count", "Length:Sarcotaces count")))

pd_count_fe <- post_fe2 |>
  ggplot(aes(x = fe_coef)) +
  ggsidekick::theme_sleek() +
  geom_density(fill = "grey90") +
  geom_vline(xintercept = 0) +
  facet_wrap(~ term, ncol = 1) +
  coord_cartesian(xlim = c(-3, 6), ylim = c(0, 1.7), expand = FALSE) +
  xlab("Coefficient estimate") + ylab("Posterior density")
pd_count_fe


(pd_length_fe + ggtitle("Mature ~ Length * Infection")) +
(pd_count_fe +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank()) +
  ggtitle("Mature ~ Length * Sarc count"))

ggsave(here::here("figures", "maturity-length-infection-count-posterior-density.png"),
  width = 7.5, height = 6)


# Mean relationship ------------
nd <- expand.grid(species = unique(ld$species),
                  fork_length = seq(min(ld$fork_length), max(ld$fork_length), length.out = 200L),
                  sarc_count = 0:10)
nd$length_std <- (nd$fork_length - mean(ld$fork_length)) / sd(ld$fork_length)

p <- brms::posterior_linpred(fit2, newdata = nd, re_formula = NA)
nd$lwr <- apply(p, 2, quantile, probs = 0.05)
nd$est <- apply(p, 2, median)
nd$upr <- apply(p, 2, quantile, probs = 0.95)

nd |>
ggplot(aes(x = fork_length, y = plogis(est),
           colour = sarc_count, fill = sarc_count, group = sarc_count)) +
  geom_ribbon(aes(ymin = plogis(lwr), ymax = plogis(upr)), colour = NA, alpha = 0.3) +
  geom_line() +
  coord_cartesian(expand = FALSE, xlim = c(50, 750), ylim = c(-0.08, 1.08)) +
  geom_segment(data = ld |> filter(mature == 1, !is.na(sarc_count)),
    mapping = aes(x = fork_length, y = 1.01 + 0.03 * sarc_presence, yend = 1.04 + 0.03 * sarc_presence),
    alpha = 0.6, position = position_dodge2(width = 0.5)) +
  geom_segment(data = ld |> filter(mature == 0, !is.na(sarc_count)),
    mapping = aes(x = fork_length, y = -0.04 - 0.03 * sarc_presence, yend = -0.01 - 0.03 * sarc_presence),
    alpha = 0.6, position = position_dodge2(width = 0.5)) +
  labs(x = "Fork length (mm)", y = "Probability of maturity",
       colour = "*Sarcotaces* sp. present", fill = "*Sarcotaces* sp. present") +
  theme(legend.title = ggtext::element_markdown()) +
  theme(legend.position = "top") +
  geom_vline(xintercept = mean(ld$fork_length))
ggsave(here::here("figures", "maturity-length-sarc-count.png"), width = 4.2, height = 4.3)

# Species level -----------------
l_bins <- ld |>
  group_by(species, fork_length) |>
  summarise(prop_mature = sum(mature == 1) / n(),
            prop_inf = sum(sarc_presence) / n(),
            .groups = "drop")

# Not very different
nd |>
ggplot(aes(x = fork_length, y = plogis(est),
           colour = sarc_count, fill = sarc_count, group = sarc_count)) +
  facet_wrap(~ species) +
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

ggsave(here::here("figures", "maturity-length-count-by-species.png"), width = 7.5, height = 10.5)

# ------------------------------------------------------------------------------

# On this centered intercept, specifying a prior is actually much easier and
# intuitive than on the original intercept, since the former represents the
# expected response value when all predictors are at their means. To treat the
# intercept as an ordinary population-level effect and avoid the centering
# parameterization, use 0 + Intercept on the right-hand side of the model formula.


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
d |> filter(specimen_sex_desc %in% c("FEMALE", "MALE")) |>
ggplot(aes(x = maturity_code)) +
  geom_histogram() +
  scale_x_continuous(breaks = 1:7) +
  facet_wrap(~ specimen_sex_desc)



# For each simulated dataset, you can bin the length and calculate the proportion
# of fish that have

# Last model: (2 models)
# Number of sarcs ~ maturity bins 1-7
# Number of sarcs ~ immature/mature
#

fit <- brm(
  formula = sarc_count ~ mo(maturity_code) + (mo(maturity_code) | species),
  family = (),
  data = d,
  iter = 1000L,
  warmup = 500L,
  chains = 4L,
  cores = 4L,
  backend = "cmdstanr"
  )
