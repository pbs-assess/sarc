library(tidyverse)
library(scales)

library(glmmTMB)
library(broom.mixed)

library(brms)
library(bayesplot)
library(tidybayes)

library(patchwork)

theme_set(ggsidekick::theme_sleek())


d <- readRDS(here::here("data-generated", "clean-data.rds"))
unique(d$year)

fig_dir <- here::here("figures")
dir.create(fig_dir, showWarnings = FALSE)

#
sort(table(d$maturity_code))
d |> count(species) |> arrange(n) |> print(n = 30)
table(d$mature)


# How does infection prevalence vary across species
# ------------------------------------------------------------------------------
# What are the sample sizes for each species
e_spp <- d |> group_by(species) |>
  summarise(prop = sum(sarc_presence) / n(),
            n = n(),
            n_encounters = sum(sarc_presence)) |>
  arrange(n, prop) |>
  print(n = 23)

e_dat <- filter(d, species %in% filter(e_spp, n > 500)$species)

# Encounter: glmmTMB
# -------------------
e_fit <- glmmTMB(
  formula = sarc_presence ~ 1 + (1 | species),
  family = binomial(),
  data = e_dat
)
summary(e_fit)

e_tmb_fixed <- tidy(e_fit, effects = "fixed")
spp_means <- tidy(e_fit, effects = "ran_vals") |>
  mutate(intercept = e_tmb_fixed[["estimate"]]) |>
  mutate(re = intercept + estimate) |>
  # Question:
  # Does this work ok-ish??? Where does the delta method we keep talking about
  # at the stock renewal discussion groups come into play???
  mutate(
    logit_lower = re - 1.96 * std.error,
    logit_upper = re + 1.96 * std.error,
    prob_estimate = plogis(re),
    prob_lower = plogis(logit_lower),
    prob_upper = plogis(logit_upper)
  )

ggplot(spp_means, aes(x = reorder(level, prob_estimate), y = prob_estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = prob_lower,
                    ymax = prob_upper), width = 0.2) +
  coord_flip() +
  labs(
    x = "Species",
    y = "Encounter probability"
  )

# Encounter: brms
# -------------------
e_fit <- brm(sarc_presence ~ 0 + Intercept + (1 | species),
  family = bernoulli(),
  data = d,
  iter = 2000L,
  warmup = 500L,
  chains = 4L,
  cores = 4L,
  backend = "cmdstanr",
  prior =
    # prior(normal(0, 5), class = b) +
    prior(student_t(3, 0, 2), class = sd) +
    prior(normal(0, 10), class = b, coef = Intercept),
  # Question:
  # increased adapt_delta as suggested by output, which I now see 0.95 is the default value
  # The higher this value is, the smaller steps, so that's why you set it to 0.9
  # in the scratch file to try and make things go faster?
  # Not totally sure what step size means. Taking too big of a step makes it
  # easier to step off track? / jump off onto another slope? (whatever that really means? leaping to different part of surface?)
  # See Betancourt 2016?
  control = list(max_treedepth = 12, adapt_delta = 0.95)
  )
beepr::beep()
# saveRDS(e_fit, here::here("data-generated", "encounter-brms-n>500.rds"))
# saveRDS(e_fit, here::here("data-generated", "encounter-brms-n>100.rds"))
# saveRDS(e_fit, here::here("data-generated", "encounter-brms-n>0.rds"))

# Question:
# should we be excluding species that have few samples?
# I thought that being able to say this would help inform minimum samples needed
# to assess things and to help talk about some of the other differences between
# species effects
# Could also help justify what species to exclude in the downstream analyses?
e_fit500 <- readRDS(here::here("data-generated", "encounter-brms-n>500.rds"))
e_fit100 <- readRDS(here::here("data-generated", "encounter-brms-n>100.rds"))
e_fit0 <- readRDS(here::here("data-generated", "encounter-brms-n>0.rds"))

e_fit <- e_fit500
e_fit <- e_fit100
plot(e_fit)

e_fit |> spread_draws(b_Intercept) |>
pivot_longer(cols = b_Intercept) |>
  mutate(prob = plogis(value)) |>
  ggplot(aes(prob)) +
  ggsidekick::theme_sleek() +
  geom_density(fill = "grey90") + facet_wrap(~name, ncol = 1) +
  geom_vline(aes(xintercept = median(prob))) +
  coord_cartesian(xlim = c(0.001, 0.10), ylim = c(0, 80), expand = FALSE) +
  xlab("Coefficient estimate") + ylab("Density")
e_spp

# I think this highlights what species have a higher sarc infection prevalence
# and that it it is not just a sampling effect.
# Though you can arguably see that this is just the number of infections / samples
# But, it does give a mean expectation which could be useful to inform sampling
# minimums?
e_fit |>
  spread_draws(b_Intercept, r_species[species, term]) |>
  mutate(species = gsub("\\.", " ", species),
         species_mean = b_Intercept + r_species,
         species_prob = plogis(species_mean)) |>
  left_join(e_spp) |>
  ggplot(aes(y = forcats::fct_reorder(factor(str_to_title(species)), n), x = species_prob)) +
  # stat_pointinterval() +
  stat_halfeye() +
  geom_text(data = e_spp |> filter(n > 0) |> mutate(species = factor(str_to_title(species))),
  aes(x = -0.005, y = species, label = n), #label = paste0("(", n, ")")),
    hjust = 1, colour = "grey50", size = 2.5) +
  scale_x_continuous(labels = label_percent(), limits = c(-0.0105, 0.125)) +
  scale_y_discrete() +
  labs(
    x = "Posterior probability",
    y = "Species"
  )

# Question: I don't really know what these are doing - more self-study
pp_check(e_fit500)
pp_check(e_fit500, type = "error_hist", ndraws = 11, binwidth = 0.1)
pp_check(e_fit500, type = "stat_2d")
pp_check(e_fit500, type = "loo_pit")
ppc_error_binned(e_fit500)


# Question: What is the difference in these model specifications:
# 1) sarc_presence ~ 1 + (1 | species) # You get the same intercept value estimated? output twice???
# 2) sarc_presence ~ 0 + Intercept + (1 | species) # this is the same as above but only outputs it only once??
# 3) sarc_presence ~ 0 + (1 | species) --> does this one even work?? is this the s
# 4) sarc_presence ~ species --> this has no information sharing at all yes?
  # really is just mean(sarc_presence / n()) per species
# 1 and 2 both estimate global mean with species means sampled from normal
# distribution around global mean?
# 3 and 4 estimate only species means, would this ever be done?
c_dat <- d |>
  filter(species %in% e_spp$species[e_spp$n_encounters > 0]) |>
  # filter(species %in% e_spp$species[e_spp$n > 0]) |> # this did not fit with hurdle negbin
  filter(!is.na(sarc_count))
c_fit <- brm(sarc_count ~ 0 + Intercept + (1 | species),
  family = brms::hurdle_negbinomial(),
  data = c_dat,
  iter = 2000L,
  warmup = 500L,
  chains = 4L,
  cores = 4L,
  backend = "cmdstanr",
  prior =
    # prior(normal(0, 5), class = b) +
    prior(student_t(3, 0, 2), class = sd) +
    prior(normal(0, 10), class = b, coef = Intercept),
  # Question:
  # increased adapt_delta as suggested by output, which I now see 0.95 is the default value
  # The higher this value is, the smaller steps, so that's why you set it to 0.9
  # in the scratch file to try and make things go faster?
  # Not totally sure what step size means. Taking too big of a step makes it
  # easier to step off track? (whatever that really means? leaping to different part of surface?)
  # See Betancourt 2016?
  control = list(max_treedepth = 12, adapt_delta = 0.95)
  )
beepr::beep()
# saveRDS(c_fit, here::here("data-generated", "count-brms-hurdle-negbin.rds"))
c_fit <- readRDS(here::here("data-generated", "count-brms-zinf-negbin.rds"))
# saveRDS(c_fit, here::here("data-generated", "count-brms-e>0.rds"))
# c_fit <- readRDS(here::here("data-generated", "count-brms-e>0.rds"))

get_prior(c_fit)

summary(c_fit)
pp_check(c_fit, type = "stat_2d", ndraws = 1000) # this looks better than with just the negbin

plot(c_fit)

c_fit |> spread_draws(b_Intercept) |>
pivot_longer(cols = b_Intercept) |>
  mutate(count = exp(value)) |>
  ggplot(aes(count)) +
  ggsidekick::theme_sleek() +
  geom_density(fill = "grey90") + facet_wrap(~name, ncol = 1) +
  geom_vline(aes(xintercept = median(count))) +
  # coord_cartesian(xlim = c(0.001, 0.10), ylim = c(0, 80), expand = FALSE) +
  xlab("Coefficient estimate") + ylab("Density")

c_post <- c_fit |>
  spread_draws(b_Intercept, r_species[species, term]) |>
  mutate(
    species = gsub("\\.", " ", species),
    species_mean = b_Intercept + r_species,
    species_count = exp(species_mean)
  )

c_levels <- c_post |>
  group_by(species) |>
  summarise(mean_est = mean(species_count)) |>
  mutate(species = forcats::fct_reorder(factor(species), mean_est))

c_post |>
  left_join(e_spp) |>
  mutate(species = factor(species, levels = levels(c_levels$species))) |>
  ggplot(aes(y = species, x = species_count)) +
  geom_jitter(data = c_dat |>
    left_join(e_spp) |>
    mutate(species = factor(species, levels = levels(c_levels$species))),
    aes(x = sarc_count), colour = "grey85", shape = 21, width = 0.3) +
  stat_pointinterval() +
  geom_text(data = e_spp |> filter(n_encounters > 0) |> mutate(species = factor(species, levels = levels(c_levels$species))),
  aes(x = 0.4, label = paste0("(", n, ")")),
    hjust = 1, colour = "grey50", size = 2.5) +
  scale_x_continuous(breaks = 1:10, limits = c(0, 10)) +
  labs(
    x = "Posterior probability",
    y = "Species"
  )
# Maybe not very informative, but maybe good to show some kind of figure that
# shows that RE/BS, in particular has a high occurence of > 3 sarcs when compared
# to the other species.

c_dat |>
  left_join(e_spp) |>
  filter(sarc_count > 0) |>
  mutate(species = factor(species, levels = levels(c_levels$species))) |>
  ggplot(aes(y = species, x = sarc_count)) +
  geom_jitter(aes(x = sarc_count), colour = "black", alpha = 0.3, shape = 21, width = 0.4, height = 0.2) +
  geom_text(data = e_spp |> filter(n_encounters > 0) |> mutate(species = factor(species, levels = levels(c_levels$species))),
  aes(x = 0.4, label = paste0("(", n, ")")),
    hjust = 1, colour = "grey50", size = 2.5) +
  scale_x_continuous(breaks = 1:10, limits = c(0, 10)) +
  labs(
    x = "Number of Sarcotaces sp.",
    y = "Species"
  )
ggsave(here::here("figures", "diff-sarc-count-by-species.png"))
ggsave(here::here("figures", "diff-sarc-count-by-species.pdf"))
# ------------------------------------------------------------------------------

nd <- tibble(species = unique(c_dat$species))

p <- brms::posterior_linpred(c_fit, newdata = nd)
nd$lwr <- apply(p, 2, quantile, probs = 0.05)
nd$est <- apply(p, 2, median)
nd$upr <- apply(p, 2, quantile, probs = 0.95)

ggplot(data = nd, aes(species, y = ))
