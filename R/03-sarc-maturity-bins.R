library(tidyverse)
library(scales)
library(brms)
library(bayesplot)
library(tidybayes)
library(patchwork)

theme_set(ggsidekick::theme_sleek())

#@ for the maturity thing (figure 6)
#@ negative binomial? poisson? zero-inflated?
#@ could answer the questions of
#@ - if there are more in males than females
#@ - are there more sarcs in different categories, namely immature compared to the mature categories
  # Number of sarcs ~ maturity bins 1-7
  # Number of sarcs ~ immature/mature (can sum up posteriors of 1 + 2, and 3:7) yes?

fit_dir <- here::here("data-generated", "models")

e_spp <- readRDS(here::here("data-generated", "clean-data-encounter-summary.rds"))
spp_filter <- e_spp |> filter(n_encounters > 0) |> pull(species)

mat_lu <- readr::read_csv(here::here("data-raw", "maturity-lookup.csv")) |>
  mutate(specimen_sex_desc = toupper(sex))


dat <- readRDS(here::here("data-generated", "clean-data.rds")) |>
  filter(sex %in% c("female", "male")) |>
  filter(species %in% spp_filter) |> # only look at species with at least 1 sarc encounter
  filter(year >= 2019) |> # the number of sarcs was only counted in recent years
  filter(maturity_code %in% 1:7) |>
  left_join(mat_lu) |>
  mutate(maturity_bin = factor(maturity_code)) |>
  drop_na(sarc_count, sex, species, year, maturity_bin)


fitf <- brm(
  # sarc_count ~ 0 + Intercept + maturity_bin * sex +
  #   (1 + maturity_bin * sex | species),
  sarc_count ~ 0 + Intercept + maturity_bin +
    (1 + maturity_bin | species),
  family = zero_inflated_poisson(),
  data = dat |> filter(sex == "female"),
  iter = 2000L,
  warmup = 500L,
  chains = 4L,
  cores = 4L,
  backend = "cmdstanr",
  prior = c(prior(normal(0, 5), class = b) +
            prior(student_t(3, 0, 2), class = sd) +
            prior(normal(0, 5), class = b, coef = Intercept) +
            prior(beta(3, 1), class = zi)),  # the brms default is beta(1, 1)
    # prior(normal(0, 5), class = b) +
    # prior(student_t(3, 0, 2), class = sd) +
    # prior(normal(0, 10), class = b, coef = Intercept),
  control = list(max_treedepth = 12, adapt_delta = 0.95)
)
beepr::beep()
saveRDS(fitf, file.path(fit_dir, "sarc-number-by-group-brms-zinp-female.rds"))
# fit <- readRDS(file.path(fit_dir, "sarc-number-by-group-brms.rds"))
fitf <- readRDS(file.path(fit_dir, "sarc-number-by-group-brms-zinp-female.rds"))

fitf_prior <- update(fitf, sample_prior = "only",
  prior = c(prior(normal(0, 2), class = b) +
            prior(student_t(3, 0, 1), class = sd) +
            prior(normal(0, 2), class = b, coef = Intercept) +
            prior(beta(5, 1), class = zi)) # zi ~ a / (a + b)
)
beepr::beep()

fitm <- update(fitf, newdata = dat |> filter(sex == "male"))
saveRDS(fitm, file.path(fit_dir, "sarc-number-by-group-brms-zinp-male.rds"))
beepr::beep()


fit3 <- update(fit,
  family = zero_inflated_negbinomial()
)
beepr::beep()

get_variables(fit)

nd <- expand.grid(
  species = unique(dat$species),
  # sex = unique(dat$sex),
  maturity_bin = unique(dat$maturity_bin)
)

# Compare the difference in sarc count between male and female maturity bins
pp <- fitf |>
  add_predicted_draws(newdata = nd) # raw count predictions
  # add_epred_draws(newdata = nd)  # expected values

# The maximum predicted values are crazy huge.
test <- fitf_prior |>
  add_predicted_draws(newdata = nd)
hist(test$.prediction)
max(test$.prediction)

pp |> filter(.prediction > 0) |>
  filter(.prediction < 30) |> # just to get rid of crazy values so I can see what is happening.
ggplot(data = _, aes(x = .prediction)) +
  geom_histogram(color = "black", fill = "blue", alpha = 0.7) +
  facet_grid(maturity_bin ~ species, scales = "free")

filter(pp, species == "rougheye/blackspotted") |> arrange(-.prediction)
filter(pp, species == "rougheye/blackspotted", .prediction < 400) |> arrange(-.prediction)


ep <- fit |> add_predicted_draws(newdata = nd)

ppd <- posterior_predict(fit, draws=50)
ppc_intervals_grouped(y = dat$sarc_count, yrep = ppd, x = as.numeric(dat$maturity_bin), group = dat$species, prob = 0.5)
beepr::beep()

pp |> filter(.prediction <= 10) |>
ggplot(data = _, aes(x = .prediction)) +
  facet_grid(sex ~ maturity_bin, scales = "free") +
  geom_density()

hist(pp$.prediction)
filter(pp, .prediction < 10 & .prediction > 7)

# -------------------------------------------------------------------------------
# Drunk monk example - to play with prior predictive check
library(ggthemes)
# define parameters
prob_drink <- 0.95  # 20% of days
rate_work  <- 1    # average 1 manuscript per day

# sample one year of production
n <- 10000

# simulate days monks drink
set.seed(11)
drink <- rbinom(n, 1, prob_drink)

# simulate manuscripts completed
y <- (1 - drink) * rpois(n, rate_work)

d <-
  tibble(Y = y) %>%
  arrange(Y) %>%
  mutate(zeros = c(rep("zeros_drink", times = sum(drink)),
                   rep("zeros_work",  times = sum(y == 0 & drink == 0)),
                   rep("nope",        times = n - sum(y == 0))
                   ))

  ggplot(data = d, aes(x = Y)) +
  geom_histogram(aes(fill = zeros),
                 binwidth = 1, size = 1/10, color = "grey92") +
  scale_fill_manual(values = c(canva_pal("Green fields")(4)[1],
                               canva_pal("Green fields")(4)[2],
                               canva_pal("Green fields")(4)[1])) +
  xlab("Manuscripts completed") +
  theme_hc() +
  theme(plot.background = element_rect(fill = "grey92"),
        legend.position = "none")

monk_fit1 <- brm(data = d, family = zero_inflated_poisson,
      Y ~ 1,
      prior = c(prior(normal(0, 10), class = Intercept),
                prior(beta(2, 2), class = zi)),  # the brms default is beta(1, 1)
      cores = 4,
      seed = 11)

monk_prior <- update(monk_fit1, sample_prior = "yes")
test <- d |> add_predicted_draws(monk_prior)

filter(test, .prediction > 0) |>
  ggplot(aes(x = .prediction)) +
  geom_histogram(binwidth = 1, size = 1/10, color = "grey92") +
  scale_x_continuous(breaks = 1:10)
filter(d, Y > 0) |>
  ggplot(aes(x = Y)) +
  geom_histogram(binwidth = 1, size = 1/10, color = "grey92") +
  scale_x_continuous(breaks = 1:10, limits = c(0.5, 9.5))




ppd <- d |> add_predicted_draws(monk_fit1)

ggplot(data = ppd, aes(x = Y)) +
  geom_histogram(aes(fill = zeros),
                 binwidth = 1, size = 1/10, color = "grey92") +
  scale_fill_manual(values = c(canva_pal("Green fields")(4)[1],
                               canva_pal("Green fields")(4)[2],
                               canva_pal("Green fields")(4)[1])) +
  xlab("Manuscripts completed") +
  theme_hc() +
  theme(plot.background = element_rect(fill = "grey92"),
        legend.position = "none")

max(ppd$.prediction)
max(y)

# Understanding effect of beta prior on zi parameter
# https://discourse.mc-stan.org/t/specifying-bernoulli-prior-for-zero-inflated-beta-response-variable-in-brm-hgam/30076
a <- 3
b <- 1
a / (a + b)
