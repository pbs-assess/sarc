library(dplyr)
library(ggplot2)

library(glmmTMB)
library(broom.mixed)

library(brms)

e_spp <- readRDS(here::here("data-generated", "encounter-spp-table.rds")) |>
  filter(n_encounters > 0) |>
  pull(species)

d <- readRDS(here::here("data-generated", "clean-data.rds")) |>
  filter(!(maturity_code %in% c("0", "8", "NULL"))) # Need complete maturity code data

m_dat <- d |>
  filter(specimen_sex_desc %in% c("FEMALE", "MALE")) |>
  filter(species %in% e_spp) |>
  drop_na(maturity_code) |>
  mutate(maturity_code = factor(maturity_code))


m1 <- glmmTMB(
  sarc_count ~ 0 + maturity_code * specimen_sex_desc + (1 | species),
  family = nbinom2(),
  data = m_dat
)
summary(m1)

m_tmb_fixed <- tidy(m1, effects = "fixed")
spp_means <- tidy(m1, effects = "ran_vals") |>
  mutate(intercept = m_tmb_fixed[1, "estimate"][[1]]) |>
  mutate(re = intercept + estimate) |>
  mutate(
    lower = exp(re - 1.96 * std.error),
    upper = exp(re + 1.96 * std.error)
  )

ggplot(spp_means, aes(x = reorder(level, estimate), y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower,
                    ymax = upper), width = 0.2) +
  coord_flip() +
  labs(
    x = "Species",
    y = "Mean count"
  )

# Encounter: brms
# -------------------
mcode_fit <- brm(sarc_count ~ 0 + Intercept + maturity_code*specimen_sex_desc + (1 | species),
  family = negbinomial(),
  data = d,
  iter = 2000L,
  warmup = 500L,
  chains = 4L,
  cores = 4L,
  backend = "cmdstanr",
  prior =
    prior(normal(0, 5), class = b) +
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
summary(mcode_fit)
