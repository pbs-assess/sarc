library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(broom)
library(brms)
library(bayesplot)
library(tidybayes)

theme_set(ggsidekick::theme_sleek())
dir.create("cache", showWarnings = FALSE)

options(brms.file_refit = "on_change") # re-fit cached models if changes
options(mc.cores = parallel::detectCores() - 2)

e_table <- readRDS(here::here("data-generated", "encounter-spp-table-systematic-years.rds"))
main_spp <- e_table |>
  filter(`0` > 10 & `1` > 10)

# Remove length-weight outliers
# ------------------------
if (!file.exists(here::here("data-generated", "length-weight-dat-clean.rds"))) {
  lwd0 <- readRDS(here::here("data-generated", "length-dat.rds")) |>
    filter(sex %in% c("female", "male")) |>
    select(specimen_id, year, species, sex, fork_length, round_weight, sarc_presence,
      survey_series_desc, mature) |> # why is this round weight?
  filter(!is.na(round_weight))

  # Identify outliers (should this also be done for the length analysis - it won't change the results there, but would be more consistent)
  # - I think it might be important to exclude outliers here becuase the mistakes will compound
  lw_fits <- lwd0 |>
    filter(!is.na(fork_length), !is.na(round_weight), fork_length > 0, round_weight > 0) |>
    mutate(
      log_length = log(fork_length),
      log_weight = log(round_weight)
    ) |>
    group_by(species) |>
    nest() |>
    mutate(model = purrr::map(data, ~ lm(log_weight ~ log_length, data = .x))) |>
    mutate(augmented = purrr::map2(model, data, ~ broom::augment(.x, data = .y))) |>
    unnest(augmented)

  outliers <- lw_fits |>
    group_by(species) |>
    filter(abs(.std.resid) > 10) # use SD to identify outliers

  # sidebar for SD/outlier justification and my understanding:
  # Understand probability of outliers expected to fall outside different SDs
  sds <- c(3, 4, 5, 10)
  two_tailed_prob <- 2 * pnorm(sds, lower.tail = FALSE)
  one_in_x <- 1 / two_tailed_prob

  tibble(
    SD = sds,
    probability = two_tailed_prob,
    percentage = scales::percent(two_tailed_prob, accuracy = 0.0001),
    one_in_x = one_in_x
  )
  # ----

  outliers |>
    count(species, name = "num_outliers")

  ggplot() +
    aes(x = fork_length, y = round_weight) +
    geom_point(data = lwd0 |> anti_join(outliers), colour = "grey80") +
    geom_point(data = outliers, colour = "red") +
    scale_y_continuous(trans = "log10") +
    scale_x_continuous(trans = "log10") +
    facet_grid(sex ~ species, scales = "free") +
    ggsidekick::theme_sleek()

  lwd <- lwd0 |> anti_join(outliers)
  saveRDS(lwd, here::here("data-generated", "length-weight-dat-clean.rds"))
} else {
  lwd <- readRDS(here::here("data-generated", "length-weight-dat-clean.rds"))
}

# Fit length-weight relationships to data used in analysis
# ------------------------
# Match `gfplot::fit_length_weight()` input format (i.e., `gfdata::get_survey_samples()`)
lwd2 <- lwd |>
    mutate(
      length = fork_length / 10,
      weight = round_weight,
      sex = ifelse(sex == "female", 2, 1))

if (!file.exists(here::here("data-generated", "lw-fits.rds"))) {

  get_lw_outputs <- function(dat, species, sex, output = "pars") {
    if (is.null(species)) {
      species <- unique(dat$species)
    } else {
      dat <- dat |> filter(species == !!species)
    }

    fit <- gfplot::fit_length_weight(dat,
      sex = sex, method = "tmb",
      too_high_quantile = 1,
      usability_codes = NULL
    )

    fit$predictions <- as_tibble(fit$predictions) |> mutate(species = species, sex = sex)
    fit$pars <- as_tibble(fit$pars) |> mutate(species = species, sex = sex)

    switch(output,
      predictions = fit$predictions,
      pars = fit$pars,
      no_tmb = fit[1:3],
      all = fit,
      stop("Invalid output type. Choose one of: 'predictions', 'pars', or 'all'.")
    )
  }

  lw_fits <- distinct(lwd2, species, sex) |>
    purrr::pmap(\(species, sex) {
      get_lw_outputs(dat = lwd2, species, sex, output = "no_tmb")
    })

  saveRDS(lw_fits, here::here("data-generated", "lw-fits.rds"))
} else {
  lw_fits <- readRDS(here::here("data-generated", "lw-fits.rds"))
}

# ------------------------------------------------------------

# Look at length-weight fits
# lw_fits <- readRDS(here::here("data-generated", "lw-fits.rds"))

# lw_preds <- purrr::map_dfr(lw_fits, "predictions")
# lw_pars <- purrr::map_dfr(lw_fits, "pars") |>
#   mutate(plot_text = glue::glue(
#       "\nln(a) = {round(log_a, 2)}\n b = {round(b, 2)}"
#     ))
# lw_dat <- purrr::map_dfr(lw_fits, "data")

# # Raw data plotted on fitted lw estimates
# ggplot(data = lw_preds) +
#   aes(x = length, y = weight, colour = factor(sex)) +
#   facet_wrap(sex ~ species, scales = "free") +
#   geom_line(aes(linetype = factor(sex))) +
#   geom_point(data = lw_dat, shape = 21, alpha = 0.5) +
#   geom_point(data = lw_dat |> filter(sarc_presence == 1), fill = "red", shape = 21, alpha = 0.8) +
#   scale_colour_manual(values = c("2" = "black", "1" = "grey50")) +
#   scale_linetype_manual(values = c("2" = "solid", "1" = "dashed")) +
#   geom_text(data = lw_pars,
#     aes(x = -Inf, y = Inf, label = paste0("ln(a) = ", round(log_a, 2), "\n  b = ", round(b, 2))),
#     hjust = -0.4, vjust = 2, size = 3)


# This was using the gfsynopsis lw pars
# ------------------------------------------------------------
# l_spp <- unique(lwd$species)
# l_spp_hyphen <- gsub(" |/", "-", l_spp)

# f <- list.files(here::here("data-raw", "length-weight-params"), full.names = TRUE)
# f <- f[grepl(paste0(l_spp_hyphen, collapse = "|"), basename(f))]

# lw_par <- purrr::map_dfr(f, readRDS) |>
#   mutate(species_common_name = gsub(" rockfish complex", "", species_common_name))

# ------------------------------------------------------------
# Body condition analysis
# ------------------------
lw_fits <- readRDS(here::here("data-generated", "lw-fits.rds"))
lw_preds <- purrr::map_dfr(lw_fits, "predictions")
lw_pars <- purrr::map_dfr(lw_fits, "pars") |>
  mutate(plot_text = glue::glue(
      "\nln(a) = {round(log_a, 2)}\n b = {round(b, 2)}"
    ))

cdat <-
  left_join(lwd2, y = lw_pars, by = c("species", "sex")) |>
  mutate(predicted_weight = exp(log_a + b * log(length))) |>
  mutate(condition = weight / predicted_weight,
         log_condition = log(condition),
         sarc_presence = factor(sarc_presence),
         sex = ifelse(sex == 2, "female", "male"),
         sex = factor(sex),
         species = factor(species),
         weight_kg = weight / 1000)

ggplot(data = cdat, mapping = aes(x = weight_kg, y = predicted_weight, color = factor(sarc_presence))) +
  geom_point(data = cdat |> filter(sarc_presence == 0), shape = 21, alpha = 0.5) +
  geom_point(data = cdat |> filter(sarc_presence == 1), shape = 21, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~ species + sex, scales = "free", labeller = labeller(.default = str_to_title), nrow = 4) +
  scale_color_manual(values = c("black", "red"))

cdat_modified <- cdat |>
  mutate(
    point_color = case_when(
      sarc_presence == 1          ~ "Infected",
      sarc_presence == 0 & sex == 'female' ~ "Female",
      sarc_presence == 0 & sex == 'male'   ~ "Male",
      TRUE ~ NA_character_
    ),
    point_color = factor(point_color, levels = c("Female", "Male", "Infected"))
  )

ggplot(data = cdat_modified, aes(x = weight_kg, y = predicted_weight, color = point_color)) +
  geom_point(shape = 21, alpha = 0.7) +
  geom_point(data = cdat_modified |> filter(sarc_presence == 1), shape = 21, alpha = 0.8, fill = "red") +
  geom_abline(intercept = 0, slope = 1) +
  ggh4x::facet_nested_wrap(
    ~ species + sex,
    scales = "free",
    labeller = labeller(.default = str_to_title),
    ncol = 2,
    strip = ggh4x::strip_nested(text_x =
        list(element_text(margin = margin(t = 5.5, b = -1, unit = "pt")),
          element_blank()
        ), by_layer_x = TRUE
    )
  ) +
  scale_color_manual(
    name = "Status", # A more informative legend title
    values = c("Female" = "black", "Male" = "grey40", "Infected" = "red")
  ) +
  ggsidekick::theme_sleek() +
  theme(ggh4x.facet.nestline = element_blank(),
        panel.spacing.y = unit(0, "pt"))

glimpse(cdat)

ggplot(data = cdat, aes(x = log_condition)) +
  geom_density(aes(fill = factor(sarc_presence), color = factor(sarc_presence)), alpha = 0.5) +
  facet_wrap(~ species + sex, scales = "free", labeller = labeller(.default = str_to_title), nrow = 4) +
  scale_fill_manual(values = c("black", "red")) +
  scale_color_manual(values = c("black", "red"))

# fit <- brm(
#   log(condition) ~ 0 + Intercept + sarc_presence * sex +
#     (1 + sarc_presence * sex | species),
#   data = cdat,
#   family = gaussian(),
#   iter = 2000L,
#   warmup = 1000L,
#   chains = 4L,
#   cores = 4L,
#   backend = "cmdstanr",
#   file = "cache/body-condition-fit1-sarc-lw-pars",
#   seed = 293829,
#   prior =
#     prior(normal(0, 2), class = b) +
#     prior(student_t(3, 0, 2), class = sd) +
#     prior(normal(0, 10), class = b, coef = Intercept),
#   control = list(max_treedepth = 12, adapt_delta = 0.8)
# )
# stop()
fit_f <- brm(
  log(condition) ~ 0 + Intercept + sarc_presence +
    (1 + sarc_presence | species),
  data = cdat |> filter(sex == "female"),
  family = gaussian(),
  iter = 2000L,
  warmup = 1000L,
  chains = 4L,
  cores = 4L,
  backend = "cmdstanr",
  file = "cache/body-condition-fit-f-sarc-lw-pars",
  seed = 293829,
  prior =
    prior(normal(0, 2), class = b) +
    prior(student_t(3, 0, 2), class = sd) +
    prior(normal(0, 10), class = b, coef = Intercept),
  control = list(max_treedepth = 12, adapt_delta = 0.98)
)

fit_m <- update(fit_f, newdata = cdat |> filter(sex == "male"), file = "cache/body-condition-fit-m-sarc-lw-pars")

summary(fit_f)
summary(fit_m)


# fit_f <- readRDS(here::here("cache", "body-condition-fit-f.rds"))
# fit_m <- readRDS(here::here("cache", "body-condition-fit-m.rds"))

mcmc_plot(
  fit_f,
  type = "areas",
  regex_pars = "^b_",
  prob = 0.95 # Shaded area is the 95% CI
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")

# Fixed effects
post_f <- fit_f |>
  spread_draws(b_Intercept, b_sarc_presence1) |>
  mutate(sex = "female")
post_m <- fit_m |>
  spread_draws(b_Intercept, b_sarc_presence1) |>
  mutate(sex = "male")
post_fe <- bind_rows(post_f, post_m) |>
  pivot_longer(
    cols = c(b_Intercept, b_sarc_presence1),
    names_to = "term",
    values_to = "fe_coef"
  ) |>
  mutate(term = gsub("^b_", "", term))

# Species effects
# post_re_f <-
# post_re_m <- fit_m |> spread_draws(r_species[species, term]) |>
#   mutate(sex = "male")
# post <- bind_rows(post_re_f, post_re_m)


# population
post_fe <- bind_rows(
  fit_f |> spread_draws(b_Intercept, b_sarc_presence1) |> mutate(sex = "female"),
  fit_m |> spread_draws(b_Intercept, b_sarc_presence1) |> mutate(sex = "male")
  ) |>
  pivot_longer(
    cols = c(b_Intercept, b_sarc_presence1),
    names_to = "term",
    values_to = "fe_coef"
  ) |>
  mutate(term = gsub("^b_", "", term))

# species effects
post_spp <- bind_rows(
  fit_f |> spread_draws(r_species[species, term]) |> mutate(sex = "female"),
  fit_m |> spread_draws(r_species[species, term]) |> mutate(sex = "male")
) |>
  ungroup() |>
  left_join(post_fe) |>
  mutate(combined = r_species + fe_coef,
         species = gsub("\\.|/", " ", species))

# combine
post_all <- bind_rows(
  post_fe |> mutate(species = "population", combined = fe_coef),
  post_spp
) |>
  mutate(term = gsub("sarc_presence1", "Infection", term),
         term = factor(term, levels = c("Intercept", "Infection"))
  )

# Coefficient posteriors
post_all |>
  mutate(species = forcats::fct_rev(species),
         species = forcats::fct_relabel(species, str_to_title),
         sex = forcats::fct_relabel(sex, str_to_title)) |>
  mutate(species = forcats::fct_relevel(species, "Population")) |>
  ggplot() +
  aes(x = combined, y = species, colour = sex, fill = sex) +
  geom_rect(aes(ymin = -Inf, ymax = 1.5, xmin = -Inf, xmax = Inf),
                      fill = "grey92", color = NA, inherit.aes = FALSE) +
  facet_wrap(~ term, scale = "free_x") +
  geom_vline(data = distinct(post_all, term) |> mutate(combined = ifelse(term == "Intercept", NA, 0)),
    aes(xintercept = combined), colour = "grey50", linetype = "dotted") +
  ggdist::stat_pointinterval(.width = c(0.5, 0.95),
    size = 2, linewidth = 1,
    position = ggstance::position_dodgev(height = -0.3)) +
  scale_color_manual(values = c("Female" = "black", "Male" = "grey60")) +
  scale_fill_manual(values = c("Female" = "black", "Male" = "grey60")) +
  guides(color = guide_legend(title = "Sex"), fill = guide_legend(title = "Sex")) +
  ggsidekick::theme_sleek(base_size = 12) +
  xlab("Coefficient estimate") +
  theme(axis.title.y.left = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 11))


# Does length affect body conditions?
# ------------------------
# Since condition is effectively a function of length, look at residuals

lw_fit_f <- brms(
  log(observed_weight) ~ 0 + Intercept + log(observed_length) +
                         (1 + log(observed_length) | species),
  data = cdat,
  family = gaussian(),
  iter = 2000L,
  warmup = 1000L,
  chains = 4L,
  cores = 4L,
  backend = "cmdstanr",
  file = "cache/lw-fit-f",
)

test <- cdat |> filter(sex == "female")
test$condition_residual <- residuals(fit_f)[, "Estimate"]



# Posterior predictive checks
# ------------------------
brms::pp_check(fit_f, ndraws = 50)

# Observations
y_obs_f <- cdat |> filter(sex == "female") |> pull(log_condition)
# Posterior predictions
yrep_f <- posterior_predict(fit_f)
ppc_stat(y = y_obs_f, yrep = yrep_f, stat = mean)
ppc_dens_overlay(y_obs_f, yrep_f)
ppc_stat_2d(y = y_obs_f, yrep = yrep_f, stat = c("mean", "sd"))

# Species means
ppc_stat_grouped(y = y_obs_f, yrep = yrep_f, group = fit_f$data$species, stat = mean)


# Prior checking
# ------------------------
sp <- filter(cdat, species == "redbanded rockfish", sex == "female")
sp <- filter(cdat, species == "yelloweye rockfish", sex == "female")
obs <- mutate(sp, .prediction = log(condition), .draw = 0)

pp <- tidybayes::predicted_draws(fit_f, newdata = sp, ndraws = 8)

ggplot(pp, aes(condition, .prediction)) +
  geom_point(aes(colour = sarc_presence), shape = 21, alpha = 0.5) +
  facet_wrap(~.draw) +
  scale_x_log10() +
  scale_color_manual(values = c("grey", "red")) +
  ggtitle("Posterior predictive simulation")


priors <- get_prior(fit_f)
fit_f_p <- update(fit_f, sample_prior = "only", prior = priors)
fit_m_p <- update(fit_m, sample_prior = "only", prior = priors)

mcmc_areas(as_draws_df(fit_f_p), regex_pars = c("^b_"))
mcmc_areas(as_draws_df(fit_f_p), regex_pars = c("^sd_")) +
  xlim(c(0, 25))

# Get prior draws
f_comb <- bind_rows(
  as_draws_df(fit_f_p) |> mutate(source = "prior"),
  as_draws_df(fit_f) |> mutate(source = "posterior")
)
f_comb |>
  pivot_longer(cols = starts_with("b_"), names_to = "parameter") |>
  ggplot(aes(x = value, fill = source)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~parameter, scales = "free")
f_comb |>
  pivot_longer(cols = starts_with("sd_"), names_to = "parameter") |>
  ggplot(aes(x = value, fill = source)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~parameter, scales = "free") +
  xlim(c(0, 20))

m_comb <- bind_rows(
  as_draws_df(fit_m_p) |> mutate(source = "prior"),
  as_draws_df(fit_m) |> mutate(source = "posterior")
)
m_comb |>
  pivot_longer(cols = starts_with("b_"), names_to = "parameter") |>
  ggplot(aes(x = value, fill = source)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~parameter, scales = "free")
m_comb |>
  pivot_longer(cols = starts_with("sd_"), names_to = "parameter") |>
  ggplot(aes(x = value, fill = source)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~parameter, scales = "free") +
  xlim(c(0, 20))

