library(tidyverse)
library(scales)
library(brms)
library(bayesplot)
library(tidybayes)
library(patchwork)
library(flextable)

options(brms.file_refit = "on_change") # re-fit cached models if changes
dir.create("cache", showWarnings = FALSE)

theme_set(ggsidekick::theme_sleek())

fit_dir <- here::here("data-generated", "models")

e_spp <- readRDS(here::here("data-generated", "clean-data-encounter-summary.rds"))
spp_filter <- e_spp |> filter(n_encounters > 0) |> pull(species)
spp_n_levels <- e_spp |> arrange(-n) |> pull(species)

mat_lu <- readr::read_csv(here::here("data-raw", "maturity-lookup.csv")) |>
  mutate(specimen_sex_desc = toupper(sex)) |>
  mutate(maturity_code2 = ifelse(maturity_code %in% 4:5, 4.5, maturity_code),
      maturity_desc2 = case_when(
        maturity_code2 == 4.5 & sex == "female" ~ "fertilized/larvae",
        maturity_code2 == 4.5 & sex == "male"   ~ "large/running ripe",
        TRUE ~ maturity_desc
    )) |>
  arrange(maturity_code) |>
  mutate(maturity_desc = forcats::fct_inorder(maturity_desc),
         maturity_desc2 = forcats::fct_inorder(maturity_desc2))

dat <- readRDS(here::here("data-generated", "clean-data.rds")) |>
  filter(sex %in% c("female", "male")) |>
  filter(species %in% spp_filter) |> # only look at species with at least 1 sarc encounter
  filter(year >= 2019) |> # the number of sarcs was only counted in recent years
  filter(maturity_code %in% 1:7) |>
  left_join(mat_lu) |>
  mutate(maturity_bin = factor(maturity_code)) |>
  mutate(species = factor(species)) |>
  drop_na(sarc_count, sex, species, year, maturity_bin) |>
  mutate(maturity_factor = factor(ifelse(mature == 1, "mature", "immature"), levels = c("immature", "mature")),
         maturity_bin2 = factor(maturity_code2))
saveRDS(dat, here::here("data-generated", "clean-data-maturity-bins.rds"))
write_csv(dat, here::here("data-generated", "clean-data-maturity-bins.csv"))

main_spp <- dat |>
  count(species, maturity_factor, sarc_presence) |>
  pivot_wider(names_from = maturity_factor, values_from = n, values_fill = 0) |>
  group_by(species) |>
  filter(sarc_presence == 1 & immature > 10 & mature > 10)
  #  print(n = 24)

main_spp_dat <- filter(dat, species %in% main_spp$species)

# Immature vs mature
# ------------------------------------------------------------------------------
# Q1: Do immature fish have higher infection rates than mature fish (sarc presence) ~ immature vs mature)
# Infection rate differences by maturity factor / sex / species
options(mc.cores = parallel::detectCores() - 2)
fit1 <- brm(
  sarc_presence ~ 0 + Intercept + maturity_factor * sex +
  (1 + maturity_factor * sex | species),
  family = bernoulli(),
  data = dat,
  iter = 2000L,
  warmup = 1000L,
  chains = 4L,
  cores = 4L,
  backend = "cmdstanr",
  file = "cache/maturity-bin-fit1",
  prior = c(prior(normal(0, 2), class = b) +
            prior(student_t(3, 0, 2), class = sd) +
            prior(normal(0, 10), class = b, coef = Intercept) +
            prior(lkj(1), class = cor)
    ),
  seed = 762914,
  control = list(max_treedepth = 12, adapt_delta = 0.9)
)
# plot(fit1)

# Q2: Do immature fish have a higher number of sarcs than mature fish (# sarcs ~ immature vs mature)
# options(mc.cores = parallel::detectCores() - 2)
fit2 <- brm(
  sarc_count ~ 0 + Intercept + maturity_factor * sex +
    (1 + maturity_factor * sex | species),
  family = negbinomial(),
  data = dat,
  iter = 2000L,
  warmup = 1000L,
  chains = 4L,
  cores = 4L,
  file = "cache/maturity-bin-fit2",
  backend = "cmdstanr",
  prior = c(prior(normal(0, 2), class = b) +
            prior(student_t(3, 0, 2), class = sd) +
            prior(normal(0, 10), class = b, coef = Intercept)),
  seed = 31997429,
  control = list(max_treedepth = 12, adapt_delta = 0.9)
)
# plot(fit2)

# Q3: Does infection rate differ across maturity bins (1 through 7)
fit3 <- brm(
  sarc_presence ~ 0 + Intercept + maturity_bin * sex,
  family = bernoulli(),
  data = dat,
  iter = 2000L,
  warmup = 1000L,
  chains = 4L,
  cores = 4L,
  backend = "cmdstanr",
  file = "cache/maturity-bin-fit3",
  prior = c(prior(normal(0, 2), class = b) +
            # prior(student_t(3, 0, 2), class = sd) +
            prior(normal(0, 10), class = b, coef = Intercept)),
  seed = 1697017132,
  control = list(max_treedepth = 12, adapt_delta = 0.9)
)
# plot(fit3)

# With combined bins 4 and 5 for both male and female (1, 2, 3, 4/5, 6, 7)
# options(mc.cores = parallel::detectCores() - 2)
fit3b <- brm(
  sarc_presence ~ 0 + Intercept + maturity_bin2 * sex +
    (1 + maturity_bin2 * sex | species),
  family = bernoulli(),
  data = dat,
  iter = 2000L,
  warmup = 1000L,
  chains = 4L,
  cores = 4L,
  file = "cache/maturity-bin-fit3b",
  backend = "cmdstanr",
  seed = 9382919,
  prior = c(prior(normal(0, 2), class = b) +
            prior(student_t(3, 0, 2), class = sd) +
            prior(normal(0, 10), class = b, coef = Intercept)),
  control = list(max_treedepth = 12, adapt_delta = 0.99)
)
# plot(fit3b)

# New data for predictions
nd <- expand.grid(
  species = unique(dat$species),
  sex = unique(dat$sex),
  maturity_factor = unique(dat$maturity_factor),
  maturity_bin = unique(dat$maturity_bin)
)

# ------------------------
# Maturity factor (immature vs mature) infection probability
# ------------------------
pop_level1 <- fit1 |>
  add_epred_draws(newdata = nd |> distinct(sex, maturity_factor), re_formula = NA)

spp_level1 <- fit1 |>
  add_epred_draws(newdata = nd |> distinct(species, sex, maturity_factor),
    re_formula = NULL)

p_inf_factor <- spp_level1 |>
  bind_rows(pop_level1 |> mutate(species = "population"))

p_inf_factor_summary <- p_inf_factor |>
  group_by(species, sex, maturity_factor) |>
  summarise(mid = median(.epred),
            lwr = quantile(.epred, probs = 0.05),
            upr = quantile(.epred, probs = 0.95))

p_inf_spp_levels <- p_inf_factor_summary |>
  filter(sex == "female", maturity_factor == "immature") |>
  filter(species != "population") |>
  arrange(-mid) |>
  pull("species") %>%
  c("population", .)


prep_plot_data <- function(data, spp_levels, main_spp_only = TRUE) {
  dat <- data
  if (main_spp_only) {
    dat <- dat |>
      filter(tolower(species) %in% c(tolower(main_spp$species), "population"))
  }
  dat <- dat |>
    mutate(species = factor(species, levels = rev(spp_levels)),
           species = forcats::fct_relabel(species, str_to_title),
           sex = forcats::fct_relabel(sex, str_to_title))
}

p1 <- p_inf_factor |>
  prep_plot_data(spp_levels = p_inf_spp_levels, main_spp_only = TRUE) |>
  ggplot() +
  aes(x = .epred, y = fct_rev(maturity_factor), colour = sex) +
  geom_rect(data = tibble(species = factor("Population", levels = str_to_title(p_inf_spp_levels)), ymin = -Inf, ymax = Inf, xmin = -Inf, xmax = Inf),
            aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            fill = "grey92", colour = NA, inherit.aes = FALSE) +
  ggdist::stat_pointinterval(
      .width = c(0.5, 0.95),
      point_size = 2,
      position = position_dodge(width = -0.5)
  ) +
  ggh4x::facet_nested(species ~ ., nest_line = TRUE, switch = "y") +
  labs(colour = "Sex", y = "Species", x = "Probability of infection") +
  scale_colour_manual(values = c("Female" = "black", "Male" = "grey60")) +
  theme(
    legend.position = "top",
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    axis.title.y = element_blank(),
    panel.spacing = unit(-1, "mm"),
  )
p1

p1
ggsave(here::here("figures", "maturity-factor-P(infected)-main-spp.pdf"), width = 4.5, height = 5.3)
# ggsave(here::here("figures", "maturity-factor-P(infected)-main-spp.png"), width = 4.5, height = 5.3)
# this crazy function to update data :o
p1 %+% prep_plot_data(p_inf_factor, spp_levels = p_inf_spp_levels, main_spp_only = FALSE)
ggsave(here::here("figures", "maturity-factor-P(infected)-all-spp.pdf"), width = 4.5, height = 7.8)
# ggsave(here::here("figures", "maturity-factor-P(infected)-all-spp.png"), width = 4.5, height = 5.4)


# Immature vs Mature ratios
# ------------------------
pop_ratios <- pop_level1 |>
  ungroup() |>
  select(.draw, sex, maturity_factor, .epred) |>
  pivot_wider(names_from = maturity_factor, values_from = .epred) |>
  mutate(imm_mat_ratio = immature / (mature)) |>
  select(.draw, sex, imm_mat_ratio) |>
  mutate(species = "population")

spp_ratios <- spp_level1 |>
  ungroup() |>
  select(species, .draw, sex, maturity_factor, .epred) |>
  pivot_wider(names_from = maturity_factor, values_from = .epred) |>
  mutate(imm_mat_ratio = immature / (mature)) |>
  select(species, .draw, sex, imm_mat_ratio)

imm_mat_ratio_df <- bind_rows(spp_ratios, pop_ratios)

ratio_spp_levels <- imm_mat_ratio_df |>
  filter(sex == "female", species != "population") |>
  group_by(species) |>
  summarise(median_ratio = median(imm_mat_ratio)) |>
  arrange(-median_ratio) |>
  pull(species) %>%
  c("population", .)

p2 <- imm_mat_ratio_df |>
  prep_plot_data(spp_levels = ratio_spp_levels, main_spp_only = TRUE) |>
  ggplot() +
  aes(x = imm_mat_ratio, y = fct_rev(species), colour = sex) +
  geom_rect(data = tibble(species = factor("Population", levels = str_to_title(ratio_spp_levels)),
    ymin = -Inf, ymax = 1.5, xmin = -Inf, xmax = Inf),
            aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            fill = "grey92", colour = NA, inherit.aes = FALSE) +
  geom_vline(xintercept = 1, linetype = "dotted", colour = "grey50") +
  ggdist::stat_pointinterval(
      .width = c(0.5, 0.95),
      point_size = 2,
      position = position_dodge(width = -0.5)
  ) +
  labs(colour = "Sex", y = "Species", x = "Ratio (immature / mature)") +
  scale_colour_manual(values = c("Female" = "black", "Male" = "grey60")) +
  theme(
    legend.position = "top",
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    axis.title.y = element_blank()
  )

p2
ggsave(here::here("figures", "maturity-ratio-main-spp.pdf"), width = 4.5, height = 5.3)
# ggsave(here::here("figures", "maturity-ratio-main-spp.png"), width = 4.5, height = 5.3)
p2 %+% prep_plot_data(imm_mat_ratio_df, spp_levels = ratio_spp_levels, main_spp_only = FALSE)
ggsave(here::here("figures", "maturity-ratio-all-spp.pdf"), width = 4.5, height = 7.8)
# ggsave(here::here("figures", "maturity-ratio-all-spp.png"), width = 4.5, height = 5.3)

# Values for text
# 35% probability that immature females have 2x higher infection rates than mature females
mean(filter(imm_mat_ratio_df, sex == "female" & species == "population")$imm_mat_ratio > 1) # 92% probability of having more infections in immature
mean(filter(imm_mat_ratio_df, sex == "female" & species == "population")$imm_mat_ratio > 2)

# 73% probability that immature males have 2x higher infection rates than mature males
mean(filter(imm_mat_ratio_df, sex == "male" & species == "population")$imm_mat_ratio > 1)
mean(filter(imm_mat_ratio_df, sex == "male" & species == "population")$imm_mat_ratio > 2)

# ------------------------------------------------------------------------------
# Q3: Does infection rate differ across maturity bins (collapsed 4 and 5)
# ------------------------------------------------------------------------------
# In natural space calculate the ratio of sarc infection in immature to mature:
nd <- expand.grid(
  # species = unique(dat$species),
  species = main_spp$species,
  sex = unique(dat$sex),
  maturity_factor = unique(dat$maturity_factor),
  maturity_bin2 = unique(dat$maturity_bin2)
) |>
droplevels()

# Start with getting infection probability for each maturity bin
pop_level3 <- fit3b |>
  add_epred_draws(newdata = nd |> distinct(sex, maturity_bin2), re_formula = NA)

spp_level3 <- fit3b |>
  add_epred_draws(newdata = nd |> distinct(species, sex, maturity_bin2), re_formula = NULL)

p_inf_bins <- spp_level3 |>
  ungroup() |>
  bind_rows(pop_level3 |> mutate(species = "population"))

p_inf_bins_summary <- p_inf_bins |>
  group_by(species, sex, maturity_bin2) |>
  summarise(mid = median(.epred),
            lwr = quantile(.epred, probs = 0.05),
            upr = quantile(.epred, probs = 0.95))

p_inf_bins_spp_levels <- p_inf_bins_summary |>
  filter(sex == "female", maturity_bin2 == "1") |>
  filter(species != "population") |>
  arrange(-mid) |>
  pull("species") %>%
  c("population", .)

# option 1 with male and female on same plot
p_inf_bins |>
  mutate(species = factor(species, levels = p_inf_bins_spp_levels),
         sex = forcats::fct_relabel(sex, str_to_title)) |>
ggplot() +
  aes(x = maturity_bin2, y = .epred, colour = sex) +
  ggdist::stat_pointinterval(
      .width = c(0.5, 0.95),
      point_size = 2,
      position = position_dodge(width = 0.5)) +
  facet_wrap(. ~ species, ncol = 3, labeller = labeller(.default = str_to_title)) +
  # facet_grid(~ species) +
  scale_colour_manual(values = c("Female" = "black", "Male" = "grey60")) +
  scale_y_continuous(limits = c(0, 0.3), oob = scales::oob_keep) +
  labs(colour = "Sex", y = "Probability of infection", x = "Maturity status") +
  theme(legend.position = "top")
ggsave(here::here("figures", "maturity-bin-P(infected)-main-spp-FM-same-plot.pdf"), width = 7.1, height = 4.4)
# ggsave(here::here("figures", "maturity-bin-P(infected)-main-spp-FM-same-plot.png"), width = 7.1, height = 4.4)

# option 2 with male and female different plots
p1 <-
p_inf_bins |>
  filter(sex == "female") |>
  prep_plot_data(spp_levels = p_inf_bins_spp_levels, main_spp_only = TRUE) |>
ggplot() +
  aes(x = maturity_bin2, y = .epred, colour = sex) +
  ggdist::stat_pointinterval(
      .width = c(0.5, 0.95),
      point_size = 2,
      position = position_dodge(width = 0.5)) +
  # facet_wrap(. ~ species, ncol = 3) +
  facet_grid(species ~ ., scales = "free_y", labeller = labeller(.default = str_to_title)) +
  scale_colour_manual(values = c("Female" = "black", "Male" = "grey60")) +
  # scale_y_continuous(limits = c(0, 0.3), oob = scales::oob_keep) +
  labs(y = "Probability of infection", x = "Maturity status") +
  theme(legend.position = "top",
        legend.title = element_blank(),
        strip.text.y = element_blank())
p1

p2 <- p1 %+%
  prep_plot_data(p_inf_bins |> filter(sex == "male"),
    spp_levels = p_inf_bins_spp_levels,
    main_spp_only = TRUE) +
    theme(axis.title.y = element_blank(),
          strip.text.y = element_text(size = 10),
          strip.clip = "off")
((p1 + ggtitle("Female")) + (p2 + ggtitle("Male"))) +
  plot_layout(guides = "collect") & theme(legend.position = "none")
ggsave(here::here("figures", "maturity-bin-P(infected)-main-spp.pdf"), width = 7.1, height = 9.6)
# ggsave(here::here("figures", "maturity-bin-P(infected)-main-spp.png"), width = 7.1, height = 9.6)

pop_only1f <- p1 %+%
  prep_plot_data(p_inf_bins |> filter(sex == "female", species == "population"),
    spp_levels = p_inf_bins_spp_levels,
    main_spp_only = TRUE) +
    theme(axis.title.y = element_blank(),
          strip.text.y = element_blank(),
          strip.clip = "off")
pop_only1m <- p1 %+%
  prep_plot_data(p_inf_bins |> filter(sex == "male", species == "population"),
    spp_levels = p_inf_bins_spp_levels,
    main_spp_only = TRUE) +
    theme(axis.title.y = element_blank(),
          strip.text.y = element_blank(),
          strip.clip = "off")
pop_only1 <- ((pop_only1f + ggtitle("Female")) + (pop_only1m + ggtitle("Male"))) +
  plot_layout(guides = "collect") & theme(legend.position = "none")

# p1 <- p1 %+%
#   prep_plot_data(p_inf_bins |> filter(sex == "male"),
#     spp_levels = p_inf_bins_spp_levels,
#     main_spp_only = FALSE)
# p2 <- p1 %+%
#   prep_plot_data(p_inf_bins |> filter(sex == "male"),
#     spp_levels = p_inf_bins_spp_levels,
#     main_spp_only = FALSE)
# ((p1 + ggtitle("Female")) / (p2 + ggtitle("Male"))) +
#   plot_layout(guides = "collect")
# ggsave(here::here("figures", "maturity-bin-P(infected)-main-spp.pdf"), width = 6.7, height = 7.2)
# ggsave(here::here("figures", "maturity-bin-P(infected)-main-spp.png"), width = 6.7, height = 7.2)


# TODO: add all spp plot
# plot_posteriors(p_inf_bins, main_spp_only = FALSE, xtitle = "Probability of infection")
# ggsave(here::here("figures", "maturity-bin-P(infected)-all-spp.pdf"), width = 4.1, height = 5.4)
# ggsave(here::here("figures", "maturity-bin-P(infected)-all-spp.png"), width = 4.1, height = 5.4)


# Bin ratios
# ------------------------
pop_ratios3 <- pop_level3 |>
  ungroup() |>
  select(.draw, sex, maturity_bin2, .epred) |>
  pivot_wider(names_from = maturity_bin2, values_from = .epred) |>
  rename(`4-5` = `4.5`) |>
  mutate(`1/2` = `1` / `2`,
         `1/3` = `1`/ `3`,
         `1/4-5` = `1` / `4-5`,
         `1/6` = `1` / `6`,
         `1/7` = `1` / `7`) |>
  mutate(`2/3` = `2` / `3`,
         `2/4-5` = `2` / `4-5`,
         `2/6` = `2` / `6`,
         `2/7` = `2` / `7`) |>
  mutate(species = "population") |>
  select(species, sex, starts_with("1"), starts_with("2"))

spp_ratios3 <- spp_level3 |>
  ungroup() |>
  select(species, .draw, sex, maturity_bin2, .epred) |>
  pivot_wider(names_from = maturity_bin2, values_from = .epred) |>
  rename(`4-5` = `4.5`) |>
  mutate(`1/2` = `1` / `2`,
         `1/3` = `1` / `3`,
         `1/4-5` = `1` / `4-5`,
         `1/6` = `1` / `6`,
         `1/7` = `1` / `7`) |>
  mutate(`2/3` = `2` / `3`,
         `2/4-5` = `2` / `4-5`,
         `2/6` = `2` / `6`,
         `2/7` = `2` / `7`) |>
  select(species, .draw, sex, starts_with("1"), starts_with("2"))

bin_ratios <- bind_rows(spp_ratios3, pop_ratios3) |>
  select(-`1`, -`2`) |>
  pivot_longer(cols = c(starts_with("1"), starts_with("2")), names_to = "bins", values_to = "ratio")

bin_ratio_spp_levels <- p_inf_bins_spp_levels

# option 1 with male and female on same plot - cluttered
p_fm <- bin_ratios |>
  prep_plot_data(spp_levels = bin_ratio_spp_levels, main_spp_only = TRUE) |>
  filter(str_detect(bins, "^1")) |>
ggplot() +
  aes(x = bins, y = ratio, colour = sex) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  ggdist::stat_pointinterval(
      .width = c(0.5, 0.95),
      point_size = 2,
      position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c("Female" = "black", "Male" = "grey60")) +
  facet_grid(species ~ ., scales = "free_y", labeller = labeller(.default = str_to_title)) +
  # facet_wrap(. ~ species, ncol = 3, labeller = labeller(.default = str_to_title)) +
  # scale_y_continuous(limits = c(0, 15), oob = scales::oob_keep) +
  labs(colour = "Sex", y = "Ratio", x = "Maturity status comparison") +
  theme(legend.position = "top")

# option 2 with male and female on different plots
pf <- p_fm %+%
  prep_plot_data(bin_ratios |> filter(str_detect(bins, "1/"), sex == "female"),
    spp_levels = bin_ratio_spp_levels,
    main_spp_only = TRUE) +
    scale_y_continuous(limits = c(0, 10), oob = scales::oob_keep) +
    theme(strip.text.y = element_blank())
pm <- p_fm %+%
  prep_plot_data(bin_ratios |> filter(str_detect(bins, "1/"), sex == "male"),
    spp_levels = bin_ratio_spp_levels, main_spp_only = TRUE) +
    scale_y_continuous(limits = c(0, 30), oob = scales::oob_keep) +
    theme(axis.title.y = element_blank(),
          strip.text.y = element_text(size = 10),
          strip.clip = "off")

(pf + ggtitle("Female")) + (pm + ggtitle("Male")) + plot_layout(guides = "collect") & theme(legend.position = "none")
ggsave(here::here("figures", "maturity-bin-ratio-main-spp-immature.pdf"), width = 7.1, height = 9.6)
# ggsave(here::here("figures", "maturity-bin-ratio-main-spp-immature.png"), width = 7.1, height = 9.6)

# Population only - immature to other bins
pop_only2f <- p_fm %+%
  prep_plot_data(bin_ratios |> filter(str_detect(bins, "1/"), sex == "female", species == "population"),
    spp_levels = bin_ratio_spp_levels,
    main_spp_only = TRUE) +
    theme(axis.title.y = element_blank(),
          strip.text.y = element_blank(),
          strip.clip = "off")
pop_only2m <- p_fm %+%
  prep_plot_data(bin_ratios |> filter(str_detect(bins, "1/"), sex == "male", species == "population"),
    spp_levels = bin_ratio_spp_levels,
    main_spp_only = TRUE) +
    theme(axis.title.y = element_blank(),
          strip.text.y = element_blank(),
          strip.clip = "off")
pop_only2 <- ((pop_only2f + ggtitle("Female")) + (pop_only2m + ggtitle("Male"))) +
  plot_layout(guides = "collect") & theme(legend.position = "none")

# Maturity status 2 comparsion
# ------------------------
# p_fm2 <- p_fm %+% # very cluttered visually
#   prep_plot_data(bin_ratios |> filter(str_detect(bins, "^2")),
#     spp_levels = bin_ratio_spp_levels, main_spp_only = TRUE)

p_f2 <- p_fm %+%
  prep_plot_data(bin_ratios |> filter(str_detect(bins, "^2"), sex == "female"),
    spp_levels = bin_ratio_spp_levels, main_spp_only = TRUE) +
    scale_y_continuous(limits = c(0, 10), oob = scales::oob_keep) +
    theme(strip.text.y = element_blank())
p_m2 <- p_f2 %+%
  prep_plot_data(bin_ratios |> filter(str_detect(bins, "^2"), sex == "male"),
    spp_levels = bin_ratio_spp_levels, main_spp_only = TRUE) +
    scale_y_continuous(limits = c(0, 30), oob = scales::oob_keep) +
    theme(axis.title.y = element_blank(),
          strip.text.y = element_text(size = 10),
          strip.clip = "off")

(p_f2 + ggtitle("Female")) + (p_m2 + ggtitle("Male")) +
  plot_layout(guides = "collect") & theme(legend.position = "none")
ggsave(here::here("figures", "maturity-bin-ratio-main-spp-maturing.pdf"), width = 7.1, height = 9.6)
# ggsave(here::here("figures", "maturity-bin-ratio-main-spp-maturing.png"), width = 7.1, height = 9.6)

pop_only3f <- p_fm %+%
  prep_plot_data(bin_ratios |> filter(str_detect(bins, "^2"), sex == "female", species == "population"),
    spp_levels = bin_ratio_spp_levels,
    main_spp_only = TRUE) +
    theme(axis.title.y = element_blank(),
          strip.text.y = element_blank(),
          strip.clip = "off")
pop_only3m <- p_fm %+%
  prep_plot_data(bin_ratios |> filter(str_detect(bins, "^2"), sex == "male", species == "population"),
    spp_levels = bin_ratio_spp_levels,
    main_spp_only = TRUE) +
    theme(axis.title.y = element_blank(),
          strip.text.y = element_blank(),
          strip.clip = "off")
pop_only3 <- ((pop_only3f + ggtitle("Female")) + (pop_only3m + ggtitle("Male"))) +
  plot_layout(guides = "collect") & theme(legend.position = "none")

# Maturity bin population only trends
# ------------------------

(pop_only1 / pop_only2 / pop_only3) + plot_annotation(tag_levels = "a", tag_suffix = ")")
ggsave(here::here("figures", "maturity-bin-comparison-population.pdf"), width = 7.1, height = 7.3)

# ------------------------------------------------------------------------------
# Q2: Do immature fish have a higher number of sarcs than mature fish (original Figure 6?)
nd <- expand.grid(
  species = unique(dat$species),
  sex = unique(dat$sex),
  maturity_factor = unique(dat$maturity_factor)
) |>
  droplevels()

pop_level2 <- fit2 |>
  add_epred_draws(newdata = nd |> distinct(sex, maturity_factor), re_formula = NA)

spp_level2 <- fit2 |>
  add_epred_draws(newdata = nd |> distinct(species, sex, maturity_factor),
    re_formula = NULL)

p_count_factor <- spp_level2 |>
  bind_rows(pop_level2 |> mutate(species = "population")) |>
  ungroup()

p_count_factor_summary <- p_count_factor |>
  group_by(species, sex, maturity_factor) |>
  summarise(mid = median(.epred),
            lwr = quantile(.epred, probs = 0.05),
            upr = quantile(.epred, probs = 0.95))

p_count_spp_levels <- p_count_factor_summary |>
  filter(sex == "female", maturity_factor == "immature") |>
  filter(species != "population") |>
  arrange(-mid) |>
  pull("species") %>%
  c("population", .)

p1 <- p_count_factor |>
  prep_plot_data(spp_levels = p_count_spp_levels, main_spp_only = TRUE) |>
  ggplot() +
  aes(x = .epred, y = fct_rev(maturity_factor), colour = sex) +
  geom_rect(data = tibble(species = factor("Population", levels = str_to_title(p_inf_spp_levels)), ymin = -Inf, ymax = Inf, xmin = -Inf, xmax = Inf),
            aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            fill = "grey92", colour = NA, inherit.aes = FALSE) +
  ggdist::stat_pointinterval(
      .width = c(0.5, 0.95),
      point_size = 2,
      position = position_dodge(width = -0.5)
  ) +
  ggh4x::facet_nested(species ~ ., nest_line = TRUE, switch = "y") +
  labs(colour = "Sex", y = "Species", x = "Expected number of cysts") +
  scale_colour_manual(values = c("Female" = "black", "Male" = "grey60")) +
  theme(
    legend.position = "top",
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    panel.spacing = unit(-1, "mm"),
  )
p1
ggsave(here::here("figures", "maturity-factor-count-main-spp.pdf"), width = 4.8, height = 5.3)
# ggsave(here::here("figures", "maturity-factor-count-main-spp.png"), width = 4.8, height = 5.3)

p1 %+%
  prep_plot_data(p_count_factor, spp_levels = p_count_spp_levels, main_spp_only = FALSE)
ggsave(here::here("figures", "maturity-factor-count-all-spp.pdf"), width = 4.8, height = 7.8)
# ggsave(here::here("figures", "maturity-factor-count-all-spp.png"), width = 4.8, height = 7.8)

# Summary tables for text
# -----------------------
# Estimates from model fit1 (sarc_presence ~ maturity_factor * sex).
# quantiles for infection probabilities for immature and mature fish (by species and sex)
# see line ~300 in this script.
p_inf_factor_summary
saveRDS(p_inf_factor_summary, here::here("data-generated", "maturity-factor-infection-probabilities.rds"))

# Estimates from model fit3b (sarc_presence ~ maturity_bin2 * sex).
# quantiles for infection probabilities for different maturity bins (1, 2, 3, 4-5, 6, 7) (by species and sex)
# see line ~464 in this script.
p_inf_bins_summary
saveRDS(p_inf_bins_summary, here::here("data-generated", "maturity-bin-infection-probabilities.rds"))

# Estimates from model fit2 (sarc_count ~ maturity_factor * sex).
# quantiles for number of cysts for immature and mature fish (by species and sex)
# see line ~400 in this script.
p_count_factor_summary
saveRDS(p_count_factor_summary, here::here("data-generated", "maturity-factor-cyst-counts.rds"))

# Ratio dataframes
# ----------------
# Posterior draws of the ratio of infection probability in immature to mature fish.
# Used for reporting the probability that immature fish have higher infection rates than mature fish.
# see line ~278 in this script.
imm_mat_ratio_df
saveRDS(imm_mat_ratio_df, here::here("data-generated", "immature-to-mature-infection-ratio.rds"))

# Posterior draws of ratios of infection probabilities between different maturity bins (e.g., 1/2, 1/3, etc.), by species and sex.
# Useful for reporting infection risk across all maturity bins.
# see line ~464 in this script.
bin_ratios
saveRDS(bin_ratios, here::here("data-generated", "maturity-bin-infection-ratios.rds"))
