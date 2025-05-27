library(tidyverse)
library(scales)
library(brms)
library(bayesplot)
library(tidybayes)
library(patchwork)
library(flextable)

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
# fit1 <- brm(
#   sarc_presence ~ 0 + Intercept + maturity_factor * sex +
#   (1 + maturity_factor * sex | species),
#   family = bernoulli(),
#   data = dat,
#   iter = 4000L,
#   warmup = 1000L,
#   chains = 4L,
#   cores = 4L,
#   backend = "cmdstanr",
#   prior = c(prior(normal(0, 5), class = b) +
#             prior(student_t(3, 0, 2), class = sd) +
#             prior(normal(0, 10), class = b, coef = Intercept)),
#   control = list(max_treedepth = 12, adapt_delta = 0.85)
# )
# saveRDS(fit1, file.path(fit_dir, "infection-by-maturity-bin-brms-by-species.rds"))
# beepr::beep()
fit1 <- readRDS(file.path(fit_dir, "infection-by-maturity-bin-brms-by-species.rds"))
# plot(fit1)

# Q2: Do immature fish have a higher number of sarcs than mature fish (# sarcs ~ immature vs mature)
# options(mc.cores = parallel::detectCores() - 2)
# fit2 <- brm(
#   sarc_count ~ 0 + Intercept + maturity_factor * sex +
#     (1 + maturity_factor * sex | species),
#   family = negbinomial(),
#   data = dat,
#   iter = 4000L,
#   warmup = 1000L,
#   chains = 4L,
#   cores = 4L,
#   backend = "cmdstanr",
#   prior = c(prior(normal(0, 5), class = b) +
#             prior(student_t(3, 0, 2), class = sd) +
#             prior(normal(0, 10), class = b, coef = Intercept)),
#   control = list(max_treedepth = 12, adapt_delta = 0.85)
# )
# beepr::beep()
# saveRDS(fit2, file.path(fit_dir, "sarc-count-by-immature-mature-brms.rds"))
fit2 <- readRDS(file.path(fit_dir, "sarc-count-by-immature-mature-brms.rds"))
# plot(fit2)

# Q3: Does infection rate differ across maturity bins (1 through 7)
# fit3 <- brm(
#   sarc_presence ~ 0 + Intercept + maturity_bin * sex,
#   family = bernoulli(),
#   data = dat,
#   iter = 4000L,
#   warmup = 1000L,
#   chains = 4L,
#   cores = 4L,
#   backend = "cmdstanr",
#   prior = c(prior(normal(0, 5), class = b) +
#             # prior(student_t(3, 0, 2), class = sd) +
#             prior(normal(0, 10), class = b, coef = Intercept)),
#   control = list(max_treedepth = 12, adapt_delta = 0.85)
# )
# beepr::beep()
# saveRDS(fit3, file.path(fit_dir, "infection-by-maturity-bin-brms.rds"))
fit3 <- readRDS(file.path(fit_dir, "infection-by-maturity-bin-brms.rds"))
# plot(fit3)

# With combined bins 4 and 5 for both male and female (1, 2, 3, 4/5, 6, 7)
# options(mc.cores = parallel::detectCores() - 2)
fit3b <- brm(
  sarc_presence ~ 0 + Intercept + maturity_bin2 * sex +
    (1 + maturity_bin2 * sex | species),
  family = bernoulli(),
  data = main_spp_dat,
  iter = 4000L,
  warmup = 1000L,
  chains = 4L,
  cores = 4L,
  backend = "cmdstanr",
  prior = c(prior(normal(0, 5), class = b) +
            prior(student_t(3, 0, 2), class = sd) +
            prior(normal(0, 10), class = b, coef = Intercept)),
  control = list(max_treedepth = 12, adapt_delta = 0.95)
)
beepr::beep()
saveRDS(fit3b, file.path(fit_dir, "infection-by-maturity-bin-collapsed4-5-brms.rds"))
fit3b <- readRDS(file.path(fit_dir, "infection-by-maturity-bin-collapsed4-5-brms.rds"))
# plot(fit3b)

# Prior checks


# Immature vs Mature
# ------------------------
#-- Q1 --#
# In natural space calculate the ratio of sarc infection in immature to mature:
nd <- expand.grid(
  species = unique(dat$species),
  sex = unique(dat$sex),
  maturity_factor = unique(dat$maturity_factor),
  maturity_bin = unique(dat$maturity_bin)
)

lpd1 <- fit1 |>
  add_linpred_draws(newdata = nd |> distinct(sex, maturity_factor), re_formula = NA, transform = TRUE)

lpd_wide1 <- lpd1 |>
  ungroup() |>
  select(-.row, -.chain, -.iteration) |>
  arrange(.draw, sex, maturity_factor) |>
  pivot_wider(names_from = c(sex, maturity_factor), values_from = .linpred, names_sep = "_")

p1 <- lpd_wide1 |>
  mcmc_intervals(pars = c("female_immature", "female_mature", "male_immature", "male_mature")) +
  scale_y_discrete(labels = c(
    "female_immature" = "female immature",
    "female_mature"   = "female mature",
    "male_immature"   = "male immature",
    "male_mature"     = "male mature"
  ))
imm_to_mat_ratio <- lpd_wide1 |>
  mutate(f_imm_mat_ratio = female_immature / female_mature,
         m_imm_mat_ratio = male_immature / male_mature)
saveRDS(imm_to_mat_ratio, here::here("data-generated", "overall-immature-to-mature-ratio.rds"))

# p2 <- imm_to_mat_ratio |>
# mcmc_intervals(pars = c("f_imm_mat_ratio", "m_imm_mat_ratio"), prob = 0.50, point_est = "median") +
#   geom_vline(xintercept = 1, linetype = "dashed") +
#   scale_x_continuous(limits = c(0, 5), breaks = 0:5) +
#   scale_y_discrete(labels = c(
#     "f_imm_mat_ratio" = "Female",
#     "m_imm_mat_ratio" = "Male"
#   )) +
#   xlab("Ratio of posterior infection probability: immature / mature")
# p2

imm_to_mat_ratio |>
  gather(-1, key = "par", value = "posterior") |>
  filter(par %in% c("f_imm_mat_ratio", "m_imm_mat_ratio")) |>
ggplot(aes(x = posterior, y = forcats::fct_rev(par))) +
  ggdist::stat_halfeye(.width = c(.5, .95)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_continuous(limits = c(0, 6)) +
  scale_y_discrete(labels = c(
    "f_imm_mat_ratio" = "Female",
    "m_imm_mat_ratio" = "Male"
  )) +
  xlab("Ratio of infection probability\n(immature / mature)") +
  theme(axis.title.y = element_blank())
ggsave(here::here("figures","maturity-factor-male-female-ratio.pdf"), width = 4.2, height = 3)
ggsave(here::here("figures","maturity-factor-male-female-ratio.png"), width = 4.2, height = 3)

# 35% probability that immature females have 2x higher infection rates than mature females
mean(imm_to_mat_ratio$f_imm_mat_ratio > 2)
mean(imm_to_mat_ratio$f_imm_mat_ratio > 1) # 92% probability of having more infections in immature
# 73% probability that immature males have 2x higher infection rates than mature males
mean(imm_to_mat_ratio$m_imm_mat_ratio > 2)

# Quantile summaries for immature:mature ratios
quantiles1 <- imm_to_mat_ratio |>
  summarise(across(
    # c(female_immature, female_mature, male_immature, male_mature),
    -.draw,
    list(
      median = ~ median(.),
      # q25 = ~ quantile(., 0.25),
      # q75 = ~ quantile(., 0.75),
      q2.5 = ~ quantile(., 0.025),
      q97.5 = ~ quantile(., 0.975)
    ),
    .names = "{.col}_{.fn}"
  )) |>
  pivot_longer(
    everything(),
    names_to = c("Group", "Quantile"),
    names_pattern = "^(.*)_(median|q25|q75|q2\\.5|q97\\.5)$"
  ) |>
  pivot_wider(names_from = Quantile, values_from = value) |>
  rename_with(~ gsub("_", " ", .x)) |>
  relocate(Group)

# ----------
# By species - comparing immature to mature
# ----------
lpd1 <- fit1 |>
  add_linpred_draws(newdata = nd |> distinct(species, sex, maturity_factor),
    re_formula = NULL, transform = TRUE)

lpd1 |>
  filter(species %in% main_spp$species) |>
  mutate(maturity_factor = forcats::fct_relevel(maturity_factor, "mature", "immature"),
    species = str_to_title(species)) |>
  ggplot(aes(x = .linpred, y = maturity_factor)) +
  stat_halfeye(.width = c(0.025, 0.975)) +
  labs(x = "Predicted probability", y = NULL) +
  # facet_grid(rows = vars(species), col = vars(sex), switch = "y", scales = "free_y") +
  ggh4x::facet_nested_wrap(. ~ species + sex, ncol = 2,
    labeller = labeller(.default = function(x) str_to_title(x))) +
  ggsidekick::theme_sleek(base_size = 10) +
  theme(
    panel.spacing.y = unit(0.1, "lines"),
    strip.placement = "outside",
    # strip.text.y.left = element_text(angle = 0, hjust = 0),
    strip.background = element_blank()
  ) +
  xlim(c(0, 0.21))
ggsave(here::here("figures", "maturity-factor-by-species.png"), width = 5.5, height = 5.4)

lpd_wide1 <- lpd1 |>
  ungroup() |>
  select(-.row, -.chain, -.iteration) |>
  # arrange(.draw, sex, maturity_factor) |>
  pivot_wider(names_from = c(sex, maturity_factor), values_from = .linpred, names_sep = "_")  |>
  mutate(f_imm_mat_ratio = female_immature / female_mature,
         m_imm_mat_ratio = male_immature / male_mature)

lpd_ratio_long1 <- lpd_wide1 |>
    pivot_longer(cols = -c(species, .draw),names_to = "group", values_to = "value")

ratio_order <- lpd_ratio_long1 |>
  filter(group == "f_imm_mat_ratio") |>
  group_by(species, group) |>
  reframe(mean_ratio = mean(value)) |>
  arrange(-mean_ratio) |>
  pull(species) |>
  as.character()


lpd_ratio_long1 |>
  filter(species %in% main_spp$species) |>
  filter(group %in% c("f_imm_mat_ratio", "m_imm_mat_ratio")) |>
  mutate(species = forcats::fct_relevel(species, ratio_order)) |>
  ggplot(aes(x = value, y = forcats::fct_rev(group))) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    # ggdist::stat_halfeye(.width = c(.5, .95)) +
    stat_pointinterval(.width = c(0.025, 0.5, 0.975), aes(colour = group)) +
  theme(panel.spacing.y = unit(0.1, "lines"),
    strip.placement = "outside",
    strip.text.y = element_text(angle = 0, hjust = 0),
    strip.text.y.left = element_text(angle = 0, hjust = 0),
    strip.background = element_blank()
  ) +
  labs(x = "Ratio", y = NULL) +
  scale_colour_manual(values = c("f_imm_mat_ratio" = "black", "m_imm_mat_ratio" = "grey"),
    labels = c("f_imm_mat_ratio" = "Female", "m_imm_mat_ratio" = "Male")) +
  guides() +
  scale_y_discrete(labels = c(
    "f_imm_mat_ratio" = "Female",
    "m_imm_mat_ratio" = "Male"
  )) +
  facet_grid(row = vars(species), switch = "y",
    labeller = labeller(.default = function(x) str_to_title(x))) +
  xlab("Ratio of infection probability\n(immature / mature)") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "top")#

ggsave(here::here("figures", "species-immature-to-mature-ratio.pdf"), width = 4, height = 5)
ggsave(here::here("figures", "species-immature-to-mature-ratio.png"), width = 4, height = 5)

# Quantile summaries for immature:mature ratios
quantiles1 <- lpd_wide1 |>
  group_by(species) |>
  summarise(across(
    c(f_imm_mat_ratio, m_imm_mat_ratio),
    list(
      median = ~ median(.),
      q2.5 = ~ quantile(., 0.025),
      q97.5 = ~ quantile(., 0.975)
    ),
    .names = "{.col}_{.fn}"
  )) |>
  pivot_longer(
    -species,
    names_to = c("Group", "Quantile"),
    names_pattern = "^(.*)_(median|q25|q75|q2\\.5|q97\\.5)$"
  ) |>
  pivot_wider(names_from = Quantile, values_from = value) |>
  rename_with(~ gsub("_", " ", .x)) |>
  relocate(Group)
quantiles1

# @Question: Similar text values like for the overall, but for each species? certain species?
# # 96% probability that immature females have 1.75x higher infection rates than mature females
# mean(imm_to_mat_ratio$f_imm_mat_ratio > 1.75)
# # 98% probability that immature males have 3x higher infection rates than mature males
# mean(imm_to_mat_ratio$m_imm_mat_ratio > 3)


# ----
# Q3: Does infection rate differ across maturity bins (collapsed 4 and 5)
# In natural space calculate the ratio of sarc infection in immature to mature:
nd <- expand.grid(
  species = unique(dat$species),
  sex = unique(dat$sex),
  maturity_factor = unique(dat$maturity_factor),
  maturity_bin2 = unique(dat$maturity_bin2)
)

lpd3 <- fit3b |>
  add_linpred_draws(newdata = nd |> distinct(sex, maturity_bin2), re_formula = NA, transform = TRUE)

lpd3 |>
  mutate(sex = as.character(sex), maturity_code = as.numeric(as.character(maturity_bin2))) |>
  left_join(mat_lu) |>
  mutate(maturity_desc = factor(maturity_desc))

lpd_wide3 <- lpd3 |>
  ungroup() |>
  select(-.row, -.chain, -.iteration) |>
  arrange(.draw, sex, maturity_bin2) |>
  pivot_wider(names_from = c(sex, maturity_bin2), values_from = .linpred, names_sep = "_")

# This gets at Table 3:
f_ratios <- lpd_wide3 |>
  select(starts_with("female")) |>
  mutate(`1-2` = female_1 / female_2,
         `1-3` = female_1 / female_3,
         `1-4.5` = female_1 / female_4.5,
         `1-6` = female_1 / female_6,
         `1-7` = female_1 / female_7) |>
  mutate(`2-3` = female_2 / female_3,
         `2-4.5` = female_2 / female_4.5,
         `2-6` = female_2 / female_6,
         `2-7` = female_2 / female_7) |>
  select(!starts_with("female")) |>
  mutate(sex = "female")

m_ratios <- lpd_wide3 |>
  select(starts_with("male")) |>
  mutate(`1-2` = male_1 / male_2,
         `1-3` = male_1 / male_3,
         `1-4.5` = male_1 / male_4.5,
         `1-6` = male_1 / male_6,
         `1-7` = male_1 / male_7) |>
  mutate(`2-3` = male_2 / male_3,
         `2-4.5` = male_2 / male_4.5,
        #  `2-5` = male_2 / male_5,
         `2-6` = male_2 / male_6,
         `2-7` = male_2 / male_7) |>
  select(!starts_with("male")) |>
  mutate(sex = "male")

bin_ratios <- bind_rows(f_ratios, m_ratios) |>
  pivot_longer(cols = -sex, names_to = "bins", values_to = "ratio") |>
  separate(bins, into = c("group1", "group2"), sep = "-") |>
  mutate(group1 = as.numeric(group1), group2 = as.numeric(group2)) |>
  left_join(mat_lu |> distinct(sex, maturity_code2, maturity_desc2), by = c("sex", "group2" = "maturity_code2")) |>
  mutate(maturity_desc = forcats::fct_inorder(maturity_desc2))

bin_ratios |>
  group_by(sex, group1, maturity_desc) |>
  summarise(`q2.5`  = quantile(ratio, 0.025),
            `q50`   = quantile(ratio, 0.50),
            `q97.5` = quantile(ratio, 0.975)
  ) |>
  mutate(
    Sex = str_to_title(sex),
    `Group 1` = ifelse(group1 == 1, "Immature", "Maturing"),
    `Group 2` = forcats::fct_relabel(maturity_desc, stringr::str_to_title),
    `Median ratio (95% CI)` = paste0(
      round(q50, 2), " (", round(q2.5, 2), " to ", round(q97.5, 2), ")"
  )) |>
  ungroup() |>
  select(Sex, `Group 1`, `Group 2`, `Median ratio (95% CI)`) |>
  group_by(Sex) |>
  mutate(
    Sex = if_else(duplicated(Sex), "", Sex),
    `Group 1` = if_else(duplicated(`Group 1`), "", `Group 1`)
  ) |>
  flextable() |>
  theme_booktabs() |>
  autofit()


table(dat$maturity_code2, dat$sarc_presence, dat$sex)

# Maybe not the most intuitive way of showing this, but it is consistent with the other ways
p1 <- bin_ratios |>
  filter(sex == "female") |>
  ggplot(aes(x = ratio, y = forcats::fct_rev(maturity_desc2))) +
  stat_pointinterval() +
  geom_vline(xintercept = 1, linetype = "dashed") +
  facet_grid(row = vars(group1), col = vars(sex)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0, 6), oob = scales::oob_keep) +
  labs(y = "Maturity status") +
  theme(strip.text.y = element_text(angle = 0))
p1
p2 <- bin_ratios |>
  filter(sex == "male") |>
  ggplot(aes(x = ratio, y = forcats::fct_rev(maturity_desc2))) +
  stat_pointinterval() +
  geom_vline(xintercept = 1, linetype = "dashed") +
  facet_grid(row = vars(group1), col = vars(sex)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0, 15), oob = scales::oob_keep) +
  theme(axis.title.y = element_blank(),
        strip.text.y = element_text(angle = 0))
  # labs(y = "Maturity status")
p1 + p2

# ------------------------------------------------------------------------------
# Q2: Do immature fish have a higher number of sarcs than mature fish (original Figure 6?)
# - yes but since this isn't splitting probability of encounter from count, it's just saying the same as the presence model
nd <- expand.grid(
  species = unique(dat$species),
  sex = unique(dat$sex),
  maturity_factor = unique(dat$maturity_factor)
)

lpd2 <- fit2 |>
  add_linpred_draws(newdata = nd |> distinct(sex, maturity_factor), re_formula = NA, transform = TRUE)

ggplot(lpd2, aes(x = .linpred, y = forcats::fct_rev(maturity_factor))) +
  stat_halfeye(.width = c(0.5, 0.8, 0.95)) +
  labs(x = "Expected number of cysts", y = "Maturity") +
  facet_wrap(~ sex)

lpd_wide2 <- lpd2 |>
  ungroup() |>
  select(-.row, -.chain, -.iteration) |>
  arrange(.draw, sex, maturity_factor) |>
  pivot_wider(names_from = c(sex, maturity_factor), values_from = .linpred, names_sep = "_")

count_ratio <- lpd_wide2 |>
  mutate(female = female_immature / female_mature,
         male = male_immature / male_mature)

# Ratio of expected NUMBER of sarc cysts
count_q <- count_ratio |>
  summarise(across(
    c(female, male),
    list(
      median = ~ median(.),
      q2.5 = ~ quantile(., 0.025),
      q97.5 = ~ quantile(., 0.975)
    ),
    .names = "{.col}_{.fn}"
  )) |>
  pivot_longer(cols = everything(), names_to = c("Sex", "Stat"), names_pattern = "(female|male)_(.*)") |>
  pivot_wider( names_from = Stat, values_from = value) |>
  mutate(Sex = str_to_title(Sex))
count_q

# -----
# better understand prior predictive checking
# -----
# fitf <- brm(
#   # sarc_count ~ 0 + Intercept + maturity_bin * sex +
#   #   (1 + maturity_bin * sex | species),
#   sarc_count ~ 0 + Intercept + maturity_bin +
#     (1 + maturity_bin | species),
#   family = zero_inflated_poisson(),
#   data = dat |> filter(sex == "female"),
#   iter = 2000L,
#   warmup = 500L,
#   chains = 4L,
#   cores = 4L,
#   backend = "cmdstanr",
#   prior = c(prior(normal(0, 5), class = b) +
#             prior(student_t(3, 0, 2), class = sd) +
#             prior(normal(0, 5), class = b, coef = Intercept) +
#             prior(beta(3, 1), class = zi)),  # the brms default is beta(1, 1)
#     # prior(normal(0, 5), class = b) +
#     # prior(student_t(3, 0, 2), class = sd) +
#     # prior(normal(0, 10), class = b, coef = Intercept),
#   control = list(max_treedepth = 12, adapt_delta = 0.95)
# )
# beepr::beep()
# saveRDS(fitf, file.path(fit_dir, "sarc-number-by-group-brms-zinp-female.rds"))
# fitf <- readRDS(file.path(fit_dir, "sarc-number-by-group-brms-zinp-female.rds"))


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
