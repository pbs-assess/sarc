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
# ld <- readRDS(here::here("data-generated", "length-dat.rds")) |>
#   filter(sex %in% c("female", "male"))

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

# -------------------
# Maturity: brms
# -------------------
ld <- readRDS("data-generated/length-dat.rds") |>
  filter(sex %in% c("female", "male"))
# fit_1f <- brm(
#   mature ~ 0 + Intercept + length_std * sarc_presence +
#     (1 + length_std * sarc_presence | species),
#   family = bernoulli(),
#   data = ld |> filter(sex == "female"),
#   iter = 1000L,
#   warmup = 500L,
#   chains = 4L,
#   cores = 4L,
#   backend = "cmdstanr",
#   prior =
#     prior(normal(0, 5), class = b) +
#     prior(student_t(3, 0, 2), class = sd) +
#     prior(normal(0, 10), class = b, coef = Intercept),
#   control = list(max_treedepth = 12, adapt_delta = 0.9)
# )
# dir.create("data-generated", showWarnings = FALSE)
# saveRDS(fit_1f, file.path(fit_dir, "maturity-length-stan-female.rds"))
# beepr::beep()
fit_1f <- readRDS(file.path(fit_dir, "maturity-length-stan-female.rds"))

# fit_1m <- update(fit_1f, newdata = ld |> filter(sex == "male"))
# saveRDS(fit_1m, file.path(fit_dir, "maturity-length-stan-male.rds"))
# beepr::beep()
fit_1m <- readRDS(file.path(fit_dir, "maturity-length-stan-male.rds"))


# Maturity ~ presence brms
# ------------------------
prior_summary(fit_1f)
get_variables(fit_1f)

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

pd_length_fe <- post_fe |>
  ggplot(aes(x = fe_coef)) +
  facet_grid(term ~ sex) +
  geom_vline(xintercept = 0) +
  ggsidekick::theme_sleek() +
  geom_density(fill = "grey90") +
  # facet_wrap(~ term, ncol = 1) +
  coord_cartesian(xlim = c(-3, 6), ylim = c(0, 1.7), expand = FALSE) +
  xlab("Coefficient estimate") + ylab("Posterior density") +
  ggtitle("Maturity ~ length: Females")
pd_length_fe
ggsave(here::here("figures", "maturity-length-posterior-densities.png"))

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


ggplot(post, aes(x = r_species, y = term)) +
  facet_wrap(~ sex) +
  geom_vline(xintercept = 0) +
  stat_halfeye()

post_spp <- left_join(post_fe, post) |>
  mutate(combined = r_species + fe_coef)

post_spp |>
  filter(term != "Intercept") |>
  group_by(term, sex, species) |>
  summarise(lwr = quantile(combined, probs = 0.0275), upr = quantile(combined, probs = 0.975), mid = median(combined)) |>
  ggplot(aes(x = mid, y = species, xmin = lwr, xmax = upr)) +
  facet_wrap(sex ~ term) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(pch = 21) +
  geom_linerange() + ggsidekick::theme_sleek() +
  xlab("Coefficient estimate") + theme(axis.title.y.left = element_blank())
ggsave(here::here("figures", "maturity-length-species-coef.png"))


get_variables(fit_1f) # get variable names
prior_summary(fit_1f) # just shows the list of priors used and what can be assigned?
# plot(fit_1f)
tidy(fit_1f)

nd <- expand.grid(species = unique(ld$species),
                  fork_length = seq(min(ld$fork_length), max(ld$fork_length), length.out = 200L),
                  sarc_presence = c(0, 1))
nd$length_std <- (nd$fork_length - mean(ld$fork_length)) / sd(ld$fork_length)
nd_f <- nd |> mutate(sex = "female")
nd_m <- nd |> mutate(sex = "male")

# Mean relationship
p_f_mean <- brms::posterior_linpred(fit_1f, newdata = nd, re_formula = NA)
nd_f_mean <- nd_f
nd_f_mean$lwr <- apply(p_f_mean, 2, quantile, probs = 0.05)
nd_f_mean$est <- apply(p_f_mean, 2, median)
nd_f_mean$upr <- apply(p_f_mean, 2, quantile, probs = 0.95)

p_m_mean <- brms::posterior_linpred(fit_1m, newdata = nd, re_formula = NA)
nd_m_mean <- nd_m
nd_m_mean$lwr <- apply(p_m_mean, 2, quantile, probs = 0.05)
nd_m_mean$est <- apply(p_m_mean, 2, median)
nd_m_mean$upr <- apply(p_m_mean, 2, quantile, probs = 0.95)

out <- bind_rows(nd_f_mean, nd_m_mean) |> as_tibble()

out |>
  mutate(sarc_pres_label = factor(sarc_presence, levels = c(0, 1), labels = c("No", "Yes"))) |>
ggplot(aes(x = fork_length, y = plogis(est),
           colour = sarc_pres_label, fill = sarc_pres_label)) +
  facet_wrap(~ sex) +
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
# Species-level
p_f <- brms::posterior_linpred(fit_1f, newdata = nd, re_formula = NULL)
nd_f_spp <- nd_f
nd_f_spp$lwr <- apply(p_f, 2, quantile, probs = 0.05)
nd_f_spp$est <- apply(p_f, 2, median)
nd_f_spp$upr <- apply(p_f, 2, quantile, probs = 0.95)

p_m <- brms::posterior_linpred(fit_1m, newdata = nd, re_formula = NULL)
nd_m_spp <- nd_m
nd_m_spp$lwr <- apply(p_m, 2, quantile, probs = 0.05)
nd_m_spp$est <- apply(p_m, 2, median)
nd_m_spp$upr <- apply(p_m, 2, quantile, probs = 0.95)

out_spp <- bind_rows(nd_f_spp, nd_m_spp) |> as_tibble()
out_spp |>
  mutate(sarc_pres_label = factor(sarc_presence, levels = c(0, 1), labels = c("No", "Yes"))) |>
ggplot(aes(x = fork_length, y = plogis(est),
           colour = sarc_pres_label, fill = sarc_pres_label)) +
  facet_grid(sex ~ stringr::str_to_title(species)) +
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
ggsave(here::here("figures", "maturity-length-infection-by-species.png"), width = 14, height = 5.5)

# --------------------------
# Compare expected length when probability of maturity > 0.5
get_p50 <- function(.int, .slope, .sd = sd(ld$fork_length), .mean = mean(ld$fork_length), p = 0.5) {
  xx <- -(log((1/p) - 1) + .int) / .slope
  (xx * .sd) + .mean
}

p50df <- group_by(post_spp, sex, species) |>
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

# p50df |> tidyr::pivot_longer(p50_sarc_0:p50_sarc_1) |>
#   group_by(species, sex, name) |>
#   summarise(lwr = quantile(value, probs = 0.0275), upr = quantile(value, probs = 0.975), mid = median(value)) |>
#   ggplot(aes(mid, y = species, colour = name, xmin = lwr, xmax = upr)) +
#   geom_pointrange(position = position_dodge(width = 0.5)) +
#   ggsidekick::theme_sleek() +
#   xlab("Length at p50 maturity") +
#   theme(legend.position = "bottom") +
#   facet_grid(. ~ sex) +
#   theme(panel.spacing.y = grid::unit(0, "mm"),
#         strip.background = element_blank(),
#         strip.text = element_text(angle = 0)) +
#   scale_colour_manual(values = c("p50_sarc_0" = "grey50", "p50_sarc_1" = "red"), labels = c("p50 No Infection", "p50 Infection"))

p50_diff <- p50df |>
  mutate(diff = p50_sarc_1 - p50_sarc_0) |>
  group_by(species, sex) |>
  summarise(
    lwr = quantile(diff, probs = 0.05),
    lwr2 = quantile(diff, probs = 0.25),
    upr = quantile(diff, probs = 0.95),
    upr2 = quantile(diff, probs = 0.75),
    mid = median(diff)
  )

p50_diff |>
  mutate(species = factor(species, levels = test |> filter(sex == "female") |> arrange(-mid) |> pull("species"))) |>
  ggplot(aes(mid, species, xmin = lwr, xmax = upr)) +
  geom_linerange(lwd = 0.4) +
  geom_linerange(aes(xmin = lwr2, xmax = upr2), lwd = .7) +
  geom_point(pch = 19) +
  geom_vline(xintercept = 0, lty = 2)+
  theme(axis.title.y.left = element_blank()) +
  ggsidekick::theme_sleek() +
  xlab("Difference in<br>length at 50% maturity<br>if infected with *Sarcotaces* sp.<br>(mm)") +
  theme(axis.title = ggtext::element_markdown()) +
  ylab("") +
  facet_grid(. ~ sex) +
  scale_x_continuous(limits = c(-100, 100), oob = scales::squish) # NOTE: copper rockfish cutoff
ggsave(here::here("figures", "maturity-length-p50-by-species.png"), width = 8, height = 7)

# ------------------------------------------------------------------------------
# Check how predictions compare to raw data
# I don't really get this
breaks <- seq(min(ld$fork_length, na.rm = TRUE), max(ld$fork_length, na.rm = TRUE), length.out = 30)

l_bins <- ld |>
  filter(sex == "female") |>
  mutate(fork_length_bin = cut(fork_length, breaks = breaks, include.lowest = TRUE)) |>
  group_by(species, fork_length_bin, sarc_presence) |>
  summarise(
    prop_mature = sum(mature) / n(),
    mean_fork_length = mean(fork_length, na.rm = TRUE),
    .groups = "drop"
  )

pp_f_population <- fit_1f |>
  add_predicted_draws(newdata = nd_f, re_formula = NULL)


breaks2 <- seq(min(ld$fork_length, na.rm = TRUE), max(ld$fork_length, na.rm = TRUE), length.out = 20)
test <- pp_f_population |>
  filter(species == "yelloweye rockfish") |>
  mutate(fork_length_bin = cut(fork_length, breaks = breaks2, include.lowest = TRUE)) |>
  group_by(species, .draw, fork_length_bin, sarc_presence) |>
  summarise(
    prop_mature = sum(.prediction) / n(),
    mean_fork_length = mean(fork_length, na.rm = TRUE),
    .groups = "drop"
  )
test |>
ggplot(data = _, aes(x = mean_fork_length, y = prop_mature, colour = factor(sarc_presence))) +
  geom_point()

test <- pp_f_population |>
  group_by(sex, species, sarc_presence, fork_length) |>
  summarise(median_p = quantile(.prediction, probs = 0.5),
            lwr = quantile(.prediction, probs = 0.05),
            upr = quantile(.prediction, probs = 0.95)) |>
ggplot(data = _, aes(x = fork_length, colour = factor(sarc_presence))) +
  geom_line(aes(y = median_p)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
  # geom_point(colour = "grey40") +
  facet_wrap(~ species)
