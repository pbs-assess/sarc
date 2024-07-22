# malcolm and christina
library(tidyverse)

# data <- read_csv("qry_sarc_maturity.csv")
full_DF <- read_csv("sarc_specimens.csv")
full_DF$year <- stringr::str_sub(full_DF$FE_BEGIN_RETRIEVAL_TIME,1,4)
# subset "full_DF" columns for ease
d <- full_DF[,c(1,5,6,8,11,12,15,18,31:33,35:37,41:44)]

# presence/ absence data in sarc_presence with values of Y, N, null -
# where null means not looked at,
# count data in sarc_count,
# sarc_comb combines sarc_presence and sarc_count: where na = not checked,
# 0 = checked and none, 1+ = checked and one or more, # = #

d <- d |>
  filter(sarc_comb != "na") |>
  mutate(sarc_presence = case_when(sarc_comb %in% c("0","N") ~ 0,
    sarc_comb %in% c("Y","1","2","3","4","5","6","7","8","9","10") ~ 1))

d$year <- as.numeric(d$year)
d$Fork_Length <- as.numeric(d$Fork_Length)
d$Round_Weight <- as.numeric(d$Round_Weight)
d$sarc_presence <- as.integer(d$sarc_presence)
names(d) <- tolower(names(d))
d$species_desc <- tolower(d$species_desc)
d$specimen_age <- as.numeric(d$specimen_age)
d$maturity_code <- as.integer(d$maturity_code)

d <- group_by(d, species_desc) |>
  mutate(n = n()) # |>
  # filter(n >= 100)

glimpse(d)

sort(table(d$maturity_code))

gfplot::maturity_assignment

d <- d |>
  filter(!is.na(maturity_code), maturity_code != 0L, maturity_code != 8L) |>
  # filter(!is.na(maturity_code)) |>
  mutate(mature = as.integer(maturity_code >= 3L))

tokeep <- d |> group_by(species_desc) |>
  summarise(
    n_pres_mat = sum(sarc_presence[mature]),
    n_abs_mat = sum(!sarc_presence[mature]),
    n_pres_imat = sum(sarc_presence[!mature]),
    n_abs_imat = sum(!sarc_presence[!mature])
  ) |>
  group_by(species_desc) |>
  filter(n_pres_mat > 0, n_abs_mat > 0, n_pres_imat > 0, n_abs_imat > 0) |>
  ungroup()

ggplot(d, aes(specimen_age, mature, colour = as.factor(sarc_presence))) +
  geom_jitter(width = 0, height = 0.1, alpha = 0.4) +
  facet_grid(species_desc~sarc_presence) +
  geom_smooth(method = "glm", method.args = list(family = binomial()), se = FALSE)

keep0 <- d |>
  filter(!is.na(fork_length), sarc_presence == 0) |>
  group_by(species_desc) |>
  summarise(total0 = sum(mature == 0) > 0, total1 = sum(mature == 1) > 0) |>
  filter(total0 & total1) |> select(species_desc)
keep1 <- d |>
  filter(!is.na(fork_length), sarc_presence == 1) |>
  group_by(species_desc) |>
  summarise(total0 = sum(mature == 0) > 0, total1 = sum(mature == 1) > 0) |>
  filter(total0 & total1) |> select(species_desc)
keep_length <- semi_join(keep0, keep1)

length_dat <- filter(d, species_desc %in% keep_length$species_desc, !is.na(fork_length))

ggplot(length_dat, aes(fork_length, mature, colour = as.factor(sarc_presence))) +
  geom_jitter(width = 0, height = 0.1, alpha = 0.9) +
  facet_grid(species_desc~sarc_presence) +
  geom_smooth(method = "glm", method.args = list(family = binomial()), se = FALSE)

library(glmmTMB)

length_dat$species_desc <- as.factor(length_dat$species_desc)

length_dat$length_std <- (length_dat$fork_length - mean(length_dat$fork_length)) / sd(length_dat$fork_length)

m1 <- glmmTMB(mature ~ length_std * sarc_presence + (1 + length_std + sarc_presence | species_desc), family = binomial(), data = length_dat)

summary(m1)

nd <- expand.grid(species_desc = unique(length_dat$species_desc), fork_length = seq(min(length_dat$fork_length), max(length_dat$fork_length), length.out = 100L), sarc_presence = c(0, 1))
nd$length_std <- (nd$fork_length - mean(length_dat$fork_length)) / sd(length_dat$fork_length)

nd$est <- predict(m1, newdata = nd)

ggplot(nd, aes(fork_length, plogis(est), colour = factor(sarc_presence))) + geom_line() +
  facet_wrap(~species_desc)

library(brms)

# pp <- brms::default_prior(
#   mature ~ 0 + Intercept + length_std * sarc_presence +
#     (1 + length_std * sarc_presence | species_desc),
#   family = bernoulli(),
#   data = length_dat,
#   iter = 500L,
#   warmup = 250L,
#   chains = 4L,
#   cores = 4L,
#   backend = "cmdstanr",
#   prior =
#     prior(normal(0, 2), class = b) +
#     prior(student_t(3, 0, 2), class = sd) +
#     prior(normal(0, 10), class = Intercept)
# )

fit1 <- brm(
  mature ~ 0 + Intercept + length_std * sarc_presence +
    (1 + length_std * sarc_presence | species_desc),
  family = bernoulli(),
  data = length_dat,
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
saveRDS(fit1, "data-generated/maturity-stan-model.rds")

prior_summary(fit1)
fit1

plot(fit1)

nd <- expand.grid(species_desc = unique(length_dat$species_desc),
  fork_length = seq(min(length_dat$fork_length), max(length_dat$fork_length), length.out = 200L),
  sarc_presence = c(0, 1))
nd$length_std <- (nd$fork_length - mean(length_dat$fork_length)) / sd(length_dat$fork_length)

p <- brms::posterior_linpred(fit1, newdata = nd)
nd$lwr <- apply(p, 2, quantile, probs = 0.05)
nd$est <- apply(p, 2, median)
nd$upr <- apply(p, 2, quantile, probs = 0.95)

set.seed(1)
length_dat$jit <- runif(nrow(length_dat), 0, 0.04)
length_dat$mature_jit <- ifelse(length_dat$mature == 1, length_dat$mature - length_dat$jit, length_dat$mature + length_dat$jit)

nd$sarc_presence_yn <- ifelse(nd$sarc_presence == 1, "Yes", "No")
length_dat$sarc_presence_yn <- ifelse(length_dat$sarc_presence == 1, "Yes", "No")

nd$species_desc <- as.character(nd$species_desc)
nd$species_desc[nd$species_desc == "rougheye/blackspotted rockfish complex"] <- "rougheye/blackspotted"
length_dat_plot <- length_dat
length_dat_plot$species_desc <- as.character(length_dat_plot$species_desc)
length_dat_plot$species_desc[length_dat_plot$species_desc == "rougheye/blackspotted rockfish complex"] <- "rougheye/blackspotted"

ggplot(nd, aes(fork_length, plogis(est), colour = sarc_presence_yn, ymin = plogis(lwr), ymax = plogis(upr), fill = sarc_presence_yn)) +
  geom_ribbon(colour = NA, alpha = 0.3) +
  geom_line() +
  facet_wrap(~stringr::str_to_title(species_desc)) +
  coord_cartesian(expand = FALSE, xlim = c(100, 750), ylim = c(-0.01, 1.01)) +
  geom_point(data = length_dat_plot, mapping = aes(fork_length, mature_jit), inherit.aes = FALSE, alpha = 0.03, pch = 21,size = 0.8) +
  ggsidekick::theme_sleek() + xlab("Fork length (mm)") + ylab("Probability of maturity") +
  theme(legend.title = ggtext::element_markdown()) +
  labs(colour = "*Sarcotaces* sp. present", fill = "*Sarcotaces* sp. present") +
  scale_fill_manual(values = c("No" = "grey50", "Yes" = "red")) +
  scale_colour_manual(values = c("No" = "grey50", "Yes" = "red")) +
  theme(legend.position = "top")

ggsave("draft_figs/prob-mature-species.pdf", width = 7, height = 4.6)
ggsave("draft_figs/prob-mature-species.png", width = 7, height = 4.6)

p <- brms::posterior_linpred(fit1, newdata = nd, re_formula = NA)
nd$lwr <- apply(p, 2, quantile, probs = 0.05)
nd$est <- apply(p, 2, median)
nd$upr <- apply(p, 2, quantile, probs = 0.95)

ggplot(nd, aes(fork_length, plogis(est), colour = factor(sarc_presence), ymin = plogis(lwr), ymax = plogis(upr), fill = factor(sarc_presence))) +
  geom_ribbon(colour = NA, alpha = 0.3) +
  geom_line() +
  coord_cartesian(expand = FALSE, xlim = c(100, 750), ylim = c(0, 1)) +
  # geom_point(data = length_dat, mapping = aes(fork_length, mature), inherit.aes = FALSE, alpha = 0.1, pch = 21,size = 1.5, position = position_jitter(width = 0, height = 0.02)) +
  ggsidekick::theme_sleek() + xlab("Fork length (mm)") + ylab("Probability of maturity") +
  theme(legend.title = ggtext::element_markdown()) +
  labs(colour = "*Sarcotaces* sp.<br>present", fill = "*Sarcotaces* sp.<br>present") +
  scale_fill_manual(values = c("0" = "grey50", "1" = "red")) +
  scale_colour_manual(values = c("0" = "grey50", "1" = "red"))+
  theme(legend.position = "inside", legend.position.inside = c(0.8, 0.5))

ggsave("draft_figs/prob-mature.pdf", width = 4, height = 3)
ggsave("draft_figs/prob-mature.png", width = 4, height = 3)

# post <- posterior_summary(fit1)

library(tidybayes)
get_variables(fit1)

p1 <- fit1 |> spread_draws(b_Intercept, b_length_std, b_sarc_presence, `b_length_std:sarc_presence`)

pivot_longer(p1, cols = b_Intercept:`b_length_std:sarc_presence`) |>
  filter(name != "b_Intercept") |>
  ggplot(aes(value)) +
  geom_vline(xintercept = 1) + ggsidekick::theme_sleek() +
  geom_density(fill = "grey90") + facet_wrap(~name, ncol = 1) +
  coord_cartesian(xlim = c(-3, 6), ylim = c(0, 1.6), expand = FALSE) +
  xlab("Coefficient estimate") + ylab("Density")
ggsave("draft_figs/maturity-fixed-coefs.png", width = 4, height = 4)

post_re <- fit1 |> spread_draws(r_species_desc[species_desc,term])
post_fe <- pivot_longer(p1, cols = b_Intercept:`b_length_std:sarc_presence`)
post_fe <- mutate(post_fe, term = gsub("^b_", "", name)) |> rename(fixed_coef_value = value)
post <- left_join(post_re, post_fe)

ggplot(post, aes(r_species_desc, species_desc)) + facet_wrap(~term)+
  geom_violin()

post |> filter(term != "Intercept") |>
  ggplot(aes(r_species_desc + fixed_coef_value, species_desc)) + facet_wrap(~term)+
  geom_vline(xintercept = 0) +
  # coord_cartesian(xlim = c(0, 5)) +
  # scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100)) +
  geom_violin(trim = TRUE)

post |> filter(term != "Intercept") |>
  mutate(combined = r_species_desc + fixed_coef_value) |>
  group_by(term, species_desc) |>
  summarise(lwr = quantile(combined, probs = 0.0275), upr = quantile(combined, probs = 0.975), mid = median(combined)) |>
  ggplot(aes(mid, species_desc, xmin = lwr, xmax = upr)) +
  facet_wrap(~term)+
  geom_vline(xintercept = 0, lty = 2) +
  # coord_cartesian(xlim = c(0, 5)) +
  # scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100)) +
  geom_point(pch = 21) +
  geom_linerange() + ggsidekick::theme_sleek() +
  xlab("Coefficient estimate") + theme(axis.title.y.left = element_blank())
ggsave("draft_figs/maturity-random-coefs.png", width = 7, height = 4)

get_p50 <- function(.int, .slope, .sd = sd(length_dat$fork_length), .mean = mean(length_dat$fork_length), p = 0.5) {
  xx <- -(log((1/p) - 1) + .int) / .slope
  (xx * .sd) + .mean
}

ggplot(nd, aes(fork_length, plogis(est), colour = factor(sarc_presence), ymin = plogis(lwr), ymax = plogis(upr), fill = factor(sarc_presence))) +
  geom_ribbon(colour = NA, alpha = 0.3) +
  geom_line() +
  coord_cartesian(expand = FALSE, xlim = c(100, 750), ylim = c(0, 1)) +
  # geom_point(data = length_dat, mapping = aes(fork_length, mature), inherit.aes = FALSE, alpha = 0.1, pch = 21,size = 1.5, position = position_jitter(width = 0, height = 0.02)) +
  ggsidekick::theme_sleek() + xlab("Fork length (mm)") + ylab("Probability of maturity") +
  theme(legend.title = ggtext::element_markdown()) +
  labs(colour = "*Sarcotaces* sp.<br>present", fill = "*Sarcotaces* sp.<br>present") +
  scale_fill_manual(values = c("0" = "grey50", "1" = "red")) +
  scale_colour_manual(values = c("0" = "grey50", "1" = "red"))+
  theme(legend.position = "inside", legend.position.inside = c(0.8, 0.5)) +
  geom_vline(xintercept = get_p50(2.85 - 1.47, 4.53)) +
  geom_hline(yintercept = 0.5)

post <- mutate(post, combined = r_species_desc + fixed_coef_value)

p50df <- group_by(post, species_desc) |>
  group_split() |>
  purrr::map_dfr(\(x) {
    intercept <- filter(x, term == "Intercept") |> pull(combined)
    length_slope <- filter(x, term == "length_std") |> pull(combined)
    sarc_adj <- filter(x, term == "sarc_presence") |> pull(combined)
    sarc_interaction <- filter(x, term == "length_std:sarc_presence") |> pull(combined)
    data.frame(
      p50_sarc_0 = get_p50(intercept, length_slope),
      p50_sarc_1 = get_p50(intercept + sarc_adj, length_slope + sarc_interaction),
      species_desc = x$species_desc[1]
    )
  })

p50df |> tidyr::pivot_longer(p50_sarc_0:p50_sarc_1) |>
  group_by(species_desc, name) |>
  summarise(lwr = quantile(value, probs = 0.0275), upr = quantile(value, probs = 0.975), mid = median(value)) |>
  ggplot(aes(mid, species_desc, colour = name, xmin = lwr, xmax = upr)) +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  ggsidekick::theme_sleek() +
  xlab("Length at p50 maturity") +
  theme(legend.position = "bottom")
ggsave("draft_figs/length-p50.png", width = 5, height = 5)

p50df |>
  mutate(diff = p50_sarc_1 - p50_sarc_0) |>
  ggplot(aes(diff, species_desc)) +
  geom_violin(trim = TRUE) +
  coord_cartesian(xlim = c(-100, 500)) +
  geom_vline(xintercept = 0, lty = 2)+
  ggsidekick::theme_sleek()

p50df |>
  mutate(diff = p50_sarc_1 - p50_sarc_0) |>
  group_by(species_desc) |>
  summarise(
    lwr = quantile(diff, probs = 0.05),
    lwr2 = quantile(diff, probs = 0.25),
    upr = quantile(diff, probs = 0.95),
    upr2 = quantile(diff, probs = 0.75),
    mid = median(diff)
  ) |>
  ggplot(aes(mid, forcats::fct_reorder(species_desc, -mid), xmin = lwr, xmax = upr)) +
  geom_linerange(lwd = 0.4) +
  geom_linerange(aes(xmin = lwr2, xmax = upr2), lwd = .7) +
  geom_point(pch = 19) +
  geom_vline(xintercept = 0, lty = 2)+
  theme(axis.title.y.left = element_blank()) +
  ggsidekick::theme_sleek() +
  xlab("Difference in<br>length at 50% maturity<br>if infected with *Sarcotaces* sp.<br>(mm)") +
  theme(axis.title = ggtext::element_markdown()) +
  ylab("")
ggsave("draft_figs/diff-p50.png", width = 4.5, height = 4)

group_by(p50df, species_desc) |>
  mutate(diff = p50_sarc_1 - p50_sarc_0) |>
  summarise(pdiff = mean(diff > 0)) |>
  ggplot(aes(pdiff, forcats::fct_reorder(species_desc, -pdiff))) +
  geom_hline(yintercept = unique(p50df$species_desc), lty = 2, col = "grey90") +
  geom_point(pch = 21) +
  ggsidekick::theme_sleek() + theme(axis.title.y.left = element_blank()) +
  coord_cartesian(xlim = c(0, 1)) +
  xlab("Probability that<br>length at 50% maturity<br>is greater for fish<br>infected with *Sarcotaces* sp.") +
  theme(axis.title.x = ggtext::element_markdown())
ggsave("draft_figs/prob-p50.png", width = 4.5, height = 4)
