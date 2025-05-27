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

e_spp <- readRDS(here::here("data-generated", "clean-data-encounter-summary.rds"))
spp_filter <- e_spp |> filter(n_encounters > 0) |> pull(species)
ld <- readRDS(here::here("data-generated", "length-dat.rds")) |>
  filter(sex %in% c("female", "male"))
cld <- ld |>
  filter(species %in% spp_filter) |>
  drop_na(sarc_count) |>
  mutate(log_sarc_count = log(sarc_count + 0.1),
         sex_f = as.factor(sex))

# dat <- cld |> filter(sarc_count > 0)
# m1 <- glmmTMB(
#   mature ~ length_std * log(sarc_count) + (1 + length_std | species),
#   family = binomial(),
#   data = dat
# )
# summary(m1)

# m2 <- glmmTMB(
#   mature ~ length_std * sarc_count + (1 + length_std | species),
#   family = binomial(),
#   data = dat
# )
# summary(m2)

# # so many zeroes I don't think it matters
# par(mfrow = c(1, 2)); plot(residuals(m1)); plot(residuals(m2))

# DHARMa::simulateResiduals(fittedModel = m1, plot = T)
# DHARMa::simulateResiduals(fittedModel = m2, plot = T)


# Maturity ~ count brms
# ------------------------
# Not really enough data to look at species effects of sarc count on maturity,
# but I think there is enough to make a general statement about what it does on
# average - that a higher intensity of infection does likely have a negative
# effect on maturity at length/age
# Done: filter out n_encounters == 0 (i.e., only species with at least one infected individual are included in this dataset)
# TODO look at species level estimates of the parameters - specifically sarc count and sarc count * length to see if there is species variation
# TODO plot residuals to make sure that the raw sarc counts are better, might need to have log(sarc_count)

# Try full model, theoretically should improve power no?
# fit_full_lc <- brm(
#   mature ~ 0 + Intercept + length_std * log_sarc_count * sex_f +
#     (1 + length_std * log_sarc_count * sex_f | species),
#   family = bernoulli(),
#   data = cld,
#   iter = 2000L,
#   warmup = 1000L,
#   chains = 4L,
#   cores = 4L,
#   backend = "cmdstanr",
#   prior =
#     prior(normal(0, 5), class = b) +
#     prior(student_t(3, 0, 2), class = sd) +
#     prior(normal(0, 10), class = b, coef = Intercept),
#   control = list(max_treedepth = 12, adapt_delta = 0.95)
# )
# saveRDS(fit_full_lc, file.path(fit_dir, "maturity-length-count-stan-sex-combined-log-count.rds"))
# beepr::beep()
fit_full_lc <- readRDS(file.path(fit_dir, "maturity-length-count-stan-sex-combined-log-count.rds"))

options(mc.cores = parallel::detectCores() - 2)
fit_fc <- brm(
  mature ~ 0 + Intercept + length_std * sarc_count +
    (1 + length_std * sarc_count | species),
  family = bernoulli(),
  data = cld |> filter(sex == "female"),
  iter = 2000L,
  warmup = 1000L,
  chains = 4L,
  cores = 4L,
  backend = "cmdstanr",
  prior =
    prior(normal(0, 5), class = b) +
    prior(student_t(3, 0, 2), class = sd) +
    prior(normal(0, 10), class = b, coef = Intercept),
  control = list(max_treedepth = 12, adapt_delta = 0.95)
)
saveRDS(fit_fc, file.path(fit_dir, "maturity-length-count-stan-female.rds"))
beepr::beep()
fit_fc <-readRDS(file.path(fit_dir, "maturity-length-count-stan-female.rds"))

# Female count using log sarc count
# fit_flc <- update(fit_fc,
#   formula = mature ~ 0 + Intercept + length_std * log_sarc_count +
#     (1 + length_std * log_sarc_count | species),
#   newdata = cld |> filter(sex == "female"))
# saveRDS(fit_flc, file.path(fit_dir, "maturity-length-count-stan-female-lc.rds"))
# beepr::beep()


# Male count model
# fit_mc <- update(fit_fc, newdata = cld |> filter(sex == "male"))
# saveRDS(fit_mc, file.path(fit_dir, "maturity-length-count-stan-male.rds"))
#beepr::beep()
fit_mc <- readRDS(file.path(fit_dir, "maturity-length-count-stan-male.rds"))

# Male count using log sarc count
# fit_mlc <- update(fit_fc,
#   formula = mature ~ 0 + Intercept + length_std * log_sarc_count +
#     (1 + length_std * log_sarc_count | species),
#   newdata = cld |> filter(sex == "male"))
# saveRDS(fit_mlc, file.path(fit_dir, "maturity-length-count-stan-male-lc.rds"))
# beepr::beep()
# fit_mlc <- readRDS(file.path(fit_dir, "maturity-length-count-stan-male-lc.rds"))


fit_fc <- readRDS(file.path(fit_dir, "maturity-length-count-stan-female.rds"))
fit_flc <- readRDS(file.path(fit_dir, "maturity-length-count-stan-female-lc.rds"))

# Compare log transform on sarc count
# Try looking at RQR? - I don't think log matters because the effect is so small? #Question
# https://cran.r-project.org/web/packages/tidybayes/vignettes/tidybayes-residuals.html
test <- cld |> filter(sex == "female") |>
  add_predicted_draws(fit_fc)
test2 <- test |> make_probability_residuals(.prediction, mature, n = 50)
p1 <- test2 |>
  ggplot(aes(sample = .p_residual)) +
  geom_qq(distribution = qunif) +
  geom_abline()

test_lc <- cld |>
  add_predicted_draws(fit_full_lc); beepr::beep()
test2_lc <- test_lc |> make_probability_residuals(.prediction, mature, n = 50)
p2 <- test2_lc |>
  ggplot(aes(sample = .p_residual)) +
  geom_qq(distribution = qunif) +
  geom_abline() +
  ggtitle("log count - both sexes in model")

p1 + p2

# https://users.aalto.fi/~ave/modelselection/diabetes.html#4_A_Bayesian_logistic_regression_model
# Predicted probabilities
fit <- fit_full_lc
linpred <- posterior_linpred(fit)
preds <- posterior_epred(fit)
pred <- colMeans(preds)
pr <- as.integer(pred >= 0.5)
loo1 <- loo(fit, save_psis = TRUE)
y <- fit$data$mature

# posterior classification accuracy
round(mean(xor(pr,as.integer(y==0))),2)
round((mean(xor(pr[y==0]>0.5,as.integer(y[y==0])))+mean(xor(pr[y==1]<0.5,as.integer(y[y==1]))))/2,2)

# LOO predictive probabilities
ploo <- loo::E_loo(preds, loo1$psis_object, type="mean", log_ratios = -log_lik(fit))$value
# LOO classification accuracy
round(mean(xor(ploo>0.5,as.integer(y==0))),2)
# LOO balanced classification accuracy
round((mean(xor(ploo[y==0]>0.5,as.integer(y[y==0])))+mean(xor(ploo[y==1]<0.5,as.integer(y[y==1]))))/2,2)
qplot(pred, ploo)

calPlotData <- caret::calibration(as.factor(y) ~ pred + loopred,
                         data = data.frame(pred = pred,loopred = ploo,y = y),
                         cuts = 10, class = "1")
ggplot(calPlotData, auto.key = list(columns = 2))+
  geom_jitter(data = data.frame(pred = pred,loopred = ploo,y = (as.numeric(y))*100), inherit.aes = FALSE,
              aes(x = loopred*100, y = y), height = 2, width = 0, alpha = 0.3) +
  scale_colour_brewer(palette = "Set1")+
  bayesplot::theme_default(base_family = "sans")

# ---
# Define a function for model evaluation
evaluate_model <- function(fit, y) {
  linpred <- posterior_linpred(fit)
  preds <- posterior_epred(fit)
  pred <- colMeans(preds)
  pr <- as.integer(pred >= 0.5)
  loo_obj <- loo(fit, save_psis = TRUE)

  # Posterior classification accuracy
  post_acc <- round(mean(xor(pr, as.integer(y == 0))), 2)
  post_bal_acc <- round((mean(xor(pr[y == 0] > 0.5, as.integer(y[y == 0]))) +
                         mean(xor(pr[y == 1] < 0.5, as.integer(y[y == 1])))) / 2, 2)

  # LOO predictive probabilities
  ploo <- loo::E_loo(preds, loo_obj$psis_object, type = "mean", log_ratios = -log_lik(fit))$value

  # LOO classification accuracy
  loo_acc <- round(mean(xor(ploo > 0.5, as.integer(y == 0))), 2)
  loo_bal_acc <- round((mean(xor(ploo[y == 0] > 0.5, as.integer(y[y == 0]))) +
                        mean(xor(ploo[y == 1] < 0.5, as.integer(y[y == 1])))) / 2, 2)

  # Calibration plot data
  cal_plot_data <- caret::calibration(as.factor(y) ~ pred + ploo,
                                      data = data.frame(pred = pred, ploo = ploo, y = y),
                                      cuts = 10, class = "1")

  return(list(
    preds = preds, pred = pred, pr = pr, ploo = ploo,
    post_acc = post_acc, post_bal_acc = post_bal_acc,
    loo_obj = loo_obj, loo_acc = loo_acc, loo_bal_acc = loo_bal_acc,
    cal_plot_data = cal_plot_data
  ))
}

# Evaluate the two models
y1 <- fit_fc$data$mature
y2 <- fit_full_lc$data$mature
result_fc <- evaluate_model(fit_fc, y); beepr::beep()
result_flc <- evaluate_model(fit_full_lc, fit_full_lc$data$mature); beepr::beep()

calPlotData1 <- with(result_fc, caret::calibration(as.factor(y) ~ pred + loopred,
                         data = data.frame(pred = pred,loopred = ploo,y = y),
                         cuts = 10, class = "1"))
calPlotData2 <- with(result_flc, caret::calibration(as.factor(y2) ~ pred + loopred,
                         data = data.frame(pred = pred,loopred = ploo,y = y2),
                         cuts = 10, class = "1"))

p1 <- ggplot(calPlotData1, auto.key = list(columns = 2))+
  geom_jitter(data = data.frame(pred = pred,loopred = ploo,y = (as.numeric(y))*100), inherit.aes = FALSE,
              aes(x = loopred*100, y = y), height = 2, width = 0, alpha = 0.3) +
  scale_colour_brewer(palette = "Set1")+
  bayesplot::theme_default(base_family = "sans") +
  ggtitle("raw count")

p2 <- ggplot(calPlotData2, auto.key = list(columns = 2))+
  geom_jitter(data = data.frame(pred = pred,loopred = ploo,y = (as.numeric(y))*100), inherit.aes = FALSE,
              aes(x = loopred*100, y = y), height = 2, width = 0, alpha = 0.3) +
  scale_colour_brewer(palette = "Set1")+
  bayesplot::theme_default(base_family = "sans") +
  ggtitle("log count")

p1 + p2


# Visualise results
# ------------------
table(cld$sarc_count) |>
enframe() |>
pivot_wider(names_from = name, values_from = value) |>
flextable::flextable()

cld |>
  filter(sarc_presence == 1) |>
  ggplot(aes(x = sarc_count, y = sex)) +
  geom_jitter(width = 0.2, height = 0.2) +
  scale_x_continuous(breaks = 1:10) +
  theme_minimal() +
  facet_grid(species ~ ., switch = "y", scales = "free_y") +
  theme(
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    strip.placement = "outside"
  ) +
  labs(x = "Sarc Count", y = "Sex")

# Fixed effects
post_fc <- fit_fc |> spread_draws(b_Intercept, b_length_std, b_sarc_count, `b_length_std:sarc_count`)
post_fc <- post_fc |>
  pivot_longer(cols = b_Intercept:`b_length_std:sarc_count`, values_to = "fe_coef") |>
  mutate(term = gsub("^b_", "", name)) #|>
  # mutate(term = gsub("length_std", "Length", term),
  #        term = gsub("sarc_count", "Sarcotaces count", term)) |>
  # mutate(term = factor(term, levels = c("Intercept", "Length", "Sarcotaces count", "Length:Sarcotaces count")))

pd_count_fe <- post_fc |>
  ggplot(aes(x = fe_coef)) +
  ggsidekick::theme_sleek() +
  geom_density(fill = "grey90") +
  geom_vline(xintercept = 0) +
  facet_wrap(~ term, ncol = 1) +
  # coord_cartesian(xlim = c(-3, 6), ylim = c(0, 2), expand = FALSE) +
  xlab("Coefficient estimate") + ylab("Posterior density")
pd_count_fe

posterior_draws <- post_flc |>
  select(term, fe_coef) |>
  pivot_wider(names_from = term, values_from = fe_coef) |>
  as.matrix()

# Plot with bayesplot
bayesplot::mcmc_areas(
  posterior_draws,
  pars = colnames(posterior_draws),
  prob = 0.8,  # 80% credible interval
  prob_outer = 0.95,  # 95% credible interval
  point_est = "mean"
) +
  ggtitle("Posterior Density for Fixed Effects Coefficients") +
  theme_sleek() +
  labs(x = "Coefficient Estimate", y = "Density")

ggsave(here::here("figures", "maturity-length-count-posterior-density.png"),
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