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
ld <- readRDS(here::here("data-generated", "length-dat.rds"))
cld <- ld |>
  filter(species %in% spp_filter) |>
  filter(year >= 2019) |> # Did not seem to count sarcs in 1999/2000
  drop_na(sarc_count) |>
  mutate(log_sarc_count = log(sarc_count + 0.01))

dat <- cld |> filter(sarc_count > 0)
m1 <- glmmTMB(
  mature ~ length_std * log(sarc_count) + (1 + length_std | species),
  family = binomial(),
  data = dat
)
summary(m1)

m2 <- glmmTMB(
  mature ~ length_std * sarc_count + (1 + length_std | species),
  family = binomial(),
  data = dat
)
summary(m2)

# so many zeroes I don't think it matters
par(mfrow = c(1, 2)); plot(residuals(m1)); plot(residuals(m2))

DHARMa::simulateResiduals(fittedModel = m1, plot = T)
DHARMa::simulateResiduals(fittedModel = m2, plot = T)


# Maturity ~ count brms
# ------------------------
# Not really enough data to look at species effects of sarc count on maturity,
# but I think there is enough to make a general statement about what it does on
# average - that a higher intensity of infection does likely have a negative
# effect on maturity at length/age
# Done: filter out n_encounters == 0
# TODO look at species level estimates of the parameters - specifically sarc count and sarc count * length to see if there is species variation
# TODO plot residuals to make sure that the raw sarc counts are better, might need to have log(sarc_count)


# fit_2f_c <- brm(
#   mature ~ 0 + Intercept + length_std * sarc_count +
#     (1 + length_std * sarc_count | species),
#   family = bernoulli(),
#   data =
#     cld |>
#       filter(sex == "female"),
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
# saveRDS(fit_2f_c, file.path(fit_dir, "maturity-length-count-stan-female.rds"))
# beepr::beep()

# saveRDS(fit_2f_lc, file.path(fit_dir, "maturity-length-count-stan-female-lc.rds"))
# beepr::beep()

fit_2m_lc <- update(fit_2f_lc, newdata = cld |> filter(sex == "male"))
saveRDS(fit_2m_lc, file.path(fit_dir, "maturity-length-count-stan-male-lc.rds"))
beepr::beep()
fit_2m_lc <- readRDS(file.path(fit_dir, "maturity-length-count-stan-male-lc.rds"))


fit_2f_c <- readRDS(file.path(fit_dir, "maturity-length-count-stan-female.rds"))
fit_2f_lc <- readRDS(file.path(fit_dir, "maturity-length-count-stan-female-lc.rds"))

# Compare log transform on sarc count
# https://users.aalto.fi/~ave/modelselection/diabetes.html#4_A_Bayesian_logistic_regression_model
# Predicted probabilities
fit <- fit_2f
linpred <- posterior_linpred(fit)
preds <- posterior_epred(fit)
pred <- colMeans(preds)
pr <- as.integer(pred >= 0.5)
loo1 <- loo(fit, save_psis = TRUE)
y <- test$mature

# posterior classification accuracy
round(mean(xor(pr,as.integer(y==0))),2)
round((mean(xor(pr[y==0]>0.5,as.integer(y[y==0])))+mean(xor(pr[y==1]<0.5,as.integer(test$mature[test$mature==1]))))/2,2)

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
    loo_obj = loo_obj,
    loo_acc = loo_acc, loo_bal_acc = loo_bal_acc,
    cal_plot_data = cal_plot_data
  ))
}

# Evaluate the two models
y <- fit_2f_lc$data$mature
result_2f_c <- evaluate_model(fit_2f_c, y); beepr::beep()
result_2f_lc <- evaluate_model(fit_2f_lc, y); beepr::beep()

loo_compare(result_2f_c$loo_obj, result_2f_lc$loo_obj)

# # fit <- fit_2f_lc
# fit <- fit_2f

# test <- fit$data

# fitted_values <- posterior_epred(fit)
# predicted_values <- posterior_predict(fit)
# observed_values <- fit$data$mature
# mean_fitted <- colMeans(fitted_values)
# mean_observed <- mean(observed_values)

# test$fitted_prob <- colMeans(fitted_values)
# p1 <- ggplot(test, aes(x = fitted_prob, y = mature)) +
#   geom_jitter(height = 0.05, width = 0) +
#   geom_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE) +
#   labs(x = "Fitted Probabilities", y = "Observed Outcomes",
#        title = paste0("Fitted vs Observed: ", tag))
# p1

# test$fitted_prob <- colMeans(fitted_values)
# test$bin <- cut(test$fitted_prob, breaks = seq(0, 1, by = 0.1))
# calibration <- test |>
#   group_by(bin) |>
#   summarize(mean_pred = mean(fitted_prob),
#             mean_obs = mean(mature))

# p2 <- ggplot(calibration, aes(x = mean_pred, y = mean_obs)) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#   labs(x = "Mean Predicted Probability", y = "Mean Observed Outcome",
#        title = paste0("Calibration Plot: ", tag))
# p2

# test$fitted_prob <- colMeans(fitted_values)
# group_summary <- test |>
#   group_by(species) |>
#   summarize(mean_fitted = mean(fitted_prob), mean_observed = mean(mature))

# p3 <- ggplot(group_summary, aes(x = mean_fitted, y = mean_observed, colour = species)) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#   labs(x = "Mean Fitted Probability", y = "Mean Observed Outcome",
#        title = paste0("Group-Level Fitted vs Observed: ", tag)) +
#   guides(colour = "none")
# p3

# p1 + p2 + p3



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