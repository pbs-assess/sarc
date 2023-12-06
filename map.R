library(dplyr)
library(ggplot2)
library(sdmTMB)

d <- readr::read_csv("sarc_specimens.csv", show_col_types = FALSE)
vroom::problems(d)
d$year <- as.numeric(stringr::str_sub(d$FE_BEGIN_RETRIEVAL_TIME, 1, 4))

names(d) <- tolower(names(d))
glimpse(d)

d$best_lat <- as.numeric(d$best_lat)
d$best_long <- as.numeric(d$best_long)

d <- filter(d, !is.na(best_lat))
d <- filter(d, !is.na(best_long))
# !?

d$lon <- d$best_long * -1
d$lat <- d$best_lat

table(d$sarc_presence)

d <- d %>%
  filter(sarc_comb != "na") |>
  mutate(sarc_presence = case_when(
    sarc_comb %in% c("0", "N") ~ 0,
    sarc_comb %in% c("Y", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10") ~ 1
  )) |>
  mutate(sarc_presence = as.factor(sarc_presence))

table(d$sarc_presence)

ggplot(d, aes(lon, lat, colour = sarc_presence, pch = sarc_presence)) +
  geom_point(alpha = 0.2) +
  geom_point(data = filter(d, sarc_presence == 1), colour = "black")

d <- filter(d, !is.na(year))

d <- add_utm_columns(d, c("lon", "lat"))

d <- filter(d, year >= 1999)

mesh <- make_mesh(d, c("X", "Y"), cutoff = 5)
plot(mesh)

d$sarc_presence <- as.numeric(as.character(d$sarc_presence))
table(d$sarc_presence)
table(d$sarc_presence)

d$log_best_depth <- log(d$best_depth + 1)

fit <- sdmTMB(
  sarc_presence ~ 1 + poly(log_best_depth, 2) + s(year, k = 5),
  # sarc_presence ~ 1,
  data = d,
  family = binomial(), spatial = "on", silent = FALSE,
  mesh = mesh)

fit

ggeffects::ggeffect(fit, "log_best_depth [all]") |> plot()
ggeffects::ggpredict(fit, terms = "year [all]") |> plot()

grid <- gfplot::synoptic_grid |> select(X, Y, best_depth = depth)
gridll <- gfplot::hbll_grid$grid |> select(lon = X, lat = Y, best_depth = depth) |>
  add_utm_columns(c("lon", "lat"))
# gridlli <- gfplot::hbll_inside_n_grid$grid |> select(X, Y, best_depth = depth)
grid <- bind_rows(grid, gridll)
grid$log_best_depth <- log(grid$best_depth + 1)

ggplot(d, aes(lon, lat, colour = sarc_presence)) +
  geom_point(data = filter(d, sarc_presence == 0), colour = "black") +
  geom_point(data = filter(d, sarc_presence == 1), colour = "red") +
  coord_fixed()

grid$year <- 2020
p <- predict(fit, newdata = grid)

ggplot(p, aes(X, Y, fill = plogis(est))) +
  geom_tile(width = 2, height = 2) +
  scale_fill_viridis_c() +
  coord_fixed() +
  theme_light()

ggplot(p, aes(X, Y, fill = omega_s)) +
  geom_tile(width = 2, height = 2) +
  scale_fill_gradient2() +
  coord_fixed() +
  theme_light()

group_by(d, year) |>
  summarise(prop = mean(as.integer(as.character(sarc_presence))), n = n()) |>
  as.data.frame()
