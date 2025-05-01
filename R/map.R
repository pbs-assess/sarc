library(dplyr)
library(ggplot2)
library(sdmTMB)

# Load data
d0 <- readRDS(here::here("data-generated", "clean-data.rds"))

d <- d0 |>
  tidyr::drop_na(lat) |>
  add_utm_columns(c("lon", "lat"), utm_crs = 32609) |>
  filter((lon > -135), fishing_event_id != 5751520) |> # weird outliers way off the coast
  filter(year >= 2019) |> # once observations become more consistent
  mutate(log_depth = log(depth))

# A few checks:
# Zero depths removed
min(d$depth)
# Only 7 fishing events in 1999 & 2000
d0 |> filter(year %in% 1999:2000) |> distinct(fishing_event_id) |> nrow()

# Fit spatial model
mesh <- make_mesh(d, c("X", "Y"), cutoff = 5)
plot(mesh)
d$year_f <- as.factor(d$year)
fit <- sdmTMB(
  sarc_presence ~ 1 + poly(log_depth, 2) + year_f,
  # sarc_presence ~ 1 + poly(log_depth, 2) + s(year, k = 5),
  # sarc_presence ~ 1,
  # time_varying = ~1,
  # time_varying_type = "rw",
  # time = "year",
  # spatiotemporal = "off",
  data = d,
  anisotropy = TRUE,
  family = binomial(), spatial = "on", silent = TRUE,
  mesh = mesh)
sanity(fit)
fit
plot_anisotropy(fit)

ggeffects::ggeffect(fit, "log_depth [all]") |> plot()
ggeffects::ggpredict(fit, terms = "year_f [all]") |> plot()

nd <- data.frame(log_depth = seq(min(d$log_depth), max(d$log_depth), length.out = 200),
                 year_f = "2022")
pp <- predict(fit, newdata = nd, se_fit = TRUE, re_form = NA)
ggplot(pp, aes(exp(log_depth), plogis(est), ymin = plogis(est - 2 * est_se), ymax = plogis(est + 2 * est_se))) +
  geom_ribbon(fill = "grey80") + geom_line() +
  scale_x_continuous(trans = "log10", expand = expansion(mult = c(0.02, 0.02), add = c(0, 0))) +
  theme_light() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02), add = c(0, 0)))

# ----------
# Prepare grid
# grid <- gfplot::synoptic_grid |> select(X, Y, depth)
# gridll <- gfplot::hbll_grid$grid |> select(lon = X, lat = Y, depth) |>
#   add_utm_columns(c("lon", "lat"))

# gridlli <- gfplot::hbll_inside_n_grid$grid |> select(X, Y, depth)
# grid <- bind_rows(grid, gridll)

# the new survey_blocks grid has a lot more cells than the old gfplot:: ones
g0 <- gfdata::survey_blocks |>
  filter(grepl("SYN|HBLL", survey_abbrev),
        active_block == TRUE) |>
  st_centroid()
g_coords <- st_coordinates(g0)
g <- g0 |>
  select(survey_abbrev, area, depth = depth_m) |>
  mutate(X = g_coords[, 1] / 1000, Y = g_coords[, 2] / 1000) |>
  st_drop_geometry() |>
  as_tibble()

# There are some weird depths in the new grids
g |> filter(depth == 0)
g |> filter(depth < 0)

g <- filter(g, depth > 0, depth < 1300) # exclude grids less than 0...
# Also filter out grids that are deeper tha 1300? The original grids were maxed out at 1300

g$log_depth <- log(g$depth)
g$year_f <- d$year_f[1]

ggplot(d, aes(lon, lat, colour = sarc_presence)) +
  geom_point(data = filter(d, sarc_presence == 0), colour = "black") +
  geom_point(data = filter(d, sarc_presence == 1), colour = "red") +
  coord_fixed()

g$year <- 2022
g$year_f <- factor(g$year, levels = levels(d$year_f))
p <- predict(fit, newdata = g)

# x <- seq(-100, 100, length.out = 100)
# b1 <- 0.05
# b2 <- -0.002
# y <- exp(1 + x*b1 + b2*x^2)
# plot(x, y)
#
# y20 <- max(y) * 0.2
# y20
#
# plot(x, y);abline(h = y20)
# first <- x[which(y >= y20)[1]]
#
# abline(v = first)
# xrev <- rev(x)
# yrev <- rev(y)
# last <- xrev[which(yrev >= y20)[1]]
# abline(v = last)
# last - first
#


qlogis_trans <- function() {
  scales::trans_new(
    name = "qlogis",
    transform = function(x) qlogis(x),
    inverse = function(x) plogis(x),
    minor_breaks = function(x) c
    domain = c(0, 1))
}

map_data <- rnaturalearth::ne_countries(
  scale = "large",
  returnclass = "sf", country = "canada")
bc_coast <- suppressWarnings(suppressMessages(
  sf::st_crop(map_data,
    c(xmin = -134, ymin = 45, xmax = -110, ymax = 57))))
utm_zone9 <- 3156
bc_coast_proj <- sf::st_transform(bc_coast, crs = utm_zone9)

library(ggtext) # for italics

g1 <- ggplot(bc_coast_proj) + geom_sf() +
  geom_tile(width = 2000, height = 2000, data = p, mapping = aes(X * 1000, Y * 1000, fill = plogis(est))) +
  scale_fill_viridis_c(option = "F", trans = "qlogis", breaks = c(0.01, 0.1, 0.5, 0.9, 0.99)) +
  coord_sf() +
  theme_light() +
  # geom_point(data = filter(d, sarc_presence == 0), mapping = aes(X, Y), inherit.aes = FALSE, colour = "white", alpha = 0.01, size = 1, pch = 4) +
  # geom_point(data = filter(d, sarc_presence == 1), mapping = aes(X * 1000, Y * 1000), inherit.aes = FALSE, colour = "white", alpha = 0.1, size = 2, pch = 4, fill = NA) +
  xlim(230957.7 - 80000, 1157991 - 220000) +
  ylim(5366427 - 25000, 6353456 - 250000) +
  theme(legend.title = element_markdown()) +
  labs(fill = "Probability of<br>*Sarcotaces* sp.<br>encounter", x = "Longitude", y = "Latitude") +
  theme(legend.position = "bottom")
# g1

d$present <- ifelse(d$sarc_presence == "1", "Yes", "No")
g2 <- ggplot(bc_coast_proj) + geom_sf() +
  geom_point(data = arrange(d, present), mapping = aes(X * 1000, Y*1000, colour = present, shape = present), inherit.aes = FALSE, alpha = 0.3, size = 1) +
  theme_light() +
  scale_shape_manual(values = c("No" = 4, "Yes" = 21)) +
  scale_colour_manual(values = c("No" = "grey60", "Yes" = "red")) +
  xlim(230957.7 - 80000, 1157991 - 220000) +
  ylim(5366427 - 25000, 6353456 - 250000)+
  theme(legend.title = element_markdown()) +
  labs(x = "Longitude", y = "Latitude", colour = "*Sarcotaces* sp.<br>encounter", shape = "*Sarcotaces* sp.<br>encounter") +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
# g2

library(patchwork)
g <- g2 + g1
ggsave("draft_figs/map.png", width = 9.5, height = 5.25)

ggplot(p, aes(X, Y, fill = est)) +
  geom_tile(width = 2, height = 2) +
  scale_fill_viridis_c(option = "C") +
  coord_fixed() +
  theme_light()

ggplot(p, aes(X, Y, fill = omega_s)) +
  geom_tile(width = 2, height = 2) +
  scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red")) +
  coord_fixed() +
  theme_light()

group_by(d, year) |>
  summarise(prop = mean(as.integer(as.character(sarc_presence))), n = n()) |>
  as.data.frame()

# TODO
# - [ ] better grid
# - [ ] early years!? all 1s
# - [ ] same depth assignment
