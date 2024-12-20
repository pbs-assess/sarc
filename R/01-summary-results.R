library(dplyr)
library(ggplot2)
library(see)
library(kableExtra)
library(flextable)

theme_set(ggsidekick::theme_sleek())

dat <- readRDS(here::here("data-generated", "clean-data-all-years.rds"))

# Summary data to include in report/notes for writing MS
sampling <- janitor::tabyl(dat, year, sarc_presence)
sampling
encounter <- janitor::tabyl(dat, species, sarc_presence) |>
  mutate(sarc_group = case_when(
    `1` == 0 ~ "No sarcs obs",
    `0` == 0 ~ "Sarcs only",
    .default = "Both")) |>
  arrange(sarc_group, -`1`, -`0`)
saveRDS(encounter, here::here("data-generated", "encounter-spp-table.rds"))

d <- dat |>
  left_join(encounter) |>
  mutate(fspecies = factor(stringr::str_to_title(species), levels = unique(stringr::str_to_title(encounter$species))))

spp_filter <- encounter |> filter(sarc_group != "No sarcs obs") |> pull(species)
# Depth distribution and summary of all specimens sampled for sarcs:
d0 <- d |> filter(sarc_presence == 0, species %in% spp_filter)
d1 <- d |> filter(sarc_presence == 1, species %in% spp_filter)


# ggplot(mapping = aes(x = fspecies, y = depth)) +
#   geom_jitter(data = d0, shape = 21, colour = "grey50", alpha = 0.3) +
#   geom_jitter(data = d1, shape = 3, colour = "black", alpha = 0.9) +
#   scale_y_reverse() +
#   theme(axis.text.x = element_text(angle = 90))
p <- ggplot(mapping = aes(x = fspecies, y = depth)) +
  geom_violinhalf(data = d0, fill = "grey50", colour = "grey30", alpha = 0.3,
    scale = "width", trim = FALSE, flip = TRUE) +
  geom_violinhalf(data = d1, colour = "black", alpha = 0.9,
    scale = "width", trim = FALSE, flip = FALSE) +
  scale_y_reverse(limits = c(max(c(d0$depth, d1$depth)), -1), expand = expansion(mult = c(0, 0)), oob = scales::censor) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_point(data = d0, aes(x = fspecies), shape = 95, colour = "grey50",
    position = position_nudge(x = -0.05)) +
  geom_point(data = d1, aes(x = fspecies), shape = 95, colour = "black",
    position = position_nudge(x = +0.05)) +
  labs(x = "Rockfish species", y = "Depth (m)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(t = 10, r = 1, b = 1, l = 10))

p
ggsave(filename = here::here("figures", paste0("depth-plot.png")),
  width = 7.6, height = 6)

