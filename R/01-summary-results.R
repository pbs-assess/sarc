library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(see)
library(kableExtra)
library(flextable)

theme_set(ggsidekick::theme_sleek())

dat <- readRDS(here::here("data-generated", "clean-data-all-years.rds"))

# Summary data to include in report/notes for writing MS
sampling <- janitor::tabyl(dat, year, sarc_presence)
# sampling

encounter <- janitor::tabyl(dat, species, sarc_presence) |>
  mutate(sarc_group = case_when(
    `1` == 0 ~ "No sarcs obs",
    `0` == 0 ~ "Sarcs only",
    .default = "Both")) |>
  arrange(sarc_group, -`1`, -`0`)
saveRDS(encounter, here::here("data-generated", "encounter-spp-table-all-years.rds"))

clean_dat <- readRDS(here::here("data-generated", "clean-data.rds"))
# Summary data to include in report/notes for writing MS
clean_encounter <- janitor::tabyl(clean_dat, species, sarc_presence) |>
  mutate(sarc_group = case_when(
    `1` == 0 ~ "No sarcs obs",
    `0` == 0 ~ "Sarcs only",
    .default = "Both")) |>
  arrange(sarc_group, -`1`, -`0`)
saveRDS(clean_encounter, here::here("data-generated", "encounter-spp-table-systematic-years.rds"))

# note difference in systematic samples and those where opportunistic sampling found sarcs
sys_opp <- left_join(
  filter(encounter, `1` != 0) |> select(-sarc_group),
  filter(clean_encounter, `1` == 0) |> select(-sarc_group),
  by = c("species"),
  suffix = c("_all", "_sys")
) |>
  mutate(category = case_when(
    `1_sys` == 0 ~ "positives all time",
    `0_all` == 0 ~ "absence not recorded",
    TRUE ~ NA
  )) |>
  drop_na(category) |>
  mutate(text = paste0(str_to_title(species), " (", `1_all`, ")"))
saveRDS(sys_opp, here::here("data-generated", "encounter-compare-systematic-opportunistic.rds"))

# Depth distribution and summary of all specimens sampled for sarcs:
d <- clean_dat |>
  left_join(clean_encounter) |>
  mutate(fspecies = factor(stringr::str_to_title(species), levels = unique(stringr::str_to_title(clean_encounter$species)))) |>
  tidyr::drop_na(depth)
spp_filter <- clean_encounter |> filter(sarc_group == "Both") |> pull(species)

d0 <- d |> filter(sarc_presence == 0, species %in% spp_filter)
d1 <- d |> filter(sarc_presence == 1, species %in% spp_filter)

p <- ggplot(mapping = aes(x = fspecies, y = depth)) +
  geom_violinhalf(data = d0, fill = "grey50", colour = "grey30", alpha = 0.3,
    scale = "width", trim = FALSE, flip = TRUE) +
  geom_violinhalf(data = d1 |> filter(species != "bocaccio"), # bocaccio only has 2 ocurrences
    colour = "black", alpha = 0.9,
    scale = "width", trim = FALSE, flip = FALSE) +
  scale_y_reverse(limits = c(max(c(d0$depth, d1$depth)), -1), expand = expansion(mult = c(0, 0)), oob = scales::censor) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_point(data = d0, aes(x = fspecies), shape = 95, colour = "grey50", size = 2,
    position = position_nudge(x = -0.05)) +
  geom_point(data = d1, aes(x = fspecies), shape = 95, colour = "black", size = 2,
    position = position_nudge(x = 0.05)) +
  labs(x = "", y = "Depth (m)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.margin = margin(t = 10, r = 1, b = 1, l = 10))

p
ggsave(filename = here::here("figures", paste0("depth-plot.pdf")),
  width = 18.2 / cm(1), height = 5.6)
# ggsave(filename = here::here("figures", paste0("depth-plot.png")),
#   width = 18.2 / cm(1), height = 5.6)

# Frequency of cyst counts across species
count_levels <- factor(clean_encounter$species)
count_text <- clean_encounter |>
  filter(`1` > 0) |>
  mutate(fspecies = factor(stringr::str_to_title(species), levels = stringr::str_to_title(count_levels)),
         n = `0` + `1`) |>
  mutate(tile_fill = ifelse(row_number() %% 2 == 0, 1, 0))
count_dat <- d |>
  filter(sarc_presence == 1, !is.na(sarc_count))


# set.seed(40)
# set.seed(42)
# set.seed(94)
set.seed(99)
ggplot() +
  aes(y = forcats::fct_rev(fspecies)) +
  geom_jitter(data = count_dat |> filter(sarc_count < 5), aes(x = sarc_count),
    colour = "black", alpha = 0.3, shape = 21, width = 0.2, height = 0.2) +
  geom_jitter(data = count_dat |> filter(sarc_count >= 5 & sarc_count < 7), aes(x = sarc_count),
    colour = "black", alpha = 0.3, shape = 21, width = 0.05, height = 0.2) +
  geom_jitter(data = count_dat |> filter(sarc_count >= 7), aes(x = sarc_count),
    colour = "black", alpha = 0.3, shape = 21, width = 0.03, height = 0.1) +
  geom_text(data = count_text, aes(x = 0.4, label = paste0("(", n, ")")),
    hjust = 1, colour = "grey50", size = 2.5
  ) +
  scale_x_continuous(breaks = 1:10, limits = c(0.2, 10), expand = c(0.1, 0)) +
  guides(fill = "none") +
  theme(aspect.ratio = 1.2,
    axis.title.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.2)) +
  labs(
    x = "Number of cysts",
    y = "Species"
  )

ggsave(here::here("figures", "diff-sarc-count-by-species.pdf"), width = 5.5, height = 4.9)
# ggsave(here::here("figures", "diff-sarc-count-by-species.png"), width = 5.5, height = 4.9)
