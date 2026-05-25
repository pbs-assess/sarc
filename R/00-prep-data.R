library(dplyr)
library(tidyr)
dir.create("data-generated", showWarnings = FALSE)

# Don't need these, they are in the samples table, but looked at this to better understand
# maturity convetions
# mat_convention <- read_csv("data-raw/maturity-conventions.csv", na = c("", "NA", "NULL")) |>
#   janitor::clean_names()

full_df <- readRDS(here::here("data-raw/sarc_specimens.rds"))
full_df$year <- stringr::str_sub(full_df$FE_BEGIN_RETRIEVAL_TIME, 1, 4)

ssid0 <- readr::read_csv(here::here("data-raw/sarc-ssids.csv"), na = c("", "NA", "NULL"))

ssid <- ssid0 |>
  distinct(SAMPLE_ID, ACTIVITY_DESC, .keep_all = TRUE) |>
  mutate(ACTIVITY_DESC = ifelse(is.na(ACTIVITY_DESC), SURVEY_SERIES_DESC, ACTIVITY_DESC)) |>
  janitor::clean_names()

# Done: get survey id/desc, fishing event id (can get unique fishing _event_id, values from samples table)
# Done: will need to check if each fishing event has both encounters and non

names(full_df)[c(1, 5, 6, 8, 11, 12, 15, 18, 31:33, 35:37, 41:44)]

# Made this lookup in case we want names for the categories
# mat_lu <- readr::read_csv(here::here("data-raw", "maturity-lookup.csv")) |>
#   mutate(specimen_sex_desc = toupper(sex))

d0 <- full_df |>
  janitor::clean_names() |>
  left_join(ssid) |>
  mutate(year = lubridate::year(fe_begin_retrieval_time),
         species = tolower(species_desc),
         depth = ifelse(best_depth == 0, NA, best_depth), # replace 0 depths with NA (this is just how it is - I asked jon and maria)
         maturity_code = as.numeric(maturity_code),
         mature = ifelse(maturity_code >= 3, 1, 0),
         lat = as.numeric(best_lat),
         lon = -1 * as.numeric(best_long),
         species = gsub(" rockfish complex", "", species), # shorten rougheye/bs name
         fspecies = as.factor(species),
         sex = tolower(specimen_sex_desc)) |>
  select(year, fishing_event_id, lat, lon, gear_desc, depth,
         survey_series_id, survey_series_desc, activity_desc, trip_sub_type_desc,
         sample_id, specimen_id, species, fspecies,
         specimen_age, sex,
         mature, maturity_code, maturity_name, maturity_convention_desc,
         fork_length, round_weight,
         sarc_presence, sarc_count, sarc_comb) |>
  mutate(
    sarc_presence = case_when( # clarify NA values
      sarc_presence == "N" ~ 0,
      sarc_presence == "NULL" & sarc_count == 0 ~ 0,
      sarc_presence == "NULL" & sarc_count == "NULL" ~ NA_real_,
      sarc_presence == "Y" ~ 1
    ),
    sarc_count = case_when( # clarify NA values
      sarc_presence == 0 ~ 0,
      sarc_presence == NA & sarc_count == 0 ~ 0,
      sarc_presence == NA & sarc_count == "NULL" ~ NA_real_,
      sarc_presence == 1 & sarc_count == "NULL" ~ NA_real_,
      .default = as.numeric(sarc_count)
    ),
  ) |>
  tidyr::drop_na(sarc_presence) |>
  mutate(across(everything(), ~ ifelse(. == "NULL", NA, .))) |>
  mutate(across(c(specimen_age, fork_length, round_weight), \(x) as.numeric(x))) |>
  mutate(sarc_pres_label = factor(sarc_presence, levels = c(0, 1), labels = c("No", "Yes"))) |>
  tidyr::drop_na(survey_series_id) # omit commercial samples - only yelloweye and not clearly systematic
  # left_join(mat_lu)

saveRDS(d0, here::here("data-generated", "clean-data-all-years.rds"))

d <- d0 |>
  filter(year %in% 2019:2022) # in the other years, only presence was recorded
  # filter(year %in% c(1999, 2000, 2019:2022)) |>
  # drop_na(survey_series_desc)
saveRDS(d, here::here("data-generated", "clean-data.rds"))

e_spp <- d |> group_by(species) |>
  reframe(prop = sum(sarc_presence) / n(),
            n = n(),
            n_encounters = sum(sarc_presence)
            ) |>
  arrange(-n, prop)
saveRDS(e_spp, here::here("data-generated", "clean-data-encounter-summary.rds"))

# trip_sub_type_desc == "OBSERVED DOMESTIC" means that there was an onboard observer
# maturity_convention_desc =="PORT SAMPLES" means using the AMR maturity scale

# Clean length data
# ------------------------------------------------------------------------------
# Some visualisations to help figure out what data need cleaning for length data
# # Ogives
# ggplot(d, aes(x = fork_length, y = mature, colour = as.factor(sarc_presence))) +
#   geom_jitter(width = 0, height = 0.1, alpha = 0.4) +
#   facet_grid(sarc_presence ~ species) +
#   geom_smooth(method = "glm", method.args = list(family = binomial()), se = FALSE) +
#   theme(legend.position="bottom")

# Keep species that have fork lengths and both mature and immature individuals
keep_length <- d |>
  filter(!is.na(fork_length)) |>
  group_by(species, sarc_presence, mature) |>
  summarise(count = n(), .groups = "drop") |>
  pivot_wider(names_from = mature, values_from = count, values_fill = 0,
    names_prefix = "mature_") |>
  filter(`mature_0` > 0, `mature_1` > 0) |>
  group_by(species) |>
  filter(all(c(0, 1) %in% sarc_presence)) |>
  ungroup() |>
  distinct(species)

length_dat <- d |>
  filter(species %in% keep_length[[1]], !is.na(fork_length)) |>
  filter(!(maturity_code %in% 0), ! is.na(mature)) |> #
  mutate(species = factor(species),
         length_std = (fork_length - mean(fork_length)) / sd(fork_length))
saveRDS(length_dat, here::here("data-generated/length-dat.rds"))

# Clean age data
# ------------------------------------------------------------------------------
# Some visualisations to help figure out what data need cleaning for age data
# Ogives
# ggplot(d, aes(x = specimen_age, y = mature, colour = as.factor(sarc_presence))) +
#   geom_jitter(width = 0, height = 0.1, alpha = 0.4) +
#   facet_grid(sarc_presence ~ species) +
#   geom_smooth(method = "glm", method.args = list(family = binomial()), se = FALSE) +
#   theme(legend.position="bottom")

# Keep species that have ages and both mature and immature individuals
keep_age <- d |>
  filter(!is.na(specimen_age)) |>
  group_by(species, sarc_presence, mature) |>
  summarise(count = n(), .groups = "drop") |>
  pivot_wider(names_from = mature, values_from = count, values_fill = 0,
    names_prefix = "mature_") |>
  filter(`mature_0` > 0, `mature_1` > 0) |>
  group_by(species) |>
  filter(all(c(0, 1) %in% sarc_presence)) |>
  ungroup() |>
  distinct(species)

age_dat <- filter(d, species %in% keep_age[[1]], !is.na(specimen_age)) |>
  mutate(species = factor(species),
         age_std = (specimen_age - mean(specimen_age)) / sd(specimen_age))
saveRDS(age_dat, here::here("data-generated/age-dat.rds"))

dir.create(here::here("tables"), showWarnings = FALSE)
dir.create(here::here("figures"), showWarnings = FALSE)
dir.create(here::here("data-generated", "models"), showWarnings = FALSE, recursive = TRUE)
