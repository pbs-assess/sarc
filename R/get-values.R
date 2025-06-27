library(dplyr)
library(stringr)

mround <- function(x, digits) {
  sprintf(paste0("%.", digits, "f"), round(x, digits))
}

write_tex_comment <- function(comment) {
  paste0("% ", comment) |>
  readr::write_lines(here::here("text", "values.tex"), append = TRUE)
}

write_tex <- function(x, macro) {
  paste0("\\newcommand{\\", macro, "}{", x, "}") |>
  readr::write_lines(here::here("text", "values.tex"), append = TRUE)
}

# Create values.tex file
dir.create(here::here("text"), showWarnings = FALSE, recursive = TRUE)
suppressWarnings(file.remove(here::here("text", "values.tex")))

# Write initial header comment
write_tex_comment("Analysis Values")
write_tex_comment("===================")

clean_species_name <- function(species) {
  species |>
    tolower() |>
    gsub(" rockfish", "", x = _) |>
    gsub("rougheye/blackspotted", "REBS", x = _) |>
    gsub("pacific ocean perch", "POP", x = _) |>
    gsub("^([a-z])", "\\U\\1", x = _, perl = TRUE) # capitalize first letter
}

# Summary data values - sample sizes
write_tex_comment("\n% Sample Size Statistics")
write_tex_comment("-------------------")
e_table <- readRDS(here::here("data-generated", "encounter-spp-table-systematic-years.rds"))
n_sys_samples <- sum(c(e_table$`0`, e_table$`1`))

e_table_all_years <- readRDS(here::here("data-generated", "encounter-spp-table-all-years.rds"))
n_all_year_samples <- sum(e_table_all_years$`0`, e_table_all_years$`1`)

write_tex(scales::comma(n_all_year_samples), "nAllYearSamples") # systematic & opportunistic
write_tex(scales::comma(n_sys_samples), "nSysSamples") # systematic only

# Account for observations that were not sampled systematically
write_tex_comment("\n% Opportunistic Sampling Statistics")
write_tex_comment("-----------------------------")
sys_opp <- readRDS(here::here("data-generated",  "encounter-compare-systematic-opportunistic.rds"))

# Species where sarc cysts have been detected in opportunistic sampling but not systematic
sys_opp |>
  filter(category == "positives all time") |>
  mutate(
    species_clean = clean_species_name(species),
    macro = paste0("n", species_clean, "Opportunistic")
  ) |>
  group_walk(~ {
    write_tex(scales::comma(.x$`1_all`), .x$macro)
  })

# Species with only opportunistic sampling - number of infections
sys_opp |>
  filter(category == "absence not recorded") |>
  mutate(
    species_clean = clean_species_name(species),
    macro = paste0(species_clean, "nOpportunisticOnly")
  ) |>
  group_walk(~ {
    write_tex(scales::comma(.x$`1_all`), .x$macro)
  })

write_tex_comment("\n% Species Never Observed with Sarcs")
write_tex_comment("-----------------------------")
never_obs <- e_table_all_years |>
  filter(`1` == 0) |>
  pull(species) |>
  str_to_title() |>
  paste(collapse = ", ")
write_tex(never_obs, "SppNeverSarcs")

# N for species never observed with sarcs
e_table_all_years |>
  filter(`1` == 0) |>
  mutate(
    species_clean = clean_species_name(species),
    species_clean = gsub(" thornyhead", "", species_clean),
    macro = paste0("n", species_clean, "NeverObs")
  ) |>
  group_walk(~ {
    write_tex(scales::comma(.x$`0`), .x$macro)
  })

write_tex_comment("\n% Age Data Statistics")
write_tex_comment("-----------------")
# Age data sampling
ad <- readRDS(here::here("data-generated", "age-dat.rds")) |>
  mutate(species = as.character(species))
ad_table <- bind_rows(
  ad |>
    group_by(species) |>
    reframe(uninfected = sum(sarc_presence == 0), infected = sum(sarc_presence)) |>
    mutate(n = uninfected + infected, species = str_to_title(species)),
  ad |>
    reframe(uninfected = sum(sarc_presence == 0), infected = sum(sarc_presence)) |>
    mutate(species = "Total")
)

write_tex(length(unique(ad$species)), "nAgeSpecies")
write_tex(filter(ad_table, species != "Total") |> pull("species") |> paste(collapse = ", "), "ageSpeciesList")

ad_table |>
  filter(species != "Total") |>
  mutate(
    species_clean = clean_species_name(species),
    macro_uninfected = paste0("n", species_clean, "AgeUninfected"),
    macro_infected = paste0("n", species_clean, "AgeInfected")
  ) |>
  group_walk(~ {
    write_tex(scales::comma(.x$uninfected), .x$macro_uninfected)
    write_tex(scales::comma(.x$infected), .x$macro_infected)
  })

# Total counts for age data
write_tex(scales::comma(filter(ad_table, species == "Total") |> pull("uninfected")), "nTotalAgeUninfected")
write_tex(scales::comma(filter(ad_table, species == "Total") |> pull("infected")), "nTotalAgeInfected")


write_tex_comment("\n% Infection Rate Statistics")
write_tex_comment("---------------------")
# Infection rates for species
e_table |>
  mutate(
    species_clean = clean_species_name(species),
    p = round(100 * (`1` / (`0` + `1`)), 1),
    macro = paste0(species_clean, "InfectionRate")
  ) |>
  filter(`1` > 0) |>  # Only include species with infections
  group_walk(~ {
    write_tex(.x$p, .x$macro)
  })

# Number of species with infections in systematic sampling
write_tex(scales::comma(sum(e_table$`1` > 0)), "nSppWithInfections")

write_tex_comment("\n% Length Ogive Statistics - difference in length at 50\\% maturity")
write_tex_comment("-----------------------------")
# Length p50 differences - species and sex
p50_diff_summary <- readRDS(here::here("data-generated", "p50-length-diff-summary.rds"))
p50_diff_summary |>
  mutate(
    species_clean = clean_species_name(species),
    species_clean = gsub(" ", "", species_clean),
    sex_clean = paste0(toupper(substr(sex, 1, 1)), substr(sex, 2, nchar(sex))), # capitalize first letter
    macro_base = paste0("LDiffFifty", species_clean, sex_clean),
    value = mround(mid, 1),
    ci = paste0(mround(lwr, 1), ", ", mround(upr, 1))
  ) |>
  group_by(species, sex) |>
  group_walk(~ {
    write_tex(.x$value, paste0(.x$macro_base, "Median"))
    write_tex(.x$ci, paste0(.x$macro_base, "CI"))
  })

write_tex_comment("\n% Age Ogive Statistics - difference in age at 50\\% maturity")
p50_diff_age <- readRDS(here::here("data-generated", "p50-age-diff-summary.rds"))
p50_diff_age |>
  mutate(
    species_clean = clean_species_name(species),
    species_clean = gsub(" ", "", species_clean),
    sex_clean = paste0(toupper(substr(sex, 1, 1)), substr(sex, 2, nchar(sex))), # capitalize first letter
    macro_base = paste0("AgeDiffFifty", species_clean, sex_clean),
    value = mround(mid, 1),
    ci = paste0(mround(lwr, 1), ", ", mround(upr, 1))
  ) |>
  group_by(species, sex) |>
  group_walk(~ {
    write_tex(.x$value, paste0(.x$macro_base, "Median"))
    write_tex(.x$ci, paste0(.x$macro_base, "CI"))
  })

write_tex_comment("\n% Difference in infection rates between immature and mature")
write_tex_comment("----------------")
# Immature to mature ratio - population
imm_to_mat_ratio <- readRDS(here::here("data-generated", "immature-to-mature-infection-ratio.rds")) |>
  mutate(sex = str_to_title(sex))
write_tex(round(100 * mean(filter(imm_to_mat_ratio, sex == "Female" & species == "population")$imm_mat_ratio > 1), 0), "probRatioImmMatFemalesGreater")
write_tex(round(100 * mean(filter(imm_to_mat_ratio, sex == "Male" & species == "population")$imm_mat_ratio > 1), 0), "probRatioImmMatMalesGreater")
write_tex(round(100 * mean(filter(imm_to_mat_ratio, sex == "Female" & species == "population")$imm_mat_ratio > 2), 0), "probRatioImmMatFemalesTwoXGreater")
write_tex(round(100 * mean(filter(imm_to_mat_ratio, sex == "Male" & species == "population")$imm_mat_ratio > 2), 0), "probRatioImmMatMalesTwoXGreater")

imm_mat_ratio_df <- imm_to_mat_ratio |>
  arrange(species, sex) |>
  group_by(species, sex) |>
  summarise(median = median(imm_mat_ratio),
            lwr = quantile(imm_mat_ratio, probs = 0.05),
            upr = quantile(imm_mat_ratio, probs = 0.95)) |>
  mutate(
    species_clean = clean_species_name(species),
    macro_base = paste0("RatioImmMat", species_clean, sex),
    value = mround(median, 1),
    ci = paste0(mround(lwr, 1), " to ", mround(upr, 1))
  ) |>
  group_walk(~ {
    write_tex(.x$value, paste0(.x$macro_base, "Median"))
    write_tex(.x$ci, paste0(.x$macro_base, "CI"))
  })

