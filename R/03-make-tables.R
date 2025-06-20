library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(flextable)
library(kableExtra)
library(scales)

# Create values.tex file
dir.create(here::here("tables"), showWarnings = FALSE, recursive = TRUE)
suppressWarnings(file.remove(here::here("tables", "tables.tex")))

sys_opp <- readRDS(here::here("data-generated", "encounter-compare-systematic-opportunistic.rds"))
clean_encounter <- readRDS(here::here("data-generated", "encounter-spp-table-systematic-years.rds"))

# Table 1 - sample sizes and infection rates
# ------------------------------------------
asterisk_spp <- sys_opp |>
  filter(category == "positives all time") |>
  mutate(Species = str_to_title(species)) |>
  pull(Species)


# HTML version
# clean_encounter |>
#   mutate(
#     Species = str_to_title(species),
#     n = `0` + `1`,
#     p = round(100 * (`1` / n), 1),
#     p_text = ifelse(`0` == 0, "-", as.character(p))
#   ) |>
#   mutate(Species = ifelse(Species %in% asterisk_spp, paste0(Species, "*"), Species)) |>
#   arrange(desc(p), desc(`1`), desc(`0`)) |>
#   select(Species, Uninfected = `0`, Infected = `1`, `Total n` = n, `Infection rate (%)` = p_text) |>
#   flextable() |>
#   # bg(i = ~ Infected == 0, bg = "grey90") |>
#   bg(i = ~ !grepl("\\*", Species) & Infected == 0, bg = "grey90") |>
#   align(j = "Infection rate (%)", align = "right", part = "all") |>
#   set_table_properties(layout = "autofit") |>
#   set_caption("Table 1. Sample sizes and infection detections for species examined for Sarcotaces sp. infection in 2019--2022. Species highlighted in grey had no observed infections in any specimens examined. Species marked with an asterisk have been found with Sarcotaces sp. cysts in non-systematic sampling events.")

# Add a column to flag rows to be highlighted
table1_df <- clean_encounter |>
  mutate(
      Species = str_to_title(species),
      Species = ifelse(Species %in% asterisk_spp, paste0(Species, "*"), Species),
      n = `0` + `1`,
      p = round(100 * (`1` / n), 1),
      p_text = ifelse(`0` == 0, "-", as.character(p)),
      # Format big numbers with commas
      Uninfected = comma(`0`),
      Infected = comma(`1`),
      `Total n` = comma(n)
    ) |>
  arrange(desc(p), desc(`1`), desc(`0`)) |>
  select(Species, Uninfected, Infected, `Total n`, `Infection rate (%)` = p_text) |>
  mutate(row_highlight = ifelse(!grepl("\\*", Species) & Infected == "0", TRUE, FALSE))

# Write to LaTeX using kableExtra
latex_table1 <- kable(
  table1_df %>% select(-row_highlight),
  format = "latex",
  booktabs = TRUE,
  longtable = FALSE,
  escape = FALSE,
  col.names = c("Species", "Uninfected", "Infected", "Total n", "Infection rate (\\%)"),
  align = c("l", "r", "r", "r", "r")
) %>%
  kable_styling(latex_options = c("hold_position"), position = "center")

# Remove \addlinespace lines
latex_lines <- strsplit(latex_table1, '\n')[[1]]
latex_lines <- latex_lines[!grepl('^\\\\addlinespace', latex_lines)]

# Replace \begin{table}[!h] with \begin{table}[h]
latex_lines[grepl('^\\\\begin\\{table\\}\\[!h\\]', latex_lines)] <- '\\begin{table}[h]'

# Remove \centering
latex_lines <- latex_lines[!grepl('^\\\\centering', latex_lines)]

# Insert caption and label after \begin{table}[h]
caption_line <- '\\caption{Sample sizes and infection detections for species examined for \\sarc infection in 2019--2022. Species highlighted in grey had no observed infections in any specimens examined. Species marked with an asterisk have been found with \\sarc cysts in non-systematic sampling events.}'
label_line <- '\\label{tab:sample-sizes}'
table_start <- which(grepl('^\\\\begin\\{table\\}\\[h\\]', latex_lines))
latex_lines <- append(latex_lines, values = c(caption_line, label_line), after = table_start)

# Insert \begin{center} before \begin{tabular}
tabular_start <- which(grepl('^\\\\begin\\{tabular\\}', latex_lines))
latex_lines <- c(
  latex_lines[1:(tabular_start - 1)],
  '\\begin{center}',
  latex_lines[tabular_start:length(latex_lines)]
)

# Insert \end{center} after \end{tabular}
tabular_end <- which(grepl('^\\\\end\\{tabular\\}', latex_lines))
latex_lines <- c(
  latex_lines[1:tabular_end],
  '\\end{center}',
  latex_lines[(tabular_end + 1):length(latex_lines)]
)

# Row highlighting: Insert \rowcolor{lightgray} before appropriate rows
midrule <- which(grepl('^\\\\midrule', latex_lines))
bottomrule <- which(grepl('^\\\\bottomrule', latex_lines))
data_rows <- (midrule + 1):(bottomrule - 1)
highlight_indices <- which(table1_df$row_highlight)

row_counter <- 0
for (i in data_rows) {
  # Only count actual data rows (not blank lines, etc.)
  if (grepl('\\\\\\\\$', latex_lines[i])) {
    row_counter <- row_counter + 1
    if (row_counter %in% highlight_indices) {
      latex_lines[i] <- paste0('\\rowcolor{lightgray}\n', latex_lines[i])
    }
  }
}

latex_table_colored <- paste(latex_lines, collapse = '\n')
# latex_table_colored <- paste0('% Add this to your LaTeX preamble: \\usepackage{colortbl} \\definecolor{lightgray}{gray}{0.9}\n', latex_table_colored)
write_lines(latex_table_colored, here::here("tables", "tables.tex"))

# Table S1 - cyst frequency
# ------------------------------------------
d <- readRDS(here::here("data-generated", "clean-data.rds"))

cyst_freq <- d |>
  drop_na(sarc_count) |>
  janitor::tabyl(sarc_count) |>
  mutate(n = comma(n)) |>
  select(-percent) |>
  pivot_wider(names_from = sarc_count, values_from = n) |>
  mutate(`# cysts` = "N") |>
  relocate(`# cysts`)

# Convert to LaTeX format
latex_cyst_table <- kable(
  cyst_freq,
  format = "latex",
  booktabs = TRUE,
  longtable = FALSE,
  escape = FALSE,
  col.names = c("\\# cysts", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"),
  align = c("c", rep("r", 10))
) %>%
  kable_styling(latex_options = c("hold_position"), position = "center")

# Replace \begin{table}[!h] with \begin{table}[h]
latex_lines <- strsplit(latex_cyst_table, '\n')[[1]]
latex_lines[grepl('^\\\\begin\\{table\\}\\[!h\\]', latex_lines)] <- '\\begin{table}[h]'

# Remove \centering
latex_lines <- latex_lines[!grepl('^\\\\centering', latex_lines)]

# Insert caption and label after \begin{table}[h]
caption_line <- '\\caption{Frequency of the number of \\sarc cysts found within each individual combined across species}'
label_line <- '\\label{tab:cyst-frequency}'
table_start <- which(grepl('^\\\\begin\\{table\\}\\[h\\]', latex_lines))
latex_lines <- append(latex_lines, values = c(caption_line, label_line), after = table_start)

# Insert \begin{center} before \begin{tabular}
tabular_start <- which(grepl('^\\\\begin\\{tabular\\}', latex_lines))
latex_lines <- c(
  latex_lines[1:(tabular_start - 1)],
  '\\begin{center}',
  latex_lines[tabular_start:length(latex_lines)]
)

# Insert \end{center} after \end{tabular}
tabular_end <- which(grepl('^\\\\end\\{tabular\\}', latex_lines))
latex_lines <- c(
  latex_lines[1:tabular_end],
  '\\end{center}',
  latex_lines[(tabular_end + 1):length(latex_lines)]
)

# Create the final LaTeX table
latex_cyst_table_final <- paste(latex_lines, collapse = '\n')

# Append to existing tables.tex
existing_content <- read_lines(here::here("tables", "tables.tex"))
new_content <- c(existing_content, "", latex_cyst_table_final)
write_lines(new_content, here::here("tables", "tables.tex"))
