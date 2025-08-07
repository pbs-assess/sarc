# Script to copy and rename figures for journal submission
# Based on figures referenced in main.tex

# Create the submission directory (go up one level from R directory)
submission_dir <- here::here("figures", "figs-submission")
dir.create(submission_dir, showWarnings = FALSE, recursive = TRUE)

# Define the mapping of original figure names to submission names
# Based on the order they appear in main.tex
figure_mapping <- list(
  "map.pdf" = "figure1.pdf",
  "diff-sarc-count-by-species.pdf" = "figure2.pdf",
  "depth-plot.pdf" = "figure3.pdf",
  "length-age-ogives.pdf" = "figure4.pdf",
  "length-p50-by-main-species.pdf" = "figure5.pdf",
  "maturity-ratio-main-spp.pdf" = "figure6.pdf",
  "condition-species-coefs.pdf" = "figure7.pdf"
)

# Source directory (go up one level from R directory)
source_dir <- here::here("sarcotaces-overleaf", "figs")

# Copy and rename each figure
for (original_name in names(figure_mapping)) {
  source_path <- file.path(source_dir, original_name)
  target_path <- file.path(submission_dir, figure_mapping[[original_name]])

  if (file.exists(source_path)) {
    file.copy(source_path, target_path, overwrite = TRUE)
    cat("Copied:", original_name, "->", figure_mapping[[original_name]], "\n")
  } else {
    cat("Warning: Source file not found:", source_path, "\n")
  }
}