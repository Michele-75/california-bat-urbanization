# preprocessing_scripts/06_build_grid_model_dataset.R
# Purpose: Merge grid-level presence (cell_id × species) with grid covariates
#   to produce the final analysis-ready modeling dataset.
#
# Inputs:
#   data/processed/gbif/grid_presence_10km.csv
#   data/processed/covariates_grid/grid_covariates_10km.csv
#
# Outputs:
#   data/processed/analysis_grid/grid_model_dataset_10km.csv

source(here::here("preprocessing_scripts", "00_setup.R"))

library(dplyr)
library(readr)

# ---- Paths ----
IN_PRES <- file.path(DIR_GBIF_PROC, "grid_presence_10km.csv")
IN_COV  <- file.path(DIR_COV_GRID_PROC, "grid_covariates_10km.csv")

DIR_ANALYSIS_GRID <- here::here("data", "processed", "analysis_grid")
if (!dir.exists(DIR_ANALYSIS_GRID)) dir.create(DIR_ANALYSIS_GRID, recursive = TRUE)

OUT_CSV <- file.path(DIR_ANALYSIS_GRID, "grid_model_dataset_10km.csv")

# ---- Early exit if output already exists ----
if (file.exists(OUT_CSV)) {
  message("Final analysis-grid dataset already exists: ", OUT_CSV)
  message("Delete the file to re-run.")
  quit(save = "no", status = 0)
}

# ---- Read inputs ----
presence <- readr::read_csv(IN_PRES, show_col_types = FALSE) %>%
  mutate(
    cell_id    = as.character(cell_id),
    species    = as.character(species),
    n_obs      = as.integer(n_obs),
    is_present = as.integer(is_present)
  )

covars <- readr::read_csv(IN_COV, show_col_types = FALSE) %>%
  mutate(cell_id = as.character(cell_id))

# ---- Merge presence panel with covariates ----
model_df <- presence %>%
  left_join(covars, by = "cell_id") %>%
  arrange(species, cell_id)

# ---- Write ----
readr::write_csv(model_df, OUT_CSV)

# ---- QA summary ----
qa_summary <- model_df %>%
  summarise(
    n_rows = n(),
    n_cells = n_distinct(cell_id),
    n_species = n_distinct(species),
    total_presences = sum(is_present),
    share_present_overall = mean(is_present),
    total_observations = sum(n_obs)
  )

qa_by_species <- model_df %>%
  group_by(species) %>%
  summarise(
    n_cells = n(),
    n_present = sum(is_present),
    share_present = mean(is_present),
    max_n_obs = max(n_obs),
    .groups = "drop"
  )

message("Saved final modeling dataset: ", OUT_CSV)
message("\nQA - Overall:")
print(qa_summary)
message("\nQA - By Species:")
print(qa_by_species)