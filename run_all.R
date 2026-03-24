# run_all.R
# Purpose:
#   Execute the full data processing pipeline in order.
#   Run this script once (from the project root) before knitting
#   capstone_analysis.Rmd.
#
# Prerequisites:
#   - GBIF credentials set in .Renviron (GBIF_USER, GBIF_PWD, GBIF_EMAIL)
#   - Raw covariate rasters placed in data/raw/ (see data/README.md for details)
#   - renv restored: renv::restore()
#
# Usage:
#   source("run_all.R")

cat("=== Starting full pipeline ===\n\n")

cat("Step 0: Project setup and directory creation...\n")
source(here::here("R", "00_setup.R"))

cat("Step 1: Download and clean GBIF bat points...\n")
source(here::here("R", "01_get_clean_bat_points.R"))

cat("Step 2: Acquire California state boundary...\n")
source(here::here("R", "02_get_ca_boundary.R"))

cat("Step 3: Build 10 km grid and accessible domain...\n")
source(here::here("R", "03_build_grid_and_domain.R"))

cat("Step 4: Compute grid-level covariates...\n")
source(here::here("R", "04_build_grid_covariates.R"))

cat("Step 5: Build grid-level species presence panel...\n")
source(here::here("R", "05_build_grid_presence.R"))

cat("Step 6: Merge presence and covariates into final dataset...\n")
source(here::here("R", "06_build_grid_model_dataset.R"))

cat("\n=== Pipeline complete ===\n")
cat("You can now knit capstone_analysis.Rmd.\n")