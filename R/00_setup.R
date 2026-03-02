# R/00_setup.R
# Purpose:
#   Central project configuration.
#   Defines directory structure and reproducibility options.
#   Safe to source at the top of every script.

# ---- Core libraries ----
library(tidyverse)
library(here)

# ---- Reproducibility ----
set.seed(123)
options(
  scipen = 999,              # avoid scientific notation
  dplyr.summarise.inform = FALSE
)

# ---- Base directories ----
DIR_RAW       <- here("data", "raw")
DIR_PROCESSED <- here("data", "processed")

# ---- Raw data directories ----
DIR_GBIF_RAW  <- here("data", "raw", "gbif")
DIR_VIIRS_RAW <- here("data", "raw", "viirs")
DIR_BOUND_RAW <- here("data", "raw", "boundaries")
DIR_COV_RAW   <- here("data", "raw", "covariates")

# ---- Processed data directories ----
DIR_GBIF_PROC     <- here("data", "processed", "gbif")
DIR_BOUND_PROC    <- here("data", "processed", "boundaries")
DIR_GRID_PROC     <- here("data", "processed", "grid")
DIR_ACCESS_PROC   <- here("data", "processed", "accessibility")
DIR_COV_GRID_PROC <- here("data", "processed", "covariates_grid")
DIR_ANALYSIS_PROC <- here("data", "processed", "analysis_grid")

# ---- Create directories (idempotent) ----
dirs <- c(
  DIR_RAW,
  DIR_PROCESSED,
  DIR_GBIF_RAW,
  DIR_VIIRS_RAW,
  DIR_BOUND_RAW,
  DIR_COV_RAW,
  DIR_GBIF_PROC,
  DIR_BOUND_PROC,
  DIR_GRID_PROC,
  DIR_ACCESS_PROC,
  DIR_COV_GRID_PROC,
  DIR_ANALYSIS_PROC
)

walk(dirs, ~ dir.create(.x, showWarnings = FALSE, recursive = TRUE))