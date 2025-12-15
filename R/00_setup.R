# R/00_setup.R

library(tidyverse)
library(here) #for directory creation

set.seed(123) #for reproducibility

# central paths and directories
DIR_RAW        <- here("data", "raw")
DIR_PROCESSED  <- here("data", "processed")

DIR_GBIF_RAW   <- here("data", "raw", "gbif")
DIR_VIIRS_RAW  <- here("data", "raw", "viirs")
DIR_BOUND_RAW  <- here("data", "raw", "boundaries")

DIR_BOUND_PROC <- here("data", "processed", "boundaries")
DIR_GBIF_PROC  <- here("data", "processed", "gbif")
DIR_VIIRS_PROC <- here("data", "processed", "viirs")

dir.create(DIR_RAW,        showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_PROCESSED,  showWarnings = FALSE, recursive = TRUE)

dir.create(DIR_GBIF_RAW,   showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_VIIRS_RAW,  showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_BOUND_RAW,  showWarnings = FALSE, recursive = TRUE)

dir.create(DIR_BOUND_PROC, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_GBIF_PROC,  showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_VIIRS_PROC, showWarnings = FALSE, recursive = TRUE)