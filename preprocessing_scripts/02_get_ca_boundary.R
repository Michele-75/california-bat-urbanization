# preprocessing_scripts/02_get_ca_boundary.R
# Purpose: Acquire a California state boundary polygon for grid construction
#   and spatial masking.
#
# Inputs:
#   US Census TIGER/Line (cartographic boundary via tigris)
#
# Outputs:
#   data/processed/boundaries/ca_boundary.gpkg  (layer: ca_boundary)

source(here::here("preprocessing_scripts", "00_setup.R"))

library(sf)
library(tigris)
library(dplyr)

options(tigris_use_cache = TRUE)

OUT_GPKG <- file.path(DIR_BOUND_PROC, "ca_boundary.gpkg")

# ---- Download CA state boundary ----
# Cartographic boundary (cb = TRUE) uses simplified geometry for faster operations.
ca <- tigris::states(cb = TRUE, year = 2022) %>%
  st_as_sf() %>%
  filter(STUSPS == "CA") %>%
  select(STUSPS, NAME) %>%
  st_make_valid()

# ---- Write ----
if (file.exists(OUT_GPKG)) file.remove(OUT_GPKG)

st_write(ca, OUT_GPKG, layer = "ca_boundary", quiet = TRUE)

message("Saved California boundary to: ", OUT_GPKG)