# R/02_get_ca_boundary.R
# Purpose:
#   Acquire a California state boundary polygon for grid construction and masking.
#
# Output:
#   data/processed/boundaries/ca_boundary.gpkg  (layer: ca_boundary)
#
# Notes:
#   Uses US Census TIGER/Line (states) via tigris.
#   This is lightweight, reproducible, and widely used.

source(here::here("R", "00_setup.R"))

library(sf)
library(tigris)
library(dplyr)

options(tigris_use_cache = TRUE)

OUT_GPKG <- file.path(DIR_BOUND_PROC, "ca_boundary.gpkg")

# ---- Download CA state boundary ----
# Use cartographic boundary ("cb") for simpler geometry and faster ops.
ca <- tigris::states(cb = TRUE, year = 2022) %>%
  st_as_sf() %>%
  filter(STUSPS == "CA") %>%
  select(STUSPS, NAME) %>%
  st_make_valid()

# ---- Write ----
if (file.exists(OUT_GPKG)) file.remove(OUT_GPKG)

st_write(ca, OUT_GPKG, layer = "ca_boundary", quiet = TRUE)

message("Saved California boundary to: ", OUT_GPKG)