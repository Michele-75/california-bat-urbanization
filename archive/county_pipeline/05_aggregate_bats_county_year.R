# R/05_aggregate_bats_county_year.R
# Purpose:
#   Aggregate cleaned GBIF bat observations to county × year × species.
#   Creates a complete panel (all counties, all years 2012–2024, all focal species),
#   filling missing combinations with 0 counts.
#
# Inputs:
#   - data/processed/gbif/gbif_bats_ca_clean.rds (preferred) OR .gpkg
#   - data/processed/boundaries/ca_counties.gpkg
#
# Outputs:
#   - data/processed/gbif/bats_county_year_species.csv
#   - data/processed/gbif/bats_county_year_species.rds
#   - data/processed/gbif/bats_agg_qa_summary.csv

#Run script with source(here::here("R", "05_aggregate_bats_county_year.R"))

source(here::here("R/00_setup.R"))

library(dplyr)
library(tidyr)
library(readr)
library(sf)
library(tibble)


# ---- Parameters ----
YEARS <- 2012:2024
TAXA_CSV <- file.path(DIR_GBIF_RAW, "gbif_backbone_taxa.csv")

# ---- Paths ----
#inputs
IN_RDS        <- file.path(DIR_GBIF_PROC, "gbif_bats_ca_clean.rds")
IN_GPKG       <- file.path(DIR_GBIF_PROC, "gbif_bats_ca_clean.gpkg")
COUNTIES_GPKG <- file.path(DIR_BOUND_PROC, "ca_counties.gpkg")

#outputs
OUT_CSV <- file.path(DIR_GBIF_PROC, "bats_county_year_species.csv")
OUT_RDS <- file.path(DIR_GBIF_PROC, "bats_county_year_species.rds")
OUT_QA  <- file.path(DIR_GBIF_PROC, "bats_agg_qa_summary.csv")


# ---- Early exit ----
if (file.exists(OUT_CSV) && file.exists(OUT_QA)) {
  message("Aggregated bats dataset already exists: ", OUT_CSV)
  message("QA summary already exists: ", OUT_QA)
  message("Skipping. Delete outputs to re-run.")
  quit(save = "no", status = 0)
}


# ---- Checks ----
if (!file.exists(COUNTIES_GPKG)) {
  stop("Missing counties file: ", COUNTIES_GPKG, "\nRun R/02_get_ca_counties.R first.")
}
if (!file.exists(IN_RDS) && !file.exists(IN_GPKG)) {
  stop(
    "Missing cleaned GBIF output.\nExpected one of:\n - ", IN_RDS, "\n - ", IN_GPKG,
    "\nRun R/04_clean_gbif_to_ca.R first."
  )
}


# ---- Load CA & Make Simple Lookup Table ----
#This step later ensures county-year panel includes all 58 counties, including those with zero bat observations
ca_counties <- st_read(COUNTIES_GPKG, quiet = TRUE) %>%
  st_drop_geometry() %>%
  select(GEOID, county_name = NAME) %>%
  distinct()


# ---- Load cleaned GBIF CA points ----
gbif_ca <- if (file.exists(IN_RDS)) {
  readRDS(IN_RDS)
} else {
  st_read(IN_GPKG, quiet = TRUE)
}


# Build table for aggregation 
gbif_tbl <- gbif_ca %>%
  st_drop_geometry() %>% #drop unneeded geometry column
  # keep only what we need for aggregation
  select(gbif_id, GEOID, county_name, year, species_key, species) %>%
  # enforce scope explicitly (defensive)
  filter(year %in% YEARS)

n_in <- nrow(gbif_tbl) #total number of our 3 bat species observations in CA counties 2012-2024

# ---- Define focal species keys from the cleaned dataset ----
# Assumes gbif_bats_ca_clean contains only your intended 3 taxa after script 01’s download filters.
FOCAL_KEYS <- gbif_tbl %>% distinct(species_key) %>% pull(species_key) %>% sort()

if (length(FOCAL_KEYS) != 3) {
  warning(
    "Expected 3 focal species_key values but found ", length(FOCAL_KEYS), ".\n",
    "This can happen if taxonomy/rank issues introduced extra keys.\n",
    "Inspect with: gbif_tbl %>% count(species_key, species, sort=TRUE)"
  )
}

gbif_tbl <- gbif_tbl %>%
  filter(species_key %in% FOCAL_KEYS) #defensive


# ---- Aggregate: county × year × species_key (unit of analysis) ----
# Use n_distinct(gbif_id) as the primary count (more robust than n() if duplicates appear).
bats_agg <- gbif_tbl %>%
  group_by(GEOID, county_name, year, species_key) %>%
  summarise( #summarize each group of county x year x species with count of unique observations
    n_obs = n_distinct(gbif_id),
    .groups = "drop"
  )


# ---- Attach a stable label per species_key (for plots/tables) ----
key_label <- gbif_tbl %>%
  distinct(species_key, species) %>%
  rename(species_label = species)


# ---- Complete the full panel (fill missing combos with 0) ----
# This is important so later joins to VIIRS (county × year) are clean and
# so "no observations" is represented as 0 rather than missing.
bats_panel <- bats_agg %>%
  select(GEOID, year, species_key, n_obs) %>%
  complete( #this function turns implicit missing values into explicit ones
    GEOID = ca_counties$GEOID, #refer to complete list of all counties
    year = YEARS,
    species_key = FOCAL_KEYS,
    fill = list(n_obs = 0) #explicitly given 0
  ) %>%
  left_join(ca_counties, by = "GEOID") %>%      # add county_name back
  left_join(key_label, by = "species_key") %>%  #attach species labels
  arrange(year, GEOID, species_key) #organized by year first, within year by county, within county by species


# ---- QA summary ----
qa_overall <- tibble(
  metric = c(
    "input_rows_point_level",
    "aggregated_nonzero_rows",
    "panel_rows_complete_grid",
    "n_counties",
    "n_years",
    "n_species_keys",
    "total_obs_sum_n_obs"
  ),
  value = c(
    n_in,
    nrow(bats_agg),
    nrow(bats_panel),
    n_distinct(bats_panel$GEOID),
    n_distinct(bats_panel$year),
    n_distinct(bats_panel$species_key),
    sum(bats_panel$n_obs, na.rm = TRUE)
  )
)

qa_by_key <- bats_panel %>%
  group_by(species_key, species_label) %>%
  summarise(total_obs = sum(n_obs), .groups = "drop") %>%
  transmute(metric = paste0("total_obs_species_key=", species_key, " label=", species_label),
            value = total_obs)

qa_by_year <- bats_panel %>%
  group_by(year) %>%
  summarise(total_obs = sum(n_obs), .groups = "drop") %>%
  transmute(metric = paste0("total_obs_year: ", year), value = total_obs)

qa_out <- bind_rows(qa_overall, qa_by_key, qa_by_year)
readr::write_csv(qa_out, OUT_QA)

# ---- Save outputs ----
readr::write_csv(bats_panel, OUT_CSV)
saveRDS(bats_panel, OUT_RDS)

message("Saved bats county×year×species panel (CSV): ", OUT_CSV)
message("Saved bats county×year×species panel (RDS): ", OUT_RDS)
message("Saved QA summary:                          ", OUT_QA)



