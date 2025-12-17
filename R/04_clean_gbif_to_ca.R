# R/04_clean_gbif_to_ca.R
# Purpose:
#   Clean GBIF bat occurrences and spatially filter to California counties.
#   Assign each occurrence to a CA county polygon (county × point).
#
# Inputs:
#   - data/raw/gbif/bats_raw_2012plus.csv          (from R/01_get_gbif_bats.R)
#   - data/processed/boundaries/ca_counties.gpkg   (from R/02_get_ca_counties.R)
#
# Outputs:
#   - data/processed/gbif/gbif_bats_ca_clean.gpkg  (sf points w/ county fields)
#   - data/processed/gbif/gbif_bats_ca_clean.rds   (same, for fast R loads)
#   - data/processed/gbif/gbif_clean_qa_summary.csv

#Run script with source(here::here("R", "04_clean_gbif_to_ca.R"))

source(here::here("R/00_setup.R"))

library(dplyr)
library(readr)
library(sf)
library(tibble)


# ---- Paths ----
#inputs
IN_GBIF_CSV   <- file.path(DIR_GBIF_RAW, "bats_raw_2012plus.csv")
COUNTIES_GPKG <- file.path(DIR_BOUND_PROC, "ca_counties.gpkg")

#outputs
OUT_GPKG <- file.path(DIR_GBIF_PROC, "gbif_bats_ca_clean.gpkg")
OUT_RDS  <- file.path(DIR_GBIF_PROC, "gbif_bats_ca_clean.rds")
OUT_QA   <- file.path(DIR_GBIF_PROC, "gbif_clean_qa_summary.csv")


# ---- Early exit ----
if (file.exists(OUT_GPKG) && file.exists(OUT_QA)) {
  message("Processed GBIF CA dataset already exists: ", OUT_GPKG)
  message("QA summary already exists: ", OUT_QA)
  message("Skipping. Delete outputs to re-run.")
  quit(save = "no", status = 0)
}


# ---- Checks ----
if (!file.exists(IN_GBIF_CSV)) {
  stop("Missing raw GBIF file: ", IN_GBIF_CSV, "\nRun R/01_get_gbif_bats.R first.")
}
if (!file.exists(COUNTIES_GPKG)) {
  stop("Missing counties file: ", COUNTIES_GPKG, "\nRun R/02_get_ca_counties.R first.")
}


# ---- Load data ----
gbif_raw <- readr::read_csv(IN_GBIF_CSV, show_col_types = FALSE)

ca_counties <- sf::st_read(COUNTIES_GPKG, quiet = TRUE) %>%
  st_transform(3310) %>% #change coordinate system for CA-wide analysis
  select(GEOID, county_name = NAME, geom)

n_raw <- nrow(gbif_raw)

# ---- Keep fields (analysis + provenance) ----
gbif_keep <- gbif_raw %>%
  transmute( #selects & creates/renames
    gbif_id       = as.character(gbif_id),
    dataset_key   = as.character(dataset_key),
    occurrence_id = as.character(occurrence_id),
    
    species         = as.character(species),
    scientific_name = as.character(scientific_name),
    
    country_code    = as.character(country_code),
    state_province  = as.character(state_province),
    locality        = as.character(locality),
    
    basis_of_record = as.character(basis_of_record),
    
    year  = suppressWarnings(as.integer(year)),
    month = suppressWarnings(as.integer(month)),
    day   = suppressWarnings(as.integer(day)),
    event_date = as.character(event_date),
    
    lat = suppressWarnings(as.numeric(decimal_latitude)),
    lon = suppressWarnings(as.numeric(decimal_longitude)),
    coord_unc_m = suppressWarnings(as.numeric(coordinate_uncertainty_in_meters)),
    
    issue = as.character(issue)
  )


# ---- Filter bat data based on broad coordinates & year ----
# CA-ish box is intentionally broad; the spatial join is the definitive CA filter.
gbif_coord_ok <- gbif_keep %>%
  filter(
    !is.na(lon), !is.na(lat),
    lon != 0, lat != 0,
    lon >= -180, lon <= 180,
    lat >=  -90, lat <=  90, #physical plausibility
    lon >= -125, lon <= -113,
    lat >=   32, lat <=   42 #Broad CA guardrail
  ) %>%
  filter(!is.na(year) & year >= 2012 & year <= 2024) %>% #Scope of project
  filter(!is.na(species) & species != "")


# ---- Convert to sf (WGS84 -> CA Albers) ----
gbif_sf <- st_as_sf( #turn rows into points stored in new column called "geometry")
  gbif_coord_ok,
  coords = c("lon", "lat"),
  crs = 4326, #these coordinates are in lat/long degrees
  remove = FALSE
) %>%
  st_transform(3310) #reproject geometry, needs to match counties' CRS


# ---- Spatial join to counties (actual CA filter) ----
gbif_ca <- st_join(
  gbif_sf,
  ca_counties, 
  join = st_within, #only keep observations completely inside CA counties
  left = FALSE
)


# ---- QA summary to track dropped data----
qa_overall <- tibble(
  metric = c(
    "raw_rows",
    "after_coord_and_year_filters",
    "after_spatial_join_to_ca_counties"
  ),
  n = c(
    n_raw,
    nrow(gbif_sf),
    nrow(gbif_ca)
  )
)

qa_species <- gbif_ca %>%
  st_drop_geometry() %>% #so behaves like regular data frame 
  count(species, sort = TRUE) %>%
  transmute(metric = paste0("final_species: ", species), n = n)

qa_year <- gbif_ca %>%
  st_drop_geometry() %>%
  count(year, sort = FALSE) %>%
  transmute(metric = paste0("final_year: ", year), n = n)

qa_out <- bind_rows(qa_overall, qa_species, qa_year)
readr::write_csv(qa_out, OUT_QA)


# ---- Save outputs ----
st_write(gbif_ca, OUT_GPKG, delete_dsn = TRUE, quiet = TRUE) #save cleaned CA bat points to .gpkg
saveRDS(gbif_ca, OUT_RDS) #Save in .rds format for fast loading in R

message("Saved processed CA GBIF points (GPKG): ", OUT_GPKG)
message("Saved processed CA GBIF points (RDS):  ", OUT_RDS)
message("Saved QA summary:                     ", OUT_QA)
