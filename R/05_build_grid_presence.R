# R/05_build_grid_presence.R
# Purpose:
#   Join cleaned GBIF bat points to the accessible 10 km grid, then build
#   cell_id × species presence table aggregated across 2012–2024.
#
# Inputs:
#   data/processed/gbif/gbif_bats_points_clean_CA_2012_2024.gpkg
#   data/processed/grid/ca_grid10km_accessible.gpkg
#
# Outputs:
#   data/processed/gbif/grid_presence_10km.csv
#   data/processed/gbif/grid_presence_10km.rds
#   data/processed/gbif/grid_presence_10km_QA.csv
#
# Run:
#   source(here::here("R", "05_build_grid_presence.R"))

source(here::here("R/00_setup.R"))

library(sf)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)

# ---- Parameters ----
CA_EPSG <- 3310

# ---- Paths ----
IN_PTS  <- file.path(DIR_GBIF_PROC, "gbif_bats_points_clean_CA_2012_2024.gpkg")
IN_GRID <- file.path(DIR_GRID_PROC, "ca_grid10km_accessible.gpkg")

OUT_CSV <- file.path(DIR_GBIF_PROC, "grid_presence_10km.csv")
OUT_RDS <- file.path(DIR_GBIF_PROC, "grid_presence_10km.rds")
OUT_QA  <- file.path(DIR_GBIF_PROC, "grid_presence_10km_QA.csv")

# ---- Early exit ----
if (file.exists(OUT_CSV) && file.exists(OUT_QA)) {
  message("Presence dataset already exists: ", OUT_CSV)
  message("QA summary already exists:       ", OUT_QA)
  message("Skipping. Delete outputs to re-run.")
  quit(save = "no", status = 0)
}

# ---- Checks ----
if (!file.exists(IN_PTS))  stop("Missing points file: ", IN_PTS,  "\nRun your cleaning script first.")
if (!file.exists(IN_GRID)) stop("Missing grid file:   ", IN_GRID, "\nRun R/03_build_grid_and_domain.R first.")

# ---- Load inputs ----
pts  <- st_read(IN_PTS,  quiet = TRUE)
grid <- st_read(IN_GRID, quiet = TRUE)

stopifnot("cell_id" %in% names(grid))
stopifnot(all(c("gbif_id", "species", "year") %in% names(pts)))

# ---- CRS harmonization ----
grid <- st_make_valid(st_transform(grid, CA_EPSG))
pts  <- st_make_valid(st_transform(pts,  CA_EPSG))

# ---- Keep only fields we need ----
grid <- grid %>%
  mutate(cell_id = as.character(cell_id)) %>%
  select(cell_id)   # keep only attributes; geometry stays attached automatically

pts <- pts %>%
  mutate(
    gbif_id = as.character(gbif_id),
    species = as.character(species),
    year    = as.integer(year)
  ) %>%
  select(gbif_id, species, year)

# ---- QA: uniqueness checks ----
if (n_distinct(grid$cell_id) != nrow(grid)) {
  stop("Grid cell_id is not unique. Fix upstream grid construction before proceeding.")
}

n_pts <- nrow(pts)
n_gbif_unique <- n_distinct(pts$gbif_id)
if (n_gbif_unique != n_pts) {
  warning(
    "gbif_id is not unique (rows=", n_pts, ", unique ids=", n_gbif_unique, "). ",
    "We will count n_distinct(gbif_id) to prevent double-counting."
  )
}

# ---- Join points -> cells ----
# st_within is strict; if you see unmatched points, switch to st_intersects.
pts_joined <- st_join(
  pts,
  grid %>% select(cell_id),
  join = st_within,
  left = TRUE
)

n_unmatched <- sum(is.na(pts_joined$cell_id))
pts_matched <- pts_joined %>% filter(!is.na(cell_id))

# ---- Species set (should be your 3 focal species) ----
species_levels <- pts_matched %>%
  st_drop_geometry() %>%
  distinct(species) %>%
  arrange(species) %>%
  pull(species)

if (length(species_levels) != 3) {
  warning(
    "Expected 3 focal species but found ", length(species_levels), ":\n- ",
    paste(species_levels, collapse = "\n- ")
  )
}

# ---- Aggregate across full study period: cell_id × species ----
agg <- pts_matched %>%
  st_drop_geometry() %>%
  group_by(cell_id, species) %>%
  summarise(
    n_obs = n_distinct(gbif_id),
    .groups = "drop"
  )

# ---- Full panel (zeros for non-detections) ----
presence <- agg %>%
  tidyr::complete(
    cell_id = grid$cell_id,
    species = species_levels,
    fill = list(n_obs = 0L)
  ) %>%
  mutate(is_present = as.integer(n_obs > 0)) %>%
  arrange(cell_id, species)

# ---- QA summary ----
qa_overall <- tibble(
  metric = c(
    "n_cells_accessible",
    "n_points_input",
    "n_points_unique_gbif_id",
    "n_points_unmatched_to_grid",
    "share_points_unmatched",
    "n_species_observed",
    "panel_rows",
    "panel_expected_rows",
    "total_obs_sum_n_obs",
    "study_year_min",
    "study_year_max"
  ),
  value = c(
    nrow(grid),
    n_pts,
    n_gbif_unique,
    n_unmatched,
    n_unmatched / n_pts,
    length(species_levels),
    nrow(presence),
    nrow(grid) * length(species_levels),
    sum(presence$n_obs),
    min(pts$year, na.rm = TRUE),
    max(pts$year, na.rm = TRUE)
  )
)

qa_by_species <- presence %>%
  group_by(species) %>%
  summarise(
    total_obs = sum(n_obs),
    n_cells_present = sum(is_present),
    share_cells_present = mean(is_present),
    .groups = "drop"
  ) %>%
  transmute(
    metric = paste0("species=", species, " | total_obs / n_cells_present / share_present"),
    value  = paste(total_obs, n_cells_present, round(share_cells_present, 4), sep = " / ")
  )

qa <- bind_rows(
  qa_overall %>% mutate(value = as.character(value)),
  qa_by_species
)

# ---- Write outputs ----
readr::write_csv(presence, OUT_CSV)
saveRDS(presence, OUT_RDS)
readr::write_csv(qa, OUT_QA)

message("Saved presence CSV: ", OUT_CSV)
message("Saved presence RDS: ", OUT_RDS)
message("Saved QA summary:   ", OUT_QA)

message("QA quick check:")
print(qa_overall)