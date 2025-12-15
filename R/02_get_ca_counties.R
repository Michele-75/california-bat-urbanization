# R/02_get_ca_counties.R
source(here::here("R/00_setup.R"))

library(tigris) #access US county boundaries
library(sf) #spatial data engine, st_*() functions interact with spatial data
library(dplyr)

options(tigris_use_cache = TRUE)

OUT_GPKG <- file.path(DIR_BOUND_PROC, "ca_counties.gpkg")

# If output already exists, don't redo work
if (file.exists(OUT_GPKG)) {
  message("Processed counties already exist: ", OUT_GPKG)
  quit(save = "no", status = 0)
}

# Download US counties (cartographic boundary = smaller & faster)
counties_us <- tigris::counties(cb = TRUE, year = 2023, class = "sf") #sf allows dplyr functions to work

# Filter to California (STATEFP == "06") and keep key fields
ca_counties <- counties_us |>
  filter(STATEFP == "06") |>
  select(GEOID, NAME, STATEFP, geometry) |>
  st_transform(3310)  # Reproject to Coordinate Ref System: NAD83 / California Albers (good statewide analysis CRS)
  #units in meters

# Write to GeoPackage (single file in place of shapefiles, GitHub-friendly)
# Output: processed California county boundaries
# - One row per county (58 total), 4 columns
# - Geometry: county polygons
# - CRS: EPSG:3310 (NAD83 / California Albers)
# - Attributes retained: GEOID, NAME, STATEFP
sf::st_write(ca_counties, OUT_GPKG, delete_dsn = TRUE)


message("Saved processed CA counties to: ", OUT_GPKG)