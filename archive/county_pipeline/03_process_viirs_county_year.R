# R/03_process_viirs_county_year.R
# Purpose:
#   Process annual VIIRS nighttime lights (EOG VNL) GeoTIFFs and summarize
#   mean radiance to California county × year.
#
# Inputs:
#   - data/processed/boundaries/ca_counties.gpkg     (from R/02_get_ca_counties.R)
#   - data/raw/viirs/                               (manually downloaded GeoTIFFs)
#       * VNL v2.1 files for years 2012–2021
#       * VNL v2.2 files for years 2022–present
#       * Product: annual average radiance ("average_masked")
#       * Units: nW/cm^2/sr
#
# Outputs:
#   - data/processed/viirs/viirs_county_year.csv
#       * One row per county per year (~58 rows per year)
#       * Columns: GEOID, county_name, year, vnl_version, mean_radiance
#
# Notes:
#   - Raw VIIRS files must be uncompressed GeoTIFFs (.tif) and filenames must
#     include a 4-digit year (e.g., 2019) for parsing.
#   - This script crops/masks each raster to California and uses county polygons
#     to extract mean radiance per county.
#   - County boundaries must be created first by running R/02_get_ca_counties.R.
#
# Run script with:
#   source(here::here("R", "03_process_viirs_county_year.R"))

source(here::here("R/00_setup.R"))

library(sf)
library(terra)
library(dplyr)
library(stringr)
library(readr)
library(purrr)
library(tibble)

# ---- Paths ----
COUNTIES_GPKG <- file.path(DIR_BOUND_PROC, "ca_counties.gpkg")
OUT_CSV       <- file.path(DIR_VIIRS_PROC, "viirs_county_year.csv")

# ---- Checks ----
if (!file.exists(COUNTIES_GPKG)) {
  stop("Missing counties file: ", COUNTIES_GPKG, "\nRun R/02_get_ca_counties.R first.")
}

viirs_files <- list.files(
  DIR_VIIRS_RAW,
  pattern = "\\.tif$",
  full.names = TRUE,
  recursive = TRUE
)

if (length(viirs_files) == 0) {
  stop("No VIIRS .tif files found under: ", DIR_VIIRS_RAW)
}


# ---- Load counties ----
ca_counties <- st_read(COUNTIES_GPKG, quiet = TRUE) %>% #creates sf data frame
  select(GEOID, NAME, geom) #keep only these attributes

ca_vect <- terra::vect(ca_counties) #Converts an sf object into a terra vector object (great for raster extraction)

# ---- Helper: summarize one file (for later use)----
summarize_one_file <- function(f) {
  # year from filename
  yr <- str_extract(basename(f), "(19|20)\\d{2}") %>% as.integer()
  if (is.na(yr)) stop("Could not extract year from filename: ", basename(f))
  
  # version from folder name (expects .../viirs/v21/... or .../viirs/v22/...)
  vnl_version <- case_when(
    str_detect(f, "[/\\\\]v21[/\\\\]") ~ "v2.1",
    str_detect(f, "[/\\\\]v22[/\\\\]") ~ "v2.2",
    TRUE ~ NA_character_
  )
  
  if (is.na(vnl_version)) {
    warning("Could not determine VNL version from path: ", f)
  }
  
  r <- rast(f) #reads the GeoTIFF file f into a SpatRaster object
  #r is global raster
  
  # Reproject counties to raster Coordinate Ref System (CRS) for correct extraction
  ca_in_r_crs <- project(ca_vect, crs(r))
  
  # Crop to rectangular extent that includes CA, then sets raster cells outside CA counties to NA
  r_ca <- crop(r, ca_in_r_crs) |> mask(ca_in_r_crs) #r_ca is raster of avg brightness values/pixels for just CA
  
  # If multiple layers exist, take the first (most files here are single-layer)
  if (nlyr(r_ca) > 1) r_ca <- r_ca[[1]]
  
  #Note: r_ca is RASTER of CA brightness values, ca_in_r_crs is SPATVECTOR of CA county polygons
  
  #Extract mean radiance by county
  ex <- terra::extract(r_ca, ca_in_r_crs, fun = mean, na.rm = TRUE) #ex is a df
  
  # ex[[1]] is polygon/county ID (row index); ex[[2]] is mean
  #now 58x5 table. Each row includes 
  tibble(
    GEOID = ca_counties$GEOID[ex[[1]]], #county ID
    county_name = ca_counties$NAME[ex[[1]]], #county name
    year = yr, #year of VNL data
    vnl_version = vnl_version, #v1 or v2
    mean_radiance = ex[[2]] #mean radiance across all pixels in the county
  )
}

# ---- Run all files ----
viirs_county_year <- purrr::map_dfr(viirs_files, summarize_one_file) %>% #(instead of a for loop,) 
  #map_dfr used in place of a for loop to iterate over each raster file f, 
  #applies helper function summarize_one_file() 
  #then stacks tibbles row wise into a df
  #For N files ,viirs_county_year dimensions= (58 × N) × 5
  arrange(year, GEOID) #sorts rows by year and then increasing ID

# ---- Save as CSV in data/processed/viirs/ ----
readr::write_csv(viirs_county_year, OUT_CSV)
message("Saved: ", OUT_CSV)

# Expected output:
# - One row per county per year (~58 rows per year)
# - mean_radiance >= 0 for masked products
# - vnl_version = v2.1 (2012–2021), v2.2 (2022+)
