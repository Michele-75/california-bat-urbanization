# R/01_get_clean_bat_points.R
# Reproducible GBIF acquisition + minimal cleaning for CA bat points (2012–2024)
#
# Outputs:
#   data/raw/gbif/gbif_bats_raw_2012_2024.csv
#   data/raw/gbif/gbif_download_meta.csv
#   data/raw/gbif/gbif_backbone_taxa.csv
#   data/processed/gbif/gbif_bats_points_clean_2012_2024.gpkg
#   data/processed/gbif/gbif_bats_points_clean_qa.csv
#
# Requirements:
#   GBIF credentials set in .Renviron:
#     GBIF_USER, GBIF_PWD, GBIF_EMAIL

source(here::here("R", "00_setup.R"))

library(rgbif)
library(dplyr)
library(readr)
library(tibble)
library(purrr)
library(janitor)
library(sf)

# ---- Parameters ----
YEAR_MIN <- 2012L
YEAR_MAX <- 2024L

# Use names you care about, but rely on GBIF backbone to resolve keys.
# Include both common taxonomies for hoary bat to avoid surprises:
FOCAL_NAMES <- c(
  "Myotis yumanensis",
  "Myotis californicus",
  "Aeorestes cinereus",  # many sources use Aeorestes
  "Lasiurus cinereus"    # many datasets still use Lasiurus
)

# Broad CA-ish bounding box (guardrail only; definitive spatial filtering happens later)
CA_BBOX <- list(lon_min = -125, lon_max = -113, lat_min = 32, lat_max = 42)

# ---- Output paths ----
OUT_RAW_CSV   <- file.path(DIR_GBIF_RAW,  "gbif_bats_raw_2012_2024.csv")
OUT_META_CSV  <- file.path(DIR_GBIF_RAW,  "gbif_download_meta.csv")
OUT_TAXA_CSV  <- file.path(DIR_GBIF_RAW,  "gbif_backbone_taxa.csv")

OUT_CLEAN_GPKG <- file.path(DIR_GBIF_PROC, "gbif_bats_points_clean_2012_2024.gpkg")
OUT_QA_CSV     <- file.path(DIR_GBIF_PROC, "gbif_bats_points_clean_qa.csv")

# ---- Helper: backbone resolution ----
resolve_backbone <- function(name) {
  bb <- rgbif::name_backbone(name = name)
  tibble(
    query_name = name,
    scientificName = bb$scientificName %||% NA_character_,
    canonicalName  = bb$canonicalName %||% NA_character_,
    speciesKey     = bb$speciesKey %||% NA_integer_,
    usageKey       = bb$usageKey %||% NA_integer_,
    rank           = bb$rank %||% NA_character_,
    status         = bb$status %||% NA_character_
  )
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

# ---- Resolve taxa + write provenance ----
taxa_tbl <- map_df(FOCAL_NAMES, resolve_backbone)

if (all(is.na(taxa_tbl$speciesKey))) {
  stop("GBIF backbone lookup failed for all species names. Check spelling or GBIF service.")
}

# Deduplicate keys (Aeorestes/Lasiurus may resolve to same speciesKey depending on backbone)
taxa_keys <- taxa_tbl %>%
  filter(!is.na(speciesKey)) %>%
  distinct(speciesKey) %>%
  pull(speciesKey)

write_csv(taxa_tbl, OUT_TAXA_CSV)

message("Resolved GBIF species keys: ", paste(taxa_keys, collapse = ", "))

# ---- Download from GBIF if raw file not present ----
if (!file.exists(OUT_RAW_CSV)) {
  
  gbif_user  <- Sys.getenv("GBIF_USER")
  gbif_pwd   <- Sys.getenv("GBIF_PWD")
  gbif_email <- Sys.getenv("GBIF_EMAIL")
  
  if (gbif_user == "" || gbif_pwd == "" || gbif_email == "") {
    stop("Missing GBIF credentials. Set GBIF_USER, GBIF_PWD, GBIF_EMAIL in your .Renviron.")
  }
  
  # Submit download request
  dl <- occ_download(
    pred_in("taxonKey", taxa_keys),
    pred("hasCoordinate", TRUE),
    pred_gte("year", YEAR_MIN),
    pred_lte("year", YEAR_MAX),
    format = "SIMPLE_CSV",
    user   = gbif_user,
    pwd    = gbif_pwd,
    email  = gbif_email
  )
  
  occ_download_wait(dl)
  
  zip_path <- occ_download_get(dl, path = DIR_GBIF_RAW, overwrite = TRUE)
  
  raw_tbl <- occ_download_import(zip_path) %>%
    as_tibble() %>%
    janitor::clean_names()
  
  write_csv(raw_tbl, OUT_RAW_CSV)
  
  # Robustly extract download key across rgbif return types
  download_key <- if (is.list(dl) && "key" %in% names(dl)) {
    dl$key
  } else {
    as.character(dl)
  }
  
  meta <- tibble(
    download_key = download_key,
    created_utc  = format(Sys.time(), tz = "UTC"),
    n_records    = nrow(raw_tbl),
    year_min     = YEAR_MIN,
    year_max     = YEAR_MAX,
    focal_names  = paste(FOCAL_NAMES, collapse = "; "),
    taxon_keys   = paste(taxa_keys, collapse = "; "),
    filters      = "hasCoordinate==TRUE; year bounded; taxonKey in focal list"
  )
  
  write_csv(meta, OUT_META_CSV)
  
  message("Downloaded and saved raw GBIF CSV: ", OUT_RAW_CSV)
  
} else {
  message("Raw GBIF CSV already exists; skipping download: ", OUT_RAW_CSV)
}

# ---- Minimal cleaning for downstream spatial workflow ----
raw <- read_csv(OUT_RAW_CSV, show_col_types = FALSE)
n_raw <- nrow(raw)

# Keep a compact, consistent set of columns
dat <- raw %>%
  transmute(
    gbif_id = as.character(gbif_id),
    dataset_key = as.character(dataset_key),
    occurrence_id = as.character(occurrence_id),
    
    species = as.character(species),
    scientific_name = as.character(scientific_name),
    
    year = suppressWarnings(as.integer(year)),
    month = suppressWarnings(as.integer(month)),
    day = suppressWarnings(as.integer(day)),
    
    basis_of_record = as.character(basis_of_record),
    
    lon = suppressWarnings(as.numeric(decimal_longitude)),
    lat = suppressWarnings(as.numeric(decimal_latitude)),
    coordinate_uncertainty_m = suppressWarnings(as.numeric(coordinate_uncertainty_in_meters)),
    
    issue = as.character(issue)
  )

# Basic filters: valid coords, year bounds, rough CA bbox (guardrail)
dat2 <- dat %>%
  filter(
    !is.na(lon), !is.na(lat),
    lon != 0, lat != 0,
    lon >= -180, lon <= 180,
    lat >=  -90, lat <=  90,
    !is.na(year), year >= YEAR_MIN, year <= YEAR_MAX,
    lon >= CA_BBOX$lon_min, lon <= CA_BBOX$lon_max,
    lat >= CA_BBOX$lat_min, lat <= CA_BBOX$lat_max
  ) %>%
  filter(!is.na(species) & species != "")

# Convert to sf points (WGS84). We'll project later as needed.
pts <- st_as_sf(dat2, coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# ---- QA summary ----
qa <- tibble(
  metric = c(
    "raw_rows",
    "rows_after_basic_filters"
  ),
  n = c(
    n_raw,
    nrow(pts)
  )
)

qa_species <- pts %>%
  st_drop_geometry() %>%
  count(species, sort = TRUE) %>%
  transmute(metric = paste0("postfilter_species: ", species), n = n)

qa_year <- pts %>%
  st_drop_geometry() %>%
  count(year, sort = TRUE) %>%
  transmute(metric = paste0("postfilter_year: ", year), n = n)

qa_out <- bind_rows(qa, qa_species, qa_year)
write_csv(qa_out, OUT_QA_CSV)

# ---- Write cleaned points ----
if (file.exists(OUT_CLEAN_GPKG)) file.remove(OUT_CLEAN_GPKG)
st_write(pts, OUT_CLEAN_GPKG, quiet = TRUE)

message("Saved cleaned GBIF bat points: ", OUT_CLEAN_GPKG)
message("Saved QA summary: ", OUT_QA_CSV)
message("Saved backbone taxa provenance: ", OUT_TAXA_CSV)