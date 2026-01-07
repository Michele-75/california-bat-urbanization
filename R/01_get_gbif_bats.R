# R/01_get_gbif_bats.R
# Purpose:
#   Download raw GBIF occurrence records for focal bat species
#   using the GBIF Occurrence Download API. Restrict records to
#   observations with coordinates and years >= 2012, and save
#   raw data and provenance metadata for reproducibility.
#
# Inputs:
#   - GBIF Occurrence Database (via rgbif API)
#   - GBIF Backbone Taxonomy (species name → speciesKey resolution)
#   - User GBIF credentials stored in .Renviron
#
# Outputs:
#   - data/raw/gbif/bats_raw_2012plus.csv        (raw occurrence records)
#   - data/raw/gbif/gbif_backbone_taxa.csv      (species → GBIF key mapping)
#   - data/raw/gbif/gbif_download_meta.csv      (download metadata + filters)
#
# Notes:
#   - Raw data are saved without spatial filtering to California;
#     spatial subsetting and cleaning occur in R/04_clean_gbif_to_ca.R.
#   - Species taxonomy is resolved via the GBIF Backbone to ensure
#     reproducible and documented taxon concepts.
#
# Run script with:
#   source(here::here("R", "01_get_gbif_bats.R"))


source(here::here("R/00_setup.R"))

library(rgbif)
library(janitor)

# ---- Parameters ----
BAT_SPECIES <- c(
  "Myotis yumanensis",
  "Myotis californicus",
  "Lasiurus cinereus"
)

OUT_CSV  <- file.path(DIR_GBIF_RAW, "bats_raw_2012plus.csv")
OUT_META <- file.path(DIR_GBIF_RAW, "gbif_download_meta.csv")

# ---- Early exit if already done ----
if (file.exists(OUT_CSV)) {
  message("Raw GBIF CSV already exists: ", OUT_CSV)
  message("Skipping download. Delete the file if you want to re-run acquisition.")
  quit(save = "no", status = 0)
}

# ---- Taxon keys ----

# GBIF (Global Biodiversity Information Facility) is an international database of biodiversity records.
# Each species in GBIF is assigned a unique numeric identifier: "speciesKey",
# which we must use when requesting occurrence data through the GBIF API.

#Note that taxonomy/backbone can map synonyms or reclassified genera
bat_taxa <- purrr::map_df(BAT_SPECIES, ~ {
  rgbif::name_backbone(name = .x) |>
    as_tibble() |>
    select(scientificName, speciesKey, rank, status)
})

# Save the GBIF backbone resolution (species names -> GBIF keys) for provenance.
# This documents exactly which taxonomic concepts were used to request occurrences.
OUT_TAXA <- file.path(DIR_GBIF_RAW, "gbif_backbone_taxa.csv")
readr::write_csv(bat_taxa, OUT_TAXA)

bat_keys <- bat_taxa$speciesKey

# ---- Credentials (from .Renviron) ----
gbif_user  <- Sys.getenv("GBIF_USER")
gbif_pwd   <- Sys.getenv("GBIF_PWD")
gbif_email <- Sys.getenv("GBIF_EMAIL")

if (gbif_user == "" || gbif_pwd == "" || gbif_email == "") {
  stop("Missing GBIF credentials. Set GBIF_USER, GBIF_PWD, GBIF_EMAIL in .Renviron.")
}

# ---- Submit download ----
#Submits a request to GBIF’s Occurrence Download API for a custom dataset that matches our filters
dl <- occ_download(
  pred_in("taxonKey", bat_keys), #look up with species key
  pred("hasCoordinate", TRUE), #exclude rows with missing coordinates
  pred_gte("year", 2012L), #only include data from 2012 onward
  format = "SIMPLE_CSV",
  user   = gbif_user,
  pwd    = gbif_pwd,
  email  = gbif_email
)

# Wait + fetch
occ_download_wait(dl)

gbif_zip <- occ_download_get(dl, path = DIR_GBIF_RAW, overwrite = TRUE)

# Import- This can take up to 5 minutes.
bats_raw <- occ_download_import(gbif_zip) |>
  as_tibble() |>
  janitor::clean_names() #convert column names to snake_case & other tidying 

# Save raw CSV in data/raw/gbif
readr::write_csv(bats_raw, OUT_CSV)

# ---- Save provenance metadata ----
download_key <- if (is.list(dl) && "key" %in% names(dl)) {
  dl$key
} else {
  as.character(dl)
}

meta <- tibble::tibble(
  download_key = download_key,
  created      = as.character(Sys.time()),
  n_records    = nrow(bats_raw),
  species      = paste(BAT_SPECIES, collapse = "; "),
  filters      = "hasCoordinate==TRUE; year>=2012; taxonKey in focal species"
)

readr::write_csv(meta, OUT_META)

message("Saved raw GBIF data to: ", OUT_CSV)
message("Saved download metadata to: ", OUT_META)
