# R/06_build_county_covariates.R
# Purpose:
#   Build county-level covariates for California counties:
#     1) pct_developed  (from NLCD land cover; developed classes 21–24)
#     2) pct_protected  (from PAD-US Vector Analysis; GAP 1–3 treated as protected)
#
# Inputs:
#   - data/processed/boundaries/ca_counties.gpkg (from R/02_get_ca_counties.R)
#   - data/raw/covariates/Annual_NLCD_LndCov_2019_CU_C1V1.tif
#   - data/raw/covariates/PADUS4_1VectorAnalysis_PADUS_Only.gdb (directory)
#
# Outputs:
#   - data/processed/covariates/ca_county_covariates.csv
#   - data/processed/covariates/ca_county_covariates.rds
#   - data/processed/covariates/county_covariates_qa.csv
#
# Run:
#   source(here::here("R", "06_build_county_covariates.R"))

source(here::here("R/00_setup.R"))

library(sf)
library(terra)
library(dplyr)
library(readr)
library(stringr)
library(tibble)

# ---- Paths ----
COUNTIES_GPKG <- file.path(DIR_BOUND_PROC, "ca_counties.gpkg")

NLCD_TIF <- file.path(DIR_COV_RAW, "Annual_NLCD_LndCov_2019_CU_C1V1.tif")
PADUS_GDB <- file.path(DIR_COV_RAW, "PADUS4_1VectorAnalysis_PADUS_Only.gdb")

OUT_CSV <- file.path(DIR_COV_PROC, "ca_county_covariates.csv")
OUT_RDS <- file.path(DIR_COV_PROC, "ca_county_covariates.rds")
OUT_QA  <- file.path(DIR_COV_PROC, "county_covariates_qa.csv")

# ---- Early exit ----
if (file.exists(OUT_CSV) && file.exists(OUT_QA)) {
  message("County covariates already exist: ", OUT_CSV)
  message("QA summary already exists:       ", OUT_QA)
  message("Skipping. Delete outputs to re-run.")
  quit(save = "no", status = 0)
}

# ---- Checks ----
if (!file.exists(COUNTIES_GPKG)) stop("Missing counties file: ", COUNTIES_GPKG, "\nRun R/02_get_ca_counties.R first.")
if (!file.exists(NLCD_TIF))      stop("Missing NLCD raster: ", NLCD_TIF)
if (!dir.exists(PADUS_GDB))      stop("Missing PAD-US geodatabase directory: ", PADUS_GDB)

message("Using NLCD raster: ", NLCD_TIF)
message("Using PAD-US GDB:  ", PADUS_GDB)

# ---- Load counties ----
# We keep geometry for intersections, but also preserve GEOID + county_name.
ca_counties <- st_read(COUNTIES_GPKG, quiet = TRUE) %>%
  st_make_valid() %>%
  select(GEOID, county_name = NAME, geom)



# ============================================================
# 1) pct_developed from NLCD (classes 21–24) — FAST VERSION
# ============================================================

# NLCD developed classes value:
#    Developed, Open Space
#    Developed, Low Intensity
#    Developed, Medium Intensity
#    Developed, High Intensity


nlcd <- terra::rast(NLCD_TIF)

# Convert counties to terra vector and project to NLCD CRS for correct raster extraction
ca_vect <- terra::vect(sf::st_transform(ca_counties, terra::crs(nlcd)))

# Crop/mask for speed: keep only California pixels
nlcd_ca <- terra::crop(nlcd, ca_vect) |>
  terra::mask(ca_vect)


# Extract *summarized* land-cover composition per county:
# - one row per county polygon
# - columns are land cover class codes with pixel counts
system.time({
  ex_tab <- terra::extract(nlcd_ca, ca_vect, fun = table, na.rm = TRUE)
})

# ex_tab structure:
#   ID  <classcode1> <classcode2> ... (counts)
# Some class code columns may be missing for some counties; treat as 0.

ex_tbl <- tibble::as_tibble(ex_tab) %>%
  mutate(ID = as.integer(ID))

# All land-cover class columns (everything except ID)
class_cols <- setdiff(names(ex_tbl), "ID")

# Developed columns are those whose label contains "Developed"
# (matches: Developed, Open Space / Low / Medium / High Intensity)
developed_cols <- class_cols[grepl("^Developed", class_cols)]

# Total pixels per county = sum across all class count columns
# Developed pixels per county = sum of columns whose class code is in DEVELOPED_CLASSES
dev_cov <- ex_tbl %>%
  mutate(
    total_cells = rowSums(dplyr::across(dplyr::all_of(class_cols)), na.rm = TRUE),
    dev_cells   = rowSums(dplyr::across(dplyr::all_of(developed_cols)), na.rm = TRUE),
    pct_developed = dplyr::if_else(total_cells > 0, 100 * dev_cells / total_cells, NA_real_)
  ) %>%
  # Map terra polygon ID (1..58) back to GEOID in the same order as ca_counties
  left_join(
    tibble::tibble(
      ID = as.integer(seq_len(nrow(ca_counties))),
      GEOID = ca_counties$GEOID
    ),
    by = "ID"
  ) %>%
  dplyr::select(GEOID, pct_developed)



# ============================================================
# 2) pct_protected from PAD-US (GAP 1–3)
# ============================================================

# PAD-US Vector Analysis is already de-overlapped and prioritized, which is ideal for area summaries.
# We treat GAP Status 1–3 as "protected" and compute the share of county area covered by protected polygons.

padus_layers <- sf::st_layers(PADUS_GDB)

# Usually the first layer is the main polygon layer for the vector analysis file.
PADUS_LAYER <- padus_layers$name[[1]]

padus <- st_read(PADUS_GDB, layer = PADUS_LAYER, quiet = TRUE) %>%
  st_make_valid() %>%
  st_transform(st_crs(ca_counties))

# Keep only PAD-US polygons that intersect CA counties (big speed-up)
padus <- padus[st_intersects(padus, ca_counties, sparse = FALSE) %>% apply(1, any), ]

# Identify the GAP status field (names vary slightly across releases)
gap_field <- "GAP_Sts"
#Use if unsure of column name
#gap_field <- names(padus)[str_detect(names(padus), regex("gap", ignore_case = TRUE))][1] 

if (is.na(gap_field) || gap_field == "") {
  stop(
    "Could not find a GAP status field in PAD-US layer.\n",
    "Inspect fields with: names(padus)\n",
    "Then update this script to point to the correct field."
  )
}

gap_vals <- suppressWarnings(as.integer(padus[[gap_field]]))

# Protected = GAP 1–3
padus_prot <- padus[!is.na(gap_vals) & gap_vals %in% c(1, 2, 3), ]

if (nrow(padus_prot) == 0) {
  warning("No polygons left after filtering GAP 1–3. pct_protected will be NA.")
  prot_cov <- tibble(GEOID = ca_counties$GEOID, pct_protected = NA_real_)
} else {
  # Intersect protected areas with counties and sum intersection areas per county
  inter <- st_intersection(
    ca_counties %>% select(GEOID),
    padus_prot %>% select(Shape) #Or whatever your geometry column is named 
  )
  
  inter_area <- inter %>%
    mutate(a = as.numeric(sf::st_area(.))) %>%
    st_drop_geometry() %>%
    group_by(GEOID) %>%
    summarise(protected_area_m2 = sum(a, na.rm = TRUE), .groups = "drop")
  
  county_area <- ca_counties %>%
    mutate(county_area_m2 = as.numeric(sf::st_area(.))) %>%
    st_drop_geometry() %>%
    select(GEOID, county_area_m2)
  
  prot_cov <- county_area %>%
    left_join(inter_area, by = "GEOID") %>%
    mutate(
      protected_area_m2 = coalesce(protected_area_m2, 0),
      pct_protected = 100 * protected_area_m2 / county_area_m2
    ) %>%
    select(GEOID, pct_protected)
}


# ============================================================
# Combine + QA + Save
# ============================================================

covs <- tibble(GEOID = ca_counties$GEOID) %>%
  left_join(ca_counties %>% st_drop_geometry() %>% select(GEOID, county_name), by = "GEOID") %>%
  left_join(dev_cov, by = "GEOID") %>%
  left_join(prot_cov, by = "GEOID")

qa <- tibble(
  metric = c(
    "nlcd_file",
    "padus_gdb",
    "padus_layer",
    "gap_field_used",
    "n_counties",
    "share_missing_pct_developed",
    "share_missing_pct_protected",
    "pct_developed_min",
    "pct_developed_max",
    "pct_protected_min",
    "pct_protected_max"
  ),
  value = c(
    NLCD_TIF,
    PADUS_GDB,
    PADUS_LAYER,
    gap_field,
    nrow(covs),
    mean(is.na(covs$pct_developed)),
    mean(is.na(covs$pct_protected)),
    suppressWarnings(min(covs$pct_developed, na.rm = TRUE)),
    suppressWarnings(max(covs$pct_developed, na.rm = TRUE)),
    suppressWarnings(min(covs$pct_protected, na.rm = TRUE)),
    suppressWarnings(max(covs$pct_protected, na.rm = TRUE))
  )
)

readr::write_csv(covs, OUT_CSV)
saveRDS(covs, OUT_RDS)
readr::write_csv(qa, OUT_QA)

message("Saved county covariates (CSV): ", OUT_CSV)
message("Saved county covariates (RDS): ", OUT_RDS)
message("Saved QA summary:              ", OUT_QA)
