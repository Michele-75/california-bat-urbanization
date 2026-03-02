# R/04_build_grid_covariates.R
# Purpose:
#   Build 10km-grid covariates on the accessible grid:
#     1) mean_radiance, log1p_radiance  (VIIRS annual composites, averaged across years)
#     2) pct_developed                  (NLCD 2019 land cover; Developed* classes)
#     3) pct_protected                  (PAD-US Vector Analysis; GAP 1–3 protected)
#     4) pop_density, log1p_pop_density (Population density raster; mean per cell)
#
# Inputs:
#   - data/processed/grid/ca_grid10km_accessible.gpkg
#   - data/raw/viirs/v21/*.tif
#   - data/raw/viirs/v22/*.tif
#   - data/raw/covariates/Annual_NLCD_LndCov_2019_CU_C1V1.tif
#   - data/raw/covariates/PADUS4_1VectorAnalysis_PADUS_Only.gdb
#   - data/raw/covariates/population_density_2020.tif
#
# Output:
#   - data/processed/covariates_grid/grid_covariates_10km.csv

source(here::here("R", "00_setup.R"))

library(sf)
library(terra)
library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(purrr)
library(exactextractr)

CA_EPSG <- 3310

# ---- Paths ----
GRID_GPKG <- file.path(DIR_GRID_PROC, "ca_grid10km_accessible.gpkg")
OUT_CSV   <- file.path(DIR_COV_GRID_PROC, "grid_covariates_10km.csv")

VIIRS_V21_DIR <- file.path(DIR_VIIRS_RAW, "v21")
VIIRS_V22_DIR <- file.path(DIR_VIIRS_RAW, "v22")

NLCD_TIF <- file.path(DIR_COV_RAW, "Annual_NLCD_LndCov_2019_CU_C1V1.tif")
PADUS_GDB <- file.path(DIR_COV_RAW, "PADUS4_1VectorAnalysis_PADUS_Only.gdb")

POP_TIF <- file.path(DIR_COV_RAW, "population_density_2020.tif")

# ---- Checks ----
if (!file.exists(GRID_GPKG)) stop("Missing grid file: ", GRID_GPKG, "\nRun R/03_build_grid_and_domain.R first.")
if (!dir.exists(VIIRS_V21_DIR) && !dir.exists(VIIRS_V22_DIR)) stop("Missing VIIRS folders in: ", DIR_VIIRS_RAW)
if (!file.exists(NLCD_TIF)) stop("Missing NLCD file: ", NLCD_TIF)
if (!dir.exists(PADUS_GDB)) stop("Missing PAD-US geodatabase directory: ", PADUS_GDB)
if (!file.exists(POP_TIF)) stop("Missing population density raster: ", POP_TIF)

# ---- Load grid ----
grid <- st_read(GRID_GPKG, quiet = TRUE) %>%
  st_make_valid() %>%
  st_transform(CA_EPSG)

stopifnot("cell_id" %in% names(grid))

grid_v <- terra::vect(grid) # terra polygons in CA_EPSG

# ============================================================
# 1) VIIRS mean radiance (avg across years)
# ============================================================

viirs_files <- c(
  list.files(VIIRS_V21_DIR, pattern = "\\.tif$", full.names = TRUE),
  list.files(VIIRS_V22_DIR, pattern = "\\.tif$", full.names = TRUE)
)

if (length(viirs_files) == 0) stop("No VIIRS .tif files found in v21/v22 folders.")

extract_year <- function(x) {
  m <- str_match(basename(x), "(201[2-9]|202[0-4])")[, 2]
  as.integer(m)
}

viirs_index <- tibble(file = viirs_files, year = map_int(viirs_files, extract_year)) %>%
  filter(!is.na(year)) %>%
  distinct(year, .keep_all = TRUE) %>%   # if duplicates exist, keep first per year
  arrange(year)

message("VIIRS years used: ", paste(viirs_index$year, collapse = ", "))

extract_viirs_mean_one <- function(tif_path, year_val) {
  r <- terra::rast(tif_path)
  
  # Crop/mask in WGS84 for speed (VIIRS rasters are typically lon/lat)
  grid_4326 <- st_transform(grid, 4326)
  r_crop <- terra::crop(r, terra::vect(grid_4326)) |> terra::mask(terra::vect(grid_4326))
  
  # Project to CA Albers and extract area-weighted mean
  r_ae <- terra::project(r_crop, paste0("EPSG:", CA_EPSG), method = "bilinear")
  
  means <- exactextractr::exact_extract(r_ae, st_transform(grid, CA_EPSG), "mean")
  
  tibble(cell_id = grid$cell_id, year = year_val, mean_radiance = means)
}

viirs_cell_year <- purrr::pmap_dfr(
  list(viirs_index$file, viirs_index$year),
  ~ extract_viirs_mean_one(..1, ..2)
)

viirs_cov <- viirs_cell_year %>%
  group_by(cell_id) %>%
  summarise(
    viirs_years = sum(!is.na(mean_radiance)),
    mean_radiance = mean(mean_radiance, na.rm = TRUE),
    log1p_radiance = log1p(mean_radiance),
    .groups = "drop"
  )

# ============================================================
# 2) pct_developed from NLCD 2019 (label columns)
#    Reuses your proven "fun=table" extraction pattern.
# ============================================================

NLCD_TIF <- file.path(DIR_COV_RAW, "Annual_NLCD_LndCov_2019_CU_C1V1.tif")
if (!file.exists(NLCD_TIF)) stop("Missing NLCD raster: ", NLCD_TIF)

message("Using NLCD raster: ", NLCD_TIF)

nlcd <- terra::rast(NLCD_TIF)

# Project grid to NLCD CRS for correct categorical extraction
grid_v_nlcd <- terra::vect(st_transform(grid, terra::crs(nlcd)))

# Crop/mask NLCD to grid extent for speed
nlcd_crop <- terra::crop(nlcd, grid_v_nlcd) |> terra::mask(grid_v_nlcd)

# Extract class composition per grid cell as a table of counts (fast & matches your county script)
ex_tab <- terra::extract(nlcd_crop, grid_v_nlcd, fun = table, na.rm = TRUE)

ex_tbl <- tibble::as_tibble(ex_tab) %>%
  mutate(ID = as.integer(ID))

class_cols <- setdiff(names(ex_tbl), "ID")
developed_cols <- class_cols[grepl("^Developed", class_cols)]

if (length(developed_cols) == 0) {
  stop(
    "No NLCD class columns start with 'Developed'.\n",
    "Inspect names(ex_tbl) to confirm class label columns were created."
  )
}

dev_cov <- ex_tbl %>%
  mutate(
    total_cells = rowSums(across(all_of(class_cols)), na.rm = TRUE),
    dev_cells   = rowSums(across(all_of(developed_cols)), na.rm = TRUE),
    pct_developed = if_else(total_cells > 0, dev_cells / total_cells, NA_real_) # proportion 0–1
  ) %>%
  left_join(
    tibble(
      ID = as.integer(seq_len(nrow(grid))),
      cell_id = grid$cell_id
    ),
    by = "ID"
  ) %>%
  select(cell_id, pct_developed)

# ============================================================
# 3) pct_protected from PAD-US (GAP 1–3)
#    Reuses your proven layer + GAP field approach.
# ============================================================

padus_layers <- sf::st_layers(PADUS_GDB)
PADUS_LAYER <- padus_layers$name[[1]]
gap_field <- "GAP_Sts"

message("Using PAD-US layer: ", PADUS_LAYER)
message("Using GAP field:    ", gap_field)

padus <- st_read(PADUS_GDB, layer = PADUS_LAYER, quiet = TRUE) %>%
  st_make_valid() %>%
  st_transform(st_crs(grid))

# Speed-up: keep only polygons that intersect the grid
padus <- padus[st_intersects(padus, grid, sparse = FALSE) %>% apply(1, any), ]

if (!(gap_field %in% names(padus))) {
  stop("Expected GAP field not found: ", gap_field, "\nInspect names(padus) to confirm.")
}

gap_vals <- suppressWarnings(as.integer(padus[[gap_field]]))
padus_prot <- padus[!is.na(gap_vals) & gap_vals %in% c(1, 2, 3), ]

if (nrow(padus_prot) == 0) {
  warning("No PAD-US polygons left after filtering GAP 1–3. pct_protected will be NA.")
  prot_cov <- tibble(cell_id = grid$cell_id, pct_protected = NA_real_)
} else {
  
  inter <- st_intersection(
    grid %>% select(cell_id),
    padus_prot 
  )
  
  inter_area <- inter %>%
    mutate(a = as.numeric(st_area(.))) %>%
    st_drop_geometry() %>%
    group_by(cell_id) %>%
    summarise(protected_area_m2 = sum(a, na.rm = TRUE), .groups = "drop")
  
  cell_area <- grid %>%
    mutate(cell_area_m2 = as.numeric(st_area(.))) %>%
    st_drop_geometry() %>%
    select(cell_id, cell_area_m2)
  
  prot_cov <- cell_area %>%
    left_join(inter_area, by = "cell_id") %>%
    mutate(
      protected_area_m2 = coalesce(protected_area_m2, 0),
      pct_protected = pmin(pmax(protected_area_m2 / cell_area_m2, 0), 1)
    ) %>%
    select(cell_id, pct_protected)
}

# ============================================================
# 4) Population density mean per grid cell
# ============================================================

pop <- terra::rast(POP_TIF)

# Project grid to pop CRS for correct extraction, then crop/mask
grid_v_pop <- terra::vect(st_transform(grid, terra::crs(pop)))

pop_crop <- terra::crop(pop, grid_v_pop) |> terra::mask(grid_v_pop)

# Mean population density per cell (area-weighted mean via exactextractr)
# Reproject pop to CA_EPSG so it's in same CRS as grid for exactextractr
pop_ae <- terra::project(pop_crop, paste0("EPSG:", CA_EPSG), method = "bilinear")

pop_mean <- exactextractr::exact_extract(pop_ae, st_transform(grid, CA_EPSG), "mean")

pop_cov <- tibble(
  cell_id = grid$cell_id,
  pop_density = pop_mean,
  log1p_pop_density = log1p(pop_mean)
)

# ============================================================
# Combine + write
# ============================================================

covars <- tibble(cell_id = grid$cell_id) %>%
  left_join(viirs_cov, by = "cell_id") %>%
  left_join(dev_cov, by = "cell_id") %>%
  left_join(prot_cov, by = "cell_id") %>%
  left_join(pop_cov, by = "cell_id") %>%
  mutate(
    pct_developed = pmin(pmax(pct_developed, 0), 1),
    pct_protected = pmin(pmax(pct_protected, 0), 1)
  ) %>%
  arrange(cell_id)

readr::write_csv(covars, OUT_CSV)

message("Saved grid covariates: ", OUT_CSV)
message("Rows (grid cells): ", nrow(covars))
message("Missingness: radiance=", sum(is.na(covars$mean_radiance)),
        " developed=", sum(is.na(covars$pct_developed)),
        " protected=", sum(is.na(covars$pct_protected)),
        " pop=", sum(is.na(covars$pop_density)))
