# R/03_build_grid_and_domain.R
# Purpose:
#   Build a 10km grid over California and define an accessible modeling domain
#   using a UNION of buffered GBIF bat points (instead of convex hull),
#   buffered by BUFFER_M.
#
# Inputs:
#   data/processed/gbif/gbif_bats_points_clean_2012_2024.gpkg
#   data/processed/boundaries/ca_boundary.gpkg
#
# Outputs:
#   data/processed/grid/ca_grid10km.gpkg
#   data/processed/accessibility/accessible_domain_pointbuffer.gpkg
#   data/processed/grid/ca_grid10km_accessible.gpkg

source(here::here("R", "00_setup.R"))

library(sf)
library(dplyr)

# ---- Parameters ----
GRID_SIZE_M <- 10000   # 10 km grid
BUFFER_M    <- 50000  # 100 km point buffer (you can reduce to 50k/75k if desired)

CA_EPSG <- 3310        # NAD83 / California Albers (meters)

# ---- Paths ----
IN_PTS  <- file.path(DIR_GBIF_PROC,  "gbif_bats_points_clean_2012_2024.gpkg")
IN_CA   <- file.path(DIR_BOUND_PROC, "ca_boundary.gpkg")

OUT_GRID_ALL <- file.path(DIR_GRID_PROC, "ca_grid10km.gpkg")
OUT_DOMAIN   <- file.path(DIR_ACCESS_PROC, "accessible_domain_pointbuffer.gpkg")
OUT_GRID_ACC <- file.path(DIR_GRID_PROC, "ca_grid10km_accessible.gpkg")

# ---- Read inputs ----
pts <- st_read(IN_PTS, quiet = TRUE)
ca  <- st_read(IN_CA, layer = "ca_boundary", quiet = TRUE)

# ---- Project to equal-area CRS ----
pts_ae <- st_make_valid(st_transform(pts, CA_EPSG))
ca_ae  <- st_make_valid(st_transform(ca,  CA_EPSG))

# ---- Build 10 km grid over CA extent, then clip to CA ----
grid <- st_make_grid(
  ca_ae,
  cellsize = GRID_SIZE_M,
  square = TRUE
) %>%
  st_as_sf() %>%
  mutate(cell_id = row_number())

# Intersect to clip cells to CA polygon
grid_ca <- st_intersection(grid, st_geometry(ca_ae)) %>%
  mutate(area_m2 = as.numeric(st_area(.))) %>%
  filter(area_m2 >= 1e6) %>%   # drop tiny slivers (< 1 km^2)
  select(cell_id)

# ---- Accessible domain: union of buffered points ----
# Buffer each point, then union buffers into one (multi)polygon
# This avoids convex hull "filling in" large unsampled triangles.
pt_buffers <- st_buffer(st_geometry(pts_ae), dist = BUFFER_M)

domain_geom <- st_union(pt_buffers)

domain <- st_as_sf(domain_geom) %>%
  mutate(domain = paste0("point_buffer_union_", BUFFER_M/1000, "km")) %>%
  st_make_valid()

# Optional cleanup: limit domain to CA plus the same buffer to avoid huge offshore areas
ca_buffer <- st_buffer(st_geometry(ca_ae), dist = BUFFER_M)
domain <- st_intersection(domain, ca_buffer) %>%
  st_make_valid()

# ---- Apply domain mask to grid ----
grid_accessible <- st_intersection(grid_ca, st_geometry(domain)) %>%
  mutate(area_m2 = as.numeric(st_area(.))) %>%
  filter(area_m2 >= 1e6) %>%  # drop tiny slivers again
  select(cell_id)

# ---- Write outputs ----
if (file.exists(OUT_GRID_ALL)) file.remove(OUT_GRID_ALL)
st_write(grid_ca, OUT_GRID_ALL, layer = "ca_grid10km", quiet = TRUE)

if (file.exists(OUT_DOMAIN)) file.remove(OUT_DOMAIN)
st_write(domain, OUT_DOMAIN, layer = "accessible_domain", quiet = TRUE)

if (file.exists(OUT_GRID_ACC)) file.remove(OUT_GRID_ACC)
st_write(grid_accessible, OUT_GRID_ACC, layer = "ca_grid10km_accessible", quiet = TRUE)

# ---- Console summary ----
message("Wrote full CA grid (10km): ", OUT_GRID_ALL, " | n_cells=", nrow(grid_ca))
message("Wrote accessible domain (point-buffer union): ", OUT_DOMAIN)
message("Wrote accessible grid subset: ", OUT_GRID_ACC, " | n_cells=", nrow(grid_accessible))
message("Share retained: ", round(nrow(grid_accessible) / nrow(grid_ca), 3))