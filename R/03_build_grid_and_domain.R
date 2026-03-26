# R/03_build_grid_and_domain.R
# Purpose: Build a regular 10 km equal-area grid for California and define an
#   accessible modeling domain using buffered GBIF bat observation points.
#
# Design choices:
#   - Bat points are filtered to the CA boundary before domain construction.
#   - The grid uses regular (unclipped) 10 km squares to avoid sliver polygons.
#   - The analysis grid includes cells with centroids inside CA, plus any cell
#     containing a bat observation ("presence rescue") to prevent edge-case loss.
#   - The accessible subset is defined by a 50 km buffer around all observations,
#     intersected with a buffered CA boundary.
#
# Inputs:
#   data/processed/gbif/gbif_bats_points_clean_2012_2024.gpkg
#   data/processed/boundaries/ca_boundary.gpkg
#
# Outputs:
#   data/processed/gbif/gbif_bats_points_clean_CA_2012_2024.gpkg
#   data/processed/grid/ca_grid10km.gpkg               (layer: ca_grid10km)
#   data/processed/accessibility/accessible_domain_pointbuffer.gpkg
#   data/processed/grid/ca_grid10km_accessible.gpkg    (layer: ca_grid10km_accessible)

source(here::here("R", "00_setup.R"))

library(sf)
library(dplyr)

# ---- Parameters ----
GRID_SIZE_M <- 10000   # 10 km grid cells
BUFFER_M    <- 50000   # 50 km point buffer for accessible domain
CA_EPSG     <- 3310    # NAD83 / California Albers (meters)

# ---- Paths ----
IN_PTS  <- file.path(DIR_GBIF_PROC,  "gbif_bats_points_clean_2012_2024.gpkg")
IN_CA   <- file.path(DIR_BOUND_PROC, "ca_boundary.gpkg")

OUT_PTS_CA   <- file.path(DIR_GBIF_PROC,   "gbif_bats_points_clean_CA_2012_2024.gpkg")
OUT_GRID_ALL <- file.path(DIR_GRID_PROC,   "ca_grid10km.gpkg")
OUT_DOMAIN   <- file.path(DIR_ACCESS_PROC, "accessible_domain_pointbuffer.gpkg")
OUT_GRID_ACC <- file.path(DIR_GRID_PROC,   "ca_grid10km_accessible.gpkg")

# ---- Read inputs ----
pts <- st_read(IN_PTS, quiet = TRUE)
ca  <- st_read(IN_CA, layer = "ca_boundary", quiet = TRUE)

# ---- Project to equal-area CRS ----
pts_ae <- st_make_valid(st_transform(pts, CA_EPSG))
ca_ae  <- st_make_valid(st_transform(ca,  CA_EPSG))

# ---- 1. Filter points to CA boundary ----
# st_intersects is tolerant on borders; st_within can miss boundary points.
inside_ca <- st_intersects(pts_ae, ca_ae, sparse = FALSE)[, 1]
pts_ca <- pts_ae[inside_ca, ]

message("Points input:      ", nrow(pts_ae))
message("Points inside CA:  ", nrow(pts_ca))
message("Points outside CA: ", sum(!inside_ca))

if (nrow(pts_ca) == 0) {
  stop("No points remain after filtering to CA boundary. Check CRS and boundary.")
}

# Write CA-only points for downstream scripts
if (file.exists(OUT_PTS_CA)) file.remove(OUT_PTS_CA)

st_write(
  pts_ca, OUT_PTS_CA,
  layer = "gbif_bats_points_clean_CA_2012_2024",
  quiet = TRUE
)

message("Wrote CA-only bat points: ", OUT_PTS_CA, " | n_points=", nrow(pts_ca))

# ---- 2. Build regular 10 km grid over CA extent ----
grid_full <- st_make_grid(
  ca_ae,
  cellsize = GRID_SIZE_M,
  square = TRUE
) %>%
  st_as_sf() %>%
  mutate(cell_id = row_number())

# Centroid-in-CA test (geometry-only centroids to avoid sf attribute warning)
grid_cent <- st_centroid(st_geometry(grid_full))
cent_in_ca <- st_within(grid_cent, ca_ae, sparse = FALSE)[, 1]
cell_ids_centroid_ca <- grid_full$cell_id[cent_in_ca]

# Presence rescue: retain any cell containing at least one bat observation
pts_cells <- st_join(
  pts_ca,
  grid_full["cell_id"],
  join = st_intersects,
  left = FALSE
) %>%
  st_drop_geometry() %>%
  distinct(cell_id)

cell_ids_presence <- pts_cells$cell_id

# Combine centroid-in-CA cells and presence-rescued cells
cell_ids_ca <- sort(unique(c(cell_ids_centroid_ca, cell_ids_presence)))

grid_ca <- grid_full %>%
  filter(cell_id %in% cell_ids_ca) %>%
  select(cell_id)

stopifnot(inherits(grid_ca, "sf"))

# ---- 3. Build accessible domain from buffered observation points ----
pt_buffers <- st_buffer(st_geometry(pts_ca), dist = BUFFER_M)
domain_geom <- st_union(pt_buffers)

domain <- st_as_sf(domain_geom) %>%
  mutate(domain = paste0("point_buffer_union_", BUFFER_M / 1000, "km")) %>%
  st_make_valid()

# Clip domain to CA + buffer to prevent large offshore extensions
ca_buffer <- st_buffer(st_geometry(ca_ae), dist = BUFFER_M)
domain <- st_intersection(domain, ca_buffer) %>%
  st_make_valid()

# ---- 4. Define accessible grid subset ----
grid_ca_cent <- st_centroid(st_geometry(grid_ca))
cent_in_domain <- st_within(grid_ca_cent, domain, sparse = FALSE)[, 1]
cell_ids_centroid_domain <- grid_ca$cell_id[cent_in_domain]

# Centroid-in-domain cells + presence-rescued cells
cell_ids_acc <- sort(unique(c(cell_ids_centroid_domain, cell_ids_presence)))

grid_accessible <- grid_ca %>%
  filter(cell_id %in% cell_ids_acc) %>%
  select(cell_id)

stopifnot(inherits(grid_accessible, "sf"))

# ---- 5. QA: confirm no CA-filtered points are lost ----
hits <- st_intersects(pts_ca, grid_accessible)
n_unmatched_ca_pts <- sum(lengths(hits) == 0)

message("CA points unmatched to accessible grid: ", n_unmatched_ca_pts)

if (n_unmatched_ca_pts > 0) {
  warning(
    "Some CA points did not match any accessible grid cell.\n",
    "Inspect with: pts_ca[lengths(st_intersects(pts_ca, grid_accessible))==0, ]"
  )
}

# ---- 6. Write outputs ----
if (file.exists(OUT_GRID_ALL)) file.remove(OUT_GRID_ALL)
st_write(grid_ca, OUT_GRID_ALL, layer = "ca_grid10km", quiet = TRUE)

if (file.exists(OUT_DOMAIN)) file.remove(OUT_DOMAIN)
st_write(domain, OUT_DOMAIN, layer = "accessible_domain", quiet = TRUE)

if (file.exists(OUT_GRID_ACC)) file.remove(OUT_GRID_ACC)
st_write(grid_accessible, OUT_GRID_ACC, layer = "ca_grid10km_accessible", quiet = TRUE)

# ---- Console summary ----
message("Wrote CA analysis grid:      ", OUT_GRID_ALL, " | n_cells=", nrow(grid_ca))
message("Wrote accessible domain:     ", OUT_DOMAIN)
message("Wrote accessible grid:       ", OUT_GRID_ACC, " | n_cells=", nrow(grid_accessible))
message("Share retained (accessible / CA): ", round(nrow(grid_accessible) / nrow(grid_ca), 3))