# R/ — Data Processing Pipeline

This directory contains the modular R scripts that transform raw data into the analysis-ready dataset used by `capstone_analysis.Rmd`. Scripts are numbered in execution order and should be run sequentially (or via `run_all.R` in the project root).

## Script overview

|  |  |  |  |  |  |  |  |  |
|----|----|----|----|----|----|----|----|----|
|  |  |  |  |  |  |  |  |  |
|  | Script |  | Purpose |  | Key inputs |  | Key outputs |  |
|  |  |  |  |  |  |  |  |  |
|  | `00_setup.R` |  | Load core libraries, define directory paths, set reproducibility options |  | — |  | Directory structure created |  |
|  |  |  |  |  |  |  |  |  |
|  | `01_get_clean_bat_points.R` |  | Download GBIF bat occurrences for 3 focal species (2012–2024), clean coordinates, resolve taxonomy |  | GBIF API (requires credentials) |  | `data/raw/gbif/gbif_bats_raw_2012_2024.csv`, `data/processed/gbif/gbif_bats_points_clean_2012_2024.gpkg` |  |
|  |  |  |  |  |  |  |  |  |
|  | `02_get_ca_boundary.R` |  | Retrieve California state boundary from US Census TIGER/Line |  | `tigris` package |  | `data/processed/boundaries/ca_boundary.gpkg` |  |
|  |  |  |  |  |  |  |  |  |
|  | `03_build_grid_and_domain.R` |  | Build 10 km equal-area grid over California, define accessible modeling domain (50 km point buffer), spatially filter bat points to CA |  | Cleaned bat points, CA boundary |  | `data/processed/grid/ca_grid10km_accessible.gpkg`, `data/processed/gbif/gbif_bats_points_clean_CA_2012_2024.gpkg` |  |
|  |  |  |  |  |  |  |  |  |
|  | `04_build_grid_covariates.R` |  | Compute grid-cell-level covariates: VIIRS radiance, NLCD % developed, PAD-US % protected, GPWv4 population density |  | Accessible grid, raw rasters/vectors |  | `data/processed/covariates_grid/grid_covariates_10km.csv` |  |
|  |  |  |  |  |  |  |  |  |
|  | `05_build_grid_presence.R` |  | Join CA bat points to the accessible grid and build a cell × species presence panel (0/1) aggregated across 2012–2024 |  | CA bat points, accessible grid |  | `data/processed/gbif/grid_presence_10km.csv` |  |
|  |  |  |  |  |  |  |  |  |
|  | `06_build_grid_model_dataset.R` |  | Merge presence panel with grid covariates to produce the final analysis dataset |  | Presence panel, grid covariates |  | `data/processed/analysis_grid/grid_model_dataset_10km.csv` |  |
|  |  |  |  |  |  |  |  |  |

## Running the pipeline

From the project root in R:

```         
source("run_all.R") 
```

Or run scripts individually:

```         
source(here::here("R", "00_setup.R")) source(here::here("R", "01_get_clean_bat_points.R")) # ... etc. 
```

## Prerequisites

-   

-   GBIF credentials in `.Renviron` (for script 01)

-   

-   Raw covariate files placed in `data/raw/` (see `data/README.md`)

-   

-   All R packages installed via `renv::restore()`

-   

## Design notes

-   

-   Every script sources `00_setup.R` first, which loads `tidyverse`, `here`, and defines all directory paths.

-   

-   Scripts use `here::here()` for project-root-relative paths, ensuring they work regardless of working directory.

-   

-   Each script checks for the existence of its inputs and writes outputs with explicit QA summaries.

-   

-   CRS is standardized to California Albers (EPSG:3310) for all spatial operations.

-   
