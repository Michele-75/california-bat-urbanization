# Data Dictionary

## Final analysis dataset: `grid_model_dataset_10km.csv`

**Location:** `data/processed/analysis_grid/grid_model_dataset_10km.csv`

This dataset is a balanced panel of 10 km grid cells × 3 bat species, representing the accessible modeling domain within California. Each row corresponds to one species in one grid cell, with cumulative presence/absence (2012–2024) and grid-level landscape covariates.

**Dimensions:** 8,433 rows × 10 columns (2,811 cells × 3 species)

### Variables

| Column | Type | Units | Description |
|----|----|----|----|
| `cell_id` | character | — | Unique identifier for each 10 km grid cell |
| `species` | character | — | Focal bat species. One of: *Aeorestes cinereus*, *Myotis californicus*, *Myotis yumanensis* |
| `n_obs` | integer | count | Number of unique GBIF occurrence records for this species in this cell, aggregated across 2012–2024. Zero indicates no records. |
| `is_present` | integer (0/1) | — | Binary presence indicator: 1 if `n_obs > 0`, 0 otherwise. Primary response variable for logistic modeling. |
| `viirs_years` | integer | count | Number of annual VIIRS composites available for this cell (out of 13 possible years). Used for QA. |
| `mean_radiance` | double | nW/cm²/sr | Mean nighttime light radiance across available annual VIIRS composites (2012–2024). Higher values indicate brighter, more urbanized areas. |
| `log1p_radiance` | double | log(nW/cm²/sr + 1) | Log-transformed radiance: `log(1 + mean_radiance)`. Reduces right skew for modeling. |
| `pct_developed` | double | proportion (0–1) | Proportion of cell area classified as developed land (NLCD 2019 "Developed" classes: Open Space, Low, Medium, and High Intensity). |
| `pct_protected` | double | proportion (0–1) | Proportion of cell area under GAP Status 1, 2, or 3 protection (PAD-US 4.1). These categories represent lands managed for biodiversity conservation. |
| `pop_density` | double | people/km² | Mean population density within the cell (WorldPop 2020, UN-adjusted). |
| `log1p_pop_density` | double | log(people/km² + 1) | Log-transformed population density: `log(1 + pop_density)`. Primary predictor in final models. |

### Notes

-   **Presence is cumulative across years.** A cell is coded as "present" if *any* observation of that species was recorded in the cell across the entire 2012–2024 study period. This is not a year-by-year panel.
-   **Covariates are spatial, not temporal.** Radiance is averaged across years; land cover and population density are from a single snapshot (NLCD 2019, WorldPop 2020). This reflects the finding that spatial variation dominates temporal variation at this scale.
-   **Grid cells are regular 10 km squares** in California Albers (EPSG:3310), restricted to the accessible modeling domain (50 km buffer around all bat observations, intersected with California).
-   **Zero observations do not necessarily mean absence.** Most grid cells have never been surveyed for bats. The accessible domain restriction mitigates this but does not eliminate it.

### Coordinate reference system

All spatial data use EPSG:3310 (NAD83 / California Albers), an equal-area projection appropriate for California-scale ecological analysis.

### Source files

The final dataset is produced by `R/06_build_grid_model_dataset.R`, which joins: - `data/processed/gbif/grid_presence_10km.csv` (from `R/05_build_grid_presence.R`) - `data/processed/covariates_grid/grid_covariates_10km.csv` (from `R/04_build_grid_covariates.R`)
