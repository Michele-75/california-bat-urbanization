# California Bats & Light Pollution

This repository contains my Environmental Data Science capstone project exploring relationships between nighttime light intensity (VIIRS) and nocturnal bat activity across California.


## Reproducible Environment

This project uses `renv` to manage package versions and ensure reproducibility.

To recreate the environment:

```r
renv::restore()


## Reproducibility

This project uses `renv` for dependency management.

To reproduce the analysis:

1. Clone the repository
2. Open the R project
3. Run `renv::restore()` to install required package versions
4. Execute scripts in `R/` to generate processed data
5. Render `analysis/01_master_analysis.qmd`