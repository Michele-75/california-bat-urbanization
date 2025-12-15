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

## How to run this project

1. Clone the repository and open the `.Rproj` file.
2. Run `renv::restore()` to install required package versions.
3. Execute scripts in the `R/` directory to download and process data.
4. Render `analysis/01_master_analysis.qmd` to reproduce figures and results.

The analysis notebook assumes all processed data files already exist.