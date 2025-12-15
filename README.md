# California Bats & Light Pollution

This repository contains my Environmental Data Science capstone project exploring relationships between nighttime light intensity (VIIRS) and nocturnal bat activity across California.


## How to run this project

This project uses `renv` to manage package versions and ensure reproducibility.

To recreate the environment:

```r
renv::restore()


## Reproducibility

This project uses `renv` for dependency management.

To reproduce the analysis:

1. Clone the repository and open the `.Rproj` file.
2. Run `renv::restore()` to install required package versions.
3. Execute scripts in the `R/` directory in numerical order to download and process data.
4. Render `analysis/01_master_analysis.qmd` to reproduce figures and results.
```

The analysis notebook assumes all processed data files already exist.

##Obtaining Data

California county boundaries were obtained from the U.S. Census TIGER/Line files (2023 release) using the `tigris` R package and saved in processed form as a GeoPackage.
