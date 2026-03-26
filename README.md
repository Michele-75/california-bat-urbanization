# How Do Human-Modified Landscapes Shape Bat Observations in California?

**Environmental Data Science Capstone Project**\
Michele Perry \| Yale School of the Environment \| 2025–2026

------------------------------------------------------------------------

## Overview

This project examines how human-modified landscapes shape bat observation patterns across California. Using 455 publicly available occurrence records from [GBIF](https://www.gbif.org/) for three focal bat species (*Aeorestes cinereus*, *Myotis yumanensis*, *Myotis californicus*), observations are mapped onto a 10 km equal-area grid and analyzed with logistic regression to assess whether species respond differently to an urban intensity gradient.

### Key findings

1.  **Population density alone** predicts bat presence nearly as well as multi-variable models (cross-validated AUC \~78%). At 10 km resolution, nighttime light and percent developed land are highly correlated with population density and add less than 1% predictive improvement.
2.  **All three species** show positive associations with increasing urban intensity, consistent with observer bias in opportunistic data.
3.  **Species differ in response strength.** *Myotis californicus* shows a notably weaker association with the urban gradient than the other two species, suggesting greater sensitivity to human-modified landscapes.

### A note on observer bias

GBIF data are opportunistic, meaning they reflect where people looked, not just where bats are. Population density captures both ecological signal and reporting effort, so results are interpreted as patterns in *observed* presence rather than confirmed ecological preference. The differences *between* species are more informative than the absolute trends.

------------------------------------------------------------------------

## Repository structure

```         
├── capstone_analysis.Rmd           # Main narrative notebook (source)
├── capstone_analysis.nb.html       # Pre-rendered HTML notebook (download to view)
├── capstone_analysis.md            # GitHub-rendered version (view on GitHub)
├── capstone_analysis_files/        # Figures for the .md rendering
├── preprocessing_scripts/          # Data processing pipeline (scripts 00–06)
│   ├── 00_setup.R
│   ├── 01_get_clean_bat_points.R
│   ├── 02_get_ca_boundary.R
│   ├── 03_build_grid_and_domain.R
│   ├── 04_build_grid_covariates.R
│   ├── 05_build_grid_presence.R
│   ├── 06_build_grid_model_dataset.R
│   ├── run_all.R
│   └── README.md
├── data/
│   ├── raw/                        # (gitignored) Downloaded source data
│   ├── processed/                  # Key processed outputs (see data/README.md)
│   ├── README.md                   # Data acquisition instructions
│   └── data_dictionary.md          # Variable descriptions for final dataset
├── communication/                  # Presentation materials
│   └── lightning_presentation_slide.png
├── renv.lock                       # R package dependencies
├── .gitignore
└── LICENSE
```

------------------------------------------------------------------------

## Reproducing this project

### Quick start (just want to read the analysis?)

**Browse on GitHub:** Open [`capstone_analysis.md`](capstone_analysis.md) to view the rendered markdown with all figures directly in your browser.

**Full interactive version:** Download `capstone_analysis.nb.html` and open it in any browser for the HTML notebook with collapsible code sections.

### Re-run the analysis locally

The repository includes the final processed dataset and key spatial files, so you can re-knit the notebook without running the full pipeline:

1.  Clone the repository
2.  Restore the R environment: `renv::restore()`
3.  Open `capstone_analysis.Rmd` in RStudio and knit

### Full pipeline (from raw data)

To reproduce everything from scratch, including downloading raw data:

1.  Clone the repository:

```         
   git clone https://github.com/Michele-75/your-repo-name.git
   cd your-repo-name
```

2.  Restore the R environment: `renv::restore()`

3.  Set GBIF credentials in `.Renviron`:

    ```         
    GBIF_USER=your_username
    GBIF_PWD=your_password
    GBIF_EMAIL=your_email
    ```

4.  Download raw covariate rasters (see `data/README.md` for sources and file placement)

5.  Run the pipeline: `source(here::here("preprocessing_scripts", "run_all.R"))`

6.  Knit `capstone_analysis.Rmd`

------------------------------------------------------------------------

## Data sources

| Dataset | Source | Resolution |
|-----------------------|--------------------|-----------------------------|
| Bat occurrences | [GBIF](https://www.gbif.org/) | Point records, 2012–2024 |
| Nighttime light | [VIIRS VNL](https://eogdata.mines.edu/products/vnl/) v2.1/v2.2 | \~500 m annual composites |
| Land cover | [NLCD 2019](https://www.mrlc.gov/) | 30 m |
| Protected areas | [PAD-US 4.1](https://www.usgs.gov/programs/gap-analysis-project/science/pad-us-data-download) | Vector polygons |
| Population density | [GPWv4 rev. 11](https://www.earthdata.nasa.gov/data/catalog/sedac-ciesin-sedac-gpwv4-popdens-r11-4.11) (SEDAC/CIESIN) | \~1 km, 2020 |
| State boundary | US Census TIGER/Line | via `tigris` R package |

------------------------------------------------------------------------

## Focal species

| Species | Common name | Ecology |
|---------------------|------------------------------|---------------------|
| *Aeorestes cinereus* | Hoary bat | Long-distance migratory, tree-roosting |
| *Myotis yumanensis* | Yuma myotis | Water-associated, often found near bridges and buildings |
| *Myotis californicus* | California myotis | Crevice-roosting, widespread across the state |

------------------------------------------------------------------------

## Related work

An interactive R Shiny application for exploring California bat occurrence patterns is available here: [California Bat Occurrence Explorer](https://michele-75.shinyapps.io/california-bat-occurrence-explorer/).

## Acknowledgements

This project was completed as part of the [Yale Environmental Data Science Certificate Program](https://environment.yale.edu/certificates/data), 2025–2026 cohort. I am grateful to the program's leaders, professors, and mentors for their guidance and support. Capstone projects from the full cohort can be viewed at the [Cohort 1 Capstones repository](https://github.com/yse-eds-cert/cohort1-capstones). Bat occurrence data were provided by the [Global Biodiversity Information Facility (GBIF)](https://www.gbif.org/).

## License

This project was created for educational purposes as part of a graduate certificate program. All code and original text are provided as-is for learning and reference. Data used in this analysis are publicly available from the sources cited above and are subject to their respective terms of use. GBIF-mediated occurrence data should be cited according to [GBIF citation guidelines](https://www.gbif.org/citation-guidelines).
