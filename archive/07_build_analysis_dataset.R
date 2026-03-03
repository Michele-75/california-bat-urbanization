# R/07_build_analysis_dataset.R
# Purpose:
#   Attach county-level covariates (pct_developed, pct_protected)
#   to the county × year × species bats + VIIRS dataset.
#
# Inputs:
#   - data/processed/analysis/bats_viirs_county_year_species.csv
#       (from R/06_build_analysis_dataset.R)
#   - data/processed/covariates/ca_county_covariates.csv 
#       (from R/06_build_county_covariates.R)
#
# Outputs:
#   - data/processed/analysis/bats_viirs_covars_county_year_species.{csv,rds}

source(here::here("R/00_setup.R"))

library(dplyr)
library(readr)

# ---- paths ----
in_bats <- here::here("data", "processed", "analysis", "bats_viirs_county_year_species.csv")
in_covs <- here::here("data", "processed", "covariates", "ca_county_covariates.csv")

out_dir <- here::here("data", "processed", "analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_csv <- file.path(out_dir, "bats_viirs_covars_county_year_species.csv")
out_rds <- file.path(out_dir, "bats_viirs_covars_county_year_species.rds")

# ---- read ----
bats <- readr::read_csv(in_bats, show_col_types = FALSE) %>%
  mutate(GEOID = as.character(GEOID))

covs <- readr::read_csv(in_covs, show_col_types = FALSE) %>%
  mutate(GEOID = as.character(GEOID)) %>%
  select(any_of(c("GEOID", "pct_developed", "pct_protected")))

# ---- minimal guardrail to prevent join duplication ----
covs <- covs %>% distinct(GEOID, .keep_all = TRUE)

# ---- join ----
eda_dat <- bats %>%
  left_join(covs, by = "GEOID") #keep every county-year-species record, even if covariate missing

# ---- write ----
readr::write_csv(eda_dat, out_csv)
saveRDS(eda_dat, out_rds)

message("Wrote EDA-ready dataset: ", out_csv)
message("Wrote EDA-ready dataset (RDS): ", out_rds)
