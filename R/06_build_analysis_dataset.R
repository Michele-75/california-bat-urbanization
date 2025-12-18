# R/06_build_analysis_dataset.R
# Purpose:
#   Join bats county×year×species panel to VIIRS county×year radiance,
#   create analysis-ready dataset + QA checks.

source(here::here("R/00_setup.R"))

library(dplyr)
library(readr)
library(tidyr)


# ---- Inputs ----
IN_BATS_RDS  <- file.path(DIR_GBIF_PROC,  "bats_county_year_species.rds")
IN_VIIRS_CSV <- file.path(DIR_VIIRS_PROC, "viirs_county_year.csv")


# ---- Outputs ----
DIR_ANALYSIS_PROC <- here::here("data", "processed", "analysis")
dir.create(DIR_ANALYSIS_PROC, showWarnings = FALSE, recursive = TRUE)

OUT_CSV <- file.path(DIR_ANALYSIS_PROC, "bats_viirs_county_year_species.csv")
OUT_RDS <- file.path(DIR_ANALYSIS_PROC, "bats_viirs_county_year_species.rds")
OUT_QA  <- file.path(DIR_ANALYSIS_PROC, "bats_viirs_join_qa_summary.csv")


#---- Checks ----
  if (!file.exists(IN_BATS_RDS))  stop("Missing bats panel: ", IN_BATS_RDS,  " (run R/05...)")
if (!file.exists(IN_VIIRS_CSV)) stop("Missing VIIRS data: ", IN_VIIRS_CSV, " (run R/03...)")


# ---- Load ----
bats <- readRDS(IN_BATS_RDS)
viirs <- readr::read_csv(IN_VIIRS_CSV, show_col_types = FALSE)


# ---- Basic variable hygiene ----
bats <- bats %>%
  mutate(
    GEOID = as.character(GEOID),
    year  = as.integer(year),
    species_key = as.integer(species_key),
    n_obs        = as.integer(n_obs),
    county_name = as.character(county_name),
    species_label = as.character(species_label)
  )

viirs <- viirs %>%
  mutate(
    GEOID = as.character(GEOID),
    year  = as.integer(year),
    county_name = as.character(county_name),
    vnl_version  = as.character(vnl_version),
    mean_radiance = as.double(mean_radiance)
  ) %>%
  select(GEOID, county_name, year, vnl_version, mean_radiance)


# ---- Join (left join keeps full bat panel incl zeros) ----
dat <- bats %>%
  left_join(viirs, by = c("GEOID", "year"), suffix = c("_bats", "_viirs")) %>% #requires GEOID and year to be the same
  mutate(
    log1p_radiance = log1p(mean_radiance), #very useful transform for later observation,
    #log1p(x)=log(1+x) --> preserves 0 values
    has_viirs = !is.na(mean_radiance) #logical indicator for where VIIRS data is missing
  )

#Only need one county name column
dat <- dat %>%
  mutate(
    county_name = coalesce(county_name_viirs, county_name_bats)
  ) %>%
  select(-county_name_bats, -county_name_viirs) 

#Reorder columns
preferred_order <- c(
  "GEOID", "county_name", "year", "species_key", "species_label",
  "n_obs", "mean_radiance", "log1p_radiance", "has_viirs", "vnl_version"
)

dat <- dat %>%
  dplyr::relocate(dplyr::any_of(preferred_order))


# ---- QA summary ----
qa <- tibble::tibble(
  metric = c(
    "rows_bats_panel",
    "rows_viirs",
    "rows_joined",
    "share_rows_missing_viirs",
    "min_year_bats", "max_year_bats",
    "min_year_viirs", "max_year_viirs"
  ),
  value = c(
    nrow(bats),
    nrow(viirs),
    nrow(dat),
    mean(is.na(dat$mean_radiance)),
    min(bats$year), max(bats$year),
    min(viirs$year), max(viirs$year)
  )
)


# Add quick “missing VIIRS by year” (helps spot a missing GeoTIFF year)
qa_missing_by_year <- dat %>%
  group_by(year) %>%
  summarise(share_missing_viirs = mean(is.na(mean_radiance)), .groups = "drop") %>%
  mutate(metric = paste0("share_missing_viirs_year: ", year)) %>%
  transmute(metric, value = share_missing_viirs)

qa_out <- bind_rows(qa, qa_missing_by_year)
readr::write_csv(qa_out, OUT_QA)


# ---- Save final processed dataset ----
readr::write_csv(dat, OUT_CSV)
saveRDS(dat, OUT_RDS)

message("Saved joined analysis dataset (CSV): ", OUT_CSV)
message("Saved joined analysis dataset (RDS): ", OUT_RDS)
message("Saved QA summary:                   ", OUT_QA)














