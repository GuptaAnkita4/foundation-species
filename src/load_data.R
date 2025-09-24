# src/R/load_data.R
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(here)
})

load_all_data <- function(path_data = here("data", "processed")) {
  list(
    s.birds.oc = read_csv(file.path(path_data, "summer_birds_occupancy.csv"), show_col_types = FALSE),
    s.birds.ab = read_csv(file.path(path_data, "summer_birds_abundance.csv"), show_col_types = FALSE),
    w.birds.oc = read_csv(file.path(path_data, "winter_birds_occupancy.csv"), show_col_types = FALSE),
    w.birds.ab = read_csv(file.path(path_data, "winter_birds_abundance.csv"), show_col_types = FALSE),
    s.habitat = read_csv(file.path(path_data, "summer_habitat.csv"), show_col_types = FALSE),
    w.habitat = read_csv(file.path(path_data, "winter_habitat.csv"), show_col_types = FALSE),
    s.indices = read_csv(file.path(path_data, "summer_indices.csv"), show_col_types = FALSE),
    w.indices = read_csv(file.path(path_data, "winter_indices.csv"), show_col_types = FALSE),
    waterbirds = read_csv(file.path(path_data, "species_waterbirds.csv"), show_col_types = FALSE),
    func.dat.diet = read_csv(file.path(path_data, "functional_data.csv"), show_col_types = FALSE)
  )
}
