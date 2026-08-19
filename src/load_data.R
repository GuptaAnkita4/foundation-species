# src/R/load_data.R
suppressPackageStartupMessages({
  library(readr)
  library(here)
})

read_csv_checked <- function(path, ...) {
  if (!file.exists(path)) {
    stop("load_all_data(): expected input file not found:\n  ", path,
         "\nCheck that data/processed/ has been built (see scripts/00_create_datasets.R) ",
         "and that path_data points at the right directory.", call. = FALSE)
  }
  read_csv(path, ...)
}

load_all_data <- function(path_data = here("data", "processed")) {
  list(
    s.birds.oc = read_csv_checked(file.path(path_data, "summer_birds_occupancy.csv"), show_col_types = FALSE),
    s.birds.ab = read_csv_checked(file.path(path_data, "summer_birds_abundance.csv"), show_col_types = FALSE),
    w.birds.oc = read_csv_checked(file.path(path_data, "winter_birds_occupancy.csv"), show_col_types = FALSE),
    w.birds.ab = read_csv_checked(file.path(path_data, "winter_birds_abundance.csv"), show_col_types = FALSE),
    s.habitat = read_csv_checked(file.path(path_data, "summer_habitat.csv"), show_col_types = FALSE),
    w.habitat = read_csv_checked(file.path(path_data, "winter_habitat.csv"), show_col_types = FALSE),
    s.indices = read_csv_checked(file.path(path_data, "summer_indices.csv"), show_col_types = FALSE),
    w.indices = read_csv_checked(file.path(path_data, "winter_indices.csv"), show_col_types = FALSE),
    waterbirds = read_csv_checked(file.path(path_data, "species_waterbirds.csv"), show_col_types = FALSE),
    func.dat.diet = read_csv_checked(file.path(path_data, "functional_data.csv"), show_col_types = FALSE)
  )
}
