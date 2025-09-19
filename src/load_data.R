# src/R/load_data.R
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(here)
})

load_all_data <- function(path_data = here("data", "processed")) {
  list(
    s_habitat = read_csv(file.path(path_data, "summer_habitat.csv"), show_col_types = FALSE),
    w_habitat = read_csv(file.path(path_data, "winter_habitat.csv"), show_col_types = FALSE),
    s_indices = read_csv(file.path(path_data, "summer_indices.csv"), show_col_types = FALSE),
    w_indices = read_csv(file.path(path_data, "winter_indices.csv"), show_col_types = FALSE),
    waterbirds = read_csv(file.path(path_data, "species_waterbirds.csv"), show_col_types = FALSE),
    func_dat_diet = read_csv(file.path(path_data, "func_dat_diet.csv"), show_col_types = FALSE)
  )
}
