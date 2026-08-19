# src/eda_utils.R
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(stringr)
  library(forcats)
})

theme_min_clean <- function() {
  theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom")
}

# join habitat + indices by grid, rename columns consistently
make_season_df <- function(habitat, indices) {
  habitat %>%
    dplyr::select(grid, wetlandTot, logArea, perWV, WVclass) %>%
    dplyr::left_join(
      indices %>%
        dplyr::select(grid, richness, FRic, count) %>%
        dplyr::rename(fric = FRic, abundance = count),
      by = "grid"
    ) %>%
    dplyr::rename(area = wetlandTot) %>%
    dplyr::filter(!is.na(area), area > 0) %>%
    dplyr::mutate(
      WVclass   = as.factor(WVclass),
      logArea   = as.numeric(logArea),
      richness  = as.numeric(richness),
      fric      = as.numeric(fric),
      abundance = as.numeric(abundance)
    )
}

plot_scatter_loess <- function(df, x, y, title, xlab, ylab) {
  ggplot(df, aes({{ x }}, {{ y }})) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "loess", se = TRUE) +
    labs(title = title, x = xlab, y = ylab) +
    theme_min_clean()
}

# Summarize abundance from a wide community matrix (first col = grid)
sum_abundance_by_traits <- function(comm_mat, species_info, season_label) {
  stopifnot(ncol(comm_mat) >= 2)
  # total individuals per species across all grids
  totals <- comm_mat %>%
    dplyr::mutate(across(-1, ~ suppressWarnings(as.numeric(.x)))) %>%
    dplyr::summarise(across(-1, ~ sum(.x, na.rm = TRUE))) %>%
    tidyr::pivot_longer(everything(), names_to = "Code", values_to = "abundance")
  
  totals %>%
    dplyr::left_join(species_info, by = "Code") %>%
    dplyr::group_by(trophLevel, residentStat, .drop = FALSE) %>%
    dplyr::summarise(total_abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(season = season_label)
}
