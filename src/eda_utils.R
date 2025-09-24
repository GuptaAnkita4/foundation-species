# eda_utils.R
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
    left_join(
      indices %>%
        dplyr::select(grid, richness, FRic, count) %>%
        rename(fric = FRic, abundance = count),
      by = "grid"
    ) %>%
    rename(area = wetlandTot) %>%
    filter(!is.na(area), area > 0) %>%
    mutate(
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
    mutate(across(-1, ~ suppressWarnings(as.numeric(.x)))) %>%
    summarise(across(-1, ~ sum(.x, na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "Code", values_to = "abundance")
  
  totals %>%
    left_join(species_info, by = "Code") %>%
    group_by(trophLevel, residentStat, .drop = FALSE) %>%
    summarise(total_abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
    mutate(season = season_label)
}
