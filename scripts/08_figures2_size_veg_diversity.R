# --- scripts/08_figureS2_size_veg_diversity.R --------------------------------
# Builds Figure S2 (supplement): mean +/- SE species richness and abundance by
# wetland size class (small/large) x FS-cover class (Low/Medium/High), summer
# and winter. 
suppressPackageStartupMessages({
  library(here); library(dplyr); library(ggplot2); library(patchwork); library(readr)
})

source(here::here("src/load_data.R"))
source(here::here("src/build_datasets.R"))

dat <- load_all_data()

s.glm.dat <- make_glm_dat(dat$s.habitat, dat$s.indices)
w.glm.dat <- make_glm_dat(dat$w.habitat, dat$w.indices)

summarize_season <- function(df, season_name) {
  df %>%
    dplyr::filter(!is.na(size_class), !is.na(WVclass)) %>%
    dplyr::group_by(size_class, WVclass) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean_rich = mean(richness, na.rm = TRUE),
      se_rich   = sd(richness, na.rm = TRUE) / sqrt(n),
      mean_abun = mean(abundance, na.rm = TRUE),
      se_abun   = sd(abundance, na.rm = TRUE) / sqrt(n),
      .groups = "drop"
    ) %>%
    dplyr::mutate(season = season_name)
}

s.summary <- summarize_season(s.glm.dat, "Summer")
w.summary <- summarize_season(w.glm.dat, "Winter")

cat("=== Summer ===\n"); print(as.data.frame(s.summary))
cat("=== Winter ===\n"); print(as.data.frame(w.summary))

STAMP <- format(Sys.Date(), "%Y-%m-%d")
OUT_DIR <- here::here("results", "figureS2_size_veg", STAMP)
fs::dir_create(OUT_DIR, recurse = TRUE)
write_csv(dplyr::bind_rows(s.summary, w.summary), file.path(OUT_DIR, "figS2_summary.csv"))

tag_theme <- theme(
  plot.tag = element_text(face = "bold", size = 18),
  plot.tag.position = c(0.02, 0.98)
)
fill_vals <- c("0" = "#009E73", "1" = "#E69F00", "2" = "#CC79A7")

bar_plot <- function(df, yvar, sevar, season_name, ylab_text, show_legend = FALSE) {
  df <- df %>% dplyr::mutate(size_class = factor(size_class, levels = c("small","large")))
  ggplot(df, aes(x = size_class, y = .data[[yvar]], fill = WVclass)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
    geom_errorbar(aes(ymin = .data[[yvar]] - .data[[sevar]], ymax = .data[[yvar]] + .data[[sevar]]),
                  position = position_dodge(width = 0.8), width = 0.15, linewidth = 0.5) +
    scale_fill_manual(values = fill_vals, name = "FS cover") +
    labs(x = "Wetland size class", y = ylab_text, title = season_name) +
    theme_bw(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 17),
      axis.text  = element_text(size = 16),
      legend.position = if (show_legend) "bottom" else "none",
      legend.title = element_text(size = 15),
      legend.text  = element_text(size = 14)
    )
}

p_a <- bar_plot(s.summary, "mean_rich", "se_rich", "Summer", "Species richness") + labs(tag = "a") + tag_theme
p_b <- bar_plot(w.summary, "mean_rich", "se_rich", "Winter", NULL) + labs(tag = "b") + tag_theme
p_c <- bar_plot(s.summary, "mean_abun", "se_abun", "Summer", "Abundance") + labs(tag = "c") + tag_theme
p_d <- bar_plot(w.summary, "mean_abun", "se_abun", "Winter", NULL, show_legend = TRUE) + labs(tag = "d") + tag_theme

final_plot <- (p_a + p_b + p_c + p_d) +
  plot_layout(ncol = 4, guides = "collect") &
  theme(legend.position = "bottom")

fs::dir_create(here::here("Figures"))
ggsave(here::here("Figures", "FigureS2_size_veg_diversity.png"), final_plot, width = 14, height = 5.5, dpi = 600)
cat("DONE\n")