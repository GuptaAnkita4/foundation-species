# --- scripts/13_hotspots_area_veg.R ------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(car)
  library(emmeans)
  library(multcompView)
  library(glmmTMB)
  library(here)
  library(fs)
  library(readr)
})

# -----------------------------------------------------------------------------
# Load data & helper to build GLM-ready frames (no FRic)
# -----------------------------------------------------------------------------
source(here::here("src/load_data.R"))

make_glm_df <- function(habitat, indices) {
  stopifnot(all(c("grid","wetlandTot","perWV","WVclass") %in% names(habitat)))
  stopifnot(all(c("grid","richness","count") %in% names(indices)))
  habitat %>%
    dplyr::select(grid, wetlandTot, perWV, WVclass) %>%
    left_join(indices %>% dplyr::select(grid, richness, count), by = "grid") %>%
    mutate(
      logArea   = log2(wetlandTot * 10 + 1),
      abundance = count
    ) %>%
    dplyr::select(grid, logArea, perWV, WVclass, richness, abundance) %>%
    na.omit()
}

dat <- load_all_data()
s.glm.dat <- make_glm_df(dat$s.habitat, dat$s.indices)
w.glm.dat <- make_glm_df(dat$w.habitat, dat$w.indices)

# Output dirs
fs::dir_create(here::here("Figures"))
STAMP    <- format(Sys.Date(), "%Y-%m-%d")
OUT_DIR  <- here::here("results", "hotspots_area_veg", STAMP)
fs::dir_create(OUT_DIR, recurse = TRUE)

# -----------------------------------------------------------------------------
# MODELS (no fric)
# -----------------------------------------------------------------------------
## Summer
s.mod_hotspot.ric <- glmmTMB(richness  ~ logArea * perWV, family = nbinom2, data = s.glm.dat)
s.mod_hotspot.ab  <- glmmTMB(abundance ~ logArea * perWV, family = nbinom2, data = s.glm.dat)

aov_s_ric <- Anova(s.mod_hotspot.ric, type = 3)
aov_s_ab  <- Anova(s.mod_hotspot.ab,  type = 3)

area_seq_s <- seq(min(s.glm.dat$logArea, na.rm = TRUE),
                  max(s.glm.dat$logArea, na.rm = TRUE), length.out = 100)
veg_seq_s  <- seq(min(s.glm.dat$perWV,   na.rm = TRUE),
                  max(s.glm.dat$perWV,   na.rm = TRUE), length.out = 100)
s.pred_grid <- expand.grid(logArea = area_seq_s, perWV = veg_seq_s)
s.pred_grid$pred_richness  <- predict(s.mod_hotspot.ric, newdata = s.pred_grid)
s.pred_grid$pred_abundance <- predict(s.mod_hotspot.ab,  newdata = s.pred_grid)

## Winter
w.mod_hotspot.ric <- glmmTMB(richness  ~ logArea * perWV, family = nbinom2, data = w.glm.dat)
w.mod_hotspot.ab  <- glmmTMB(abundance ~ logArea * perWV, family = nbinom2, data = w.glm.dat)

aov_w_ric <- Anova(w.mod_hotspot.ric, type = 3)
aov_w_ab  <- Anova(w.mod_hotspot.ab,  type = 3)

area_seq_w <- seq(min(w.glm.dat$logArea, na.rm = TRUE),
                  max(w.glm.dat$logArea, na.rm = TRUE), length.out = 100)
veg_seq_w  <- seq(min(w.glm.dat$perWV,   na.rm = TRUE),
                  max(w.glm.dat$perWV,   na.rm = TRUE), length.out = 100)
w.pred_grid <- expand.grid(logArea = area_seq_w, perWV = veg_seq_w)
w.pred_grid$pred_richness  <- predict(w.mod_hotspot.ric, newdata = w.pred_grid)
w.pred_grid$pred_abundance <- predict(w.mod_hotspot.ab,  newdata = w.pred_grid)

# -----------------------------------------------------------------------------
# Binning, reshaping, summaries
# -----------------------------------------------------------------------------
breaks  <- c(-Inf, 2, 4, Inf)
labels  <- c("Small", "Medium", "Large")

s.pred_grid_binned <- s.pred_grid %>%
  mutate(area_bin = cut(logArea, breaks = breaks, labels = labels),
         season   = "Summer")

w.pred_grid_binned <- w.pred_grid %>%
  mutate(area_bin = cut(logArea, breaks = breaks, labels = labels),
         season   = "Winter")

combined_pred <- bind_rows(s.pred_grid_binned, w.pred_grid_binned) %>%
  pivot_longer(cols = c(pred_richness, pred_abundance),
               names_to = "metric", values_to = "value") %>%
  mutate(
    metric   = factor(metric, levels = c("pred_richness", "pred_abundance"),
                      labels = c("Richness", "Abundance")),
    season   = factor(season, levels = c("Summer", "Winter")),
    area_bin = factor(area_bin, levels = c("Small", "Medium", "Large"))
  )

line_summary_sd <- combined_pred %>%
  group_by(area_bin, perWV, season, metric) %>%
  summarize(
    mean_value = mean(value, na.rm = TRUE),
    sd         = sd(value,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(lower = mean_value - sd,
         upper = mean_value + sd)

# -----------------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------------
color_vals <- c("Small" = "#56B4E9", "Medium" = "#F0E442", "Large" = "#009E73")

range_ric <- range(combined_pred %>% filter(metric == "Richness")  %>% pull(value), na.rm = TRUE)
range_ab  <- range(combined_pred %>% filter(metric == "Abundance") %>% pull(value), na.rm = TRUE)

plot_metric_season <- function(df, metric_name, season_name) {
  df %>%
    filter(metric == metric_name, season == season_name) %>%
    ggplot(aes(x = perWV, y = value, group = logArea, color = area_bin)) +
    geom_line(alpha = 0.1, linewidth = 0.6) +
    scale_color_manual(values = color_vals) +
    labs(x = "% Vegetation Cover", y = NULL, color = NULL,
         title = paste(season_name)) +
    theme_classic(base_size = 12) +
    theme(
      axis.line = element_line(color = "black"),
      strip.background = element_blank(),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5)
    )
}

# “many-lines” view
p1 <- plot_metric_season(combined_pred, "Richness",  "Summer") + ylim(range_ric)
p2 <- plot_metric_season(combined_pred, "Abundance", "Summer") + ylim(range_ab)
p3 <- plot_metric_season(combined_pred, "Richness",  "Winter") + ylim(range_ric)
p4 <- plot_metric_season(combined_pred, "Abundance", "Winter") + ylim(range_ab)

final_plot <- (p1 + p2) / (p3 + p4) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Mean±SD ribbons
plot_line_sd_metric <- function(df, metric_name, season_name) {
  df %>%
    filter(metric == metric_name, season == season_name) %>%
    ggplot(aes(x = perWV, y = mean_value, color = area_bin, fill = area_bin)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = color_vals) +
    scale_fill_manual(values = color_vals) +
    labs(
      x = "% Vegetation Cover",
      y = if (metric_name == "Richness") "Log(Richness)" else "Log(Abundance)",
      title = paste(season_name),
      color = "Wetland Size", fill = "Wetland Size"
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11)
    )
}

range_ric_ci <- range(line_summary_sd %>% filter(metric == "Richness")  %>% dplyr::select(lower, upper), na.rm = TRUE)
range_ab_ci  <- range(line_summary_sd %>% filter(metric == "Abundance") %>% dplyr::select(lower, upper), na.rm = TRUE)

q1 <- plot_line_sd_metric(line_summary_sd, "Richness",  "Summer") + ylim(range_ric_ci)
q2 <- plot_line_sd_metric(line_summary_sd, "Abundance", "Summer") + ylim(range_ab_ci)
q3 <- plot_line_sd_metric(line_summary_sd, "Richness",  "Winter")  + ylim(range_ric_ci)
q4 <- plot_line_sd_metric(line_summary_sd, "Abundance", "Winter") + ylim(range_ab_ci)

final_line_ci_plot <- ((q1 + q2) / (q3 + q4)) + plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9)
  )

# -----------------------------------------------------------------------------
# Save outputs
# -----------------------------------------------------------------------------
ggsave(here::here("Figures", "biointeraction_lineplot_manylines.png"),
       final_plot, width = 7, height = 6, dpi = 300)

ggsave(here::here("Figures", "biointeraction_lineplot_with_ci.png"),
       final_line_ci_plot, width = 7, height = 6, dpi = 300)

# Save predictions & ANOVAs
readr::write_csv(s.pred_grid, file.path(OUT_DIR, "summer_pred_grid.csv"))
readr::write_csv(w.pred_grid, file.path(OUT_DIR, "winter_pred_grid.csv"))
readr::write_csv(line_summary_sd, file.path(OUT_DIR, "line_summary_sd.csv"))

dump_anova <- function(aov_obj, path) {
  tb <- as.data.frame(aov_obj)
  tb$term <- rownames(tb); rownames(tb) <- NULL
  readr::write_csv(tb |> dplyr::relocate(term), path)
}
dump_anova(aov_s_ric, file.path(OUT_DIR, "anova_summer_richness.csv"))
dump_anova(aov_s_ab,  file.path(OUT_DIR, "anova_summer_abundance.csv"))
dump_anova(aov_w_ric, file.path(OUT_DIR, "anova_winter_richness.csv"))
dump_anova(aov_w_ab,  file.path(OUT_DIR, "anova_winter_abundance.csv"))

# Save model objects (no fric models)
saveRDS(
  list(
    s_mods = list(ric = s.mod_hotspot.ric, ab = s.mod_hotspot.ab),
    w_mods = list(ric = w.mod_hotspot.ric, ab = w.mod_hotspot.ab)
  ),
  file.path(OUT_DIR, "hotspot_models.rds")
)

message("✅ Outputs written to ", OUT_DIR)
