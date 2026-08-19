# --- scripts/02b_figureS1_render.R -------------------------------------------
# Renders Figure S1 (supplement) as a standalone a-d tagged 4-panel figure,
# reusing the exact SAR/AAR power-law analysis from scripts/02_SAR_AAR.R (Fig.
# 1, panels d-g there). 

suppressPackageStartupMessages({
  library(here); library(dplyr); library(ggplot2); library(patchwork); library(readr)
})

find_latest_sar_aar_dir <- function() {
  root <- here::here("results", "sar_aar")
  dirs <- list.dirs(root, recursive = FALSE)
  if (!length(dirs)) {
    stop("find_latest_sar_aar_dir(): no results/sar_aar/<date>/ folder found -- ",
         "run scripts/02_SAR_AAR.R first (or run via run_all.R, which always runs ",
         "02 before 02b) to produce sar_aar_coefficients.csv and the prediction CSVs.",
         call. = FALSE)
  }

  
  sort(dirs, decreasing = TRUE)[1]
}


source(here::here("src/load_data.R"))
source(here::here("src/build_datasets.R"))

dat <- load_all_data()
s.sparea <- make_sparea(dat$s.habitat, dat$s.indices)
w.sparea <- make_sparea(dat$w.habitat, dat$w.indices)

safe_log <- function(x) ifelse(x > 0, log(x), NA_real_)

s.sparea <- as.data.frame(s.sparea) %>%
  dplyr::mutate(logArea = log2(aT), logRich = safe_log(s), logAbund = safe_log(a))
w.sparea <- as.data.frame(w.sparea) %>%
  dplyr::mutate(logArea = log2(aT), logRich = safe_log(s), logAbund = safe_log(a))

sar_aar_dir <- find_latest_sar_aar_dir()
coefs <- read_csv(file.path(sar_aar_dir, "sar_aar_coefficients.csv"), show_col_types = FALSE)
getc <- function(model, col) coefs[[col]][coefs$model == model]

pred_s_sar <- read_csv(file.path(sar_aar_dir, "pred_summer_richness.csv"), show_col_types = FALSE)
pred_w_sar <- read_csv(file.path(sar_aar_dir, "pred_winter_richness.csv"), show_col_types = FALSE)
pred_s_aar <- read_csv(file.path(sar_aar_dir, "pred_summer_abundance.csv"), show_col_types = FALSE)
pred_w_aar <- read_csv(file.path(sar_aar_dir, "pred_winter_abundance.csv"), show_col_types = FALSE)

p_s_sar <- getc("summer_richness","p");  r2_s_sar <- getc("summer_richness","R2")
p_w_sar <- getc("winter_richness","p");  r2_w_sar <- getc("winter_richness","R2")
p_s_aar <- getc("summer_abundance","p"); r2_s_aar <- getc("summer_abundance","R2")
p_w_aar <- getc("winter_abundance","p"); r2_w_aar <- getc("winter_abundance","R2")

# ---- plot ranges (identical to scripts/02_SAR_AAR.R) ----
xlim_rich <- range(c(s.sparea$logArea, w.sparea$logArea), na.rm = TRUE)
ylim_rich <- range(c(s.sparea$logRich, w.sparea$logRich,
                     pred_s_sar$logPred, pred_w_sar$logPred), na.rm = TRUE)
xlim_abund <- range(c(s.sparea$logArea, w.sparea$logArea), na.rm = TRUE)
ylim_abund <- range(c(s.sparea$logAbund, w.sparea$logAbund,
                      pred_s_aar$logPred, pred_w_aar$logPred), na.rm = TRUE)

okabe_ito <- c(blue = "#0072B2", vermillion = "#D55E00")
my_theme <- theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 17),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5))

fmt_p3 <- function(p, thresh = 0.001) {
  if (is.na(p)) return("NA")
  if (p < thresh) return(paste0("<", format(thresh, nsmall = 3)))
  paste0("=", sprintf("%.3f", p))
}
lab_stats <- function(p, r2) {
  paste0("p", fmt_p3(p), "\nR² = ", sprintf("%.3f", r2))
}
ann_xy <- function(xlim, ylim) {
  list(x = min(xlim) + 0.02 * diff(xlim), y = max(ylim) - 0.02 * diff(ylim))
}
tag_theme <- theme(
  plot.tag = element_text(face = "bold", size = 18),
  plot.tag.position = c(0.02, 0.98)
)

a1 <- ann_xy(xlim_rich, ylim_rich)
p1 <- ggplot(s.sparea, aes(x = logArea, y = logRich)) +
  geom_point(shape = 21, fill = okabe_ito["vermillion"], color = "black", size = 2, stroke = 0.4, na.rm = TRUE) +
  geom_line(data = pred_s_sar, aes(x = logArea, y = logPred), color = okabe_ito["vermillion"], linewidth = 1) +
  labs(title = "Summer", x = expression(log[2]("Wetland Area")), y = expression(log("Richness"))) +
  scale_x_continuous(limits = xlim_rich) + scale_y_continuous(limits = ylim_rich) +
  annotate("text", x = a1$x, y = a1$y, label = lab_stats(p_s_sar, r2_s_sar), hjust = 0, vjust = 1, size = 5) +
  my_theme

a2 <- ann_xy(xlim_rich, ylim_rich)
p2 <- ggplot(w.sparea, aes(x = logArea, y = logRich)) +
  geom_point(shape = 21, fill = okabe_ito["blue"], color = "black", size = 2, stroke = 0.4, na.rm = TRUE) +
  geom_line(data = pred_w_sar, aes(x = logArea, y = logPred), color = okabe_ito["blue"], linewidth = 1) +
  labs(title = "Winter", x = expression(log[2]("Wetland Area")), y = NULL) +
  scale_x_continuous(limits = xlim_rich) + scale_y_continuous(limits = ylim_rich) +
  annotate("text", x = a2$x, y = a2$y, label = lab_stats(p_w_sar, r2_w_sar), hjust = 0, vjust = 1, size = 5) +
  my_theme

b1 <- ann_xy(xlim_abund, ylim_abund)
p3 <- ggplot(s.sparea, aes(x = logArea, y = logAbund)) +
  geom_point(shape = 21, fill = okabe_ito["vermillion"], color = "black", size = 2, stroke = 0.4, na.rm = TRUE) +
  geom_line(data = pred_s_aar, aes(x = logArea, y = logPred), color = okabe_ito["vermillion"], linewidth = 1) +
  labs(title = "Summer", x = expression(log[2]("Wetland Area")), y = expression(log("Abundance"))) +
  scale_x_continuous(limits = xlim_abund) + scale_y_continuous(limits = ylim_abund) +
  annotate("text", x = b1$x, y = b1$y, label = lab_stats(p_s_aar, r2_s_aar), hjust = 0, vjust = 1, size = 5) +
  my_theme

b2 <- ann_xy(xlim_abund, ylim_abund)
p4 <- ggplot(w.sparea, aes(x = logArea, y = logAbund)) +
  geom_point(shape = 21, fill = okabe_ito["blue"], color = "black", size = 2, stroke = 0.4, na.rm = TRUE) +
  geom_line(data = pred_w_aar, aes(x = logArea, y = logPred), color = okabe_ito["blue"], linewidth = 1) +
  labs(title = "Winter", x = expression(log[2]("Wetland Area")), y = NULL) +
  scale_x_continuous(limits = xlim_abund) + scale_y_continuous(limits = ylim_abund) +
  annotate("text", x = b2$x, y = b2$y, label = lab_stats(p_w_aar, r2_w_aar), hjust = 0, vjust = 1, size = 5) +
  my_theme

p1 <- p1 + labs(tag = "a") + tag_theme
p2 <- p2 + labs(tag = "b") + tag_theme
p3 <- p3 + labs(tag = "c") + tag_theme
p4 <- p4 + labs(tag = "d") + tag_theme

final_plot <- (p1 + p2 + plot_spacer() + p3 + p4) +
  plot_layout(ncol = 5, widths = c(1, 1, 0.05, 1, 1))

fs::dir_create(here::here("Figures"))
ggsave(here::here("Figures", "FigureS1_SAR_AAR.png"), final_plot, width = 14, height = 5, dpi = 600)
cat("Saved. p/R2 check:\n")
cat(sprintf("summer_richness  p=%s R2=%.3f\n", fmt_p3(p_s_sar), r2_s_sar))
cat(sprintf("winter_richness  p=%s R2=%.3f\n", fmt_p3(p_w_sar), r2_w_sar))
cat(sprintf("summer_abundance p=%s R2=%.3f\n", fmt_p3(p_s_aar), r2_s_aar))
cat(sprintf("winter_abundance p=%s R2=%.3f\n", fmt_p3(p_w_aar), r2_w_aar))
cat("DONE\n")
