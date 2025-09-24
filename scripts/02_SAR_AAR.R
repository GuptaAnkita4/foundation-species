# --- scripts/02_sar_aar.R ----------------------------------------------------

# Load required libraries
pacman::p_load(
  ggplot2, dplyr, readr, sars, cowplot, glmmTMB, patchwork,
  here, fs
)


######################################################
##-------Testing Species-Area Relationships---------##
######################################################

# ---- Analyzing Winter Species-Area Relationships ----
# Winter
w.sar <- sar_power(as.matrix(w.sparea[, c("aT", "s")]))  # SAR
w.aar <- sar_power(as.matrix(w.sparea[, c("aT", "a")]))  # AAR

summary(w.sar)
summary(w.aar)

# ---- Analyzing Summer Species-Area Relationships ----
# Summer
s.sar <- sar_power(as.matrix(s.sparea[, c("aT", "s")]))
s.aar <- sar_power(as.matrix(s.sparea[, c("aT", "a")]))

summary(s.sar)
summary(s.aar)

# Data prep for plotting
s.sparea <- data.frame(s.sparea)
w.sparea <- data.frame(w.sparea)

s.sparea$logArea  <- log(s.sparea$aT)
s.sparea$logRich  <- log(s.sparea$s)
s.sparea$logAbund <- log(s.sparea$a)

w.sparea$logArea  <- log(w.sparea$aT)
w.sparea$logRich  <- log(w.sparea$s)
w.sparea$logAbund <- log(w.sparea$a)

# Extract coefficients
c_s_sar <- s.sar$par[1]; z_s_sar <- s.sar$par[2]; p_s_sar <- s.sar$sigConf[8]; r2_s_sar <- s.sar$R2
c_w_sar <- w.sar$par[1]; z_w_sar <- w.sar$par[2]; p_w_sar <- w.sar$sigConf[8]; r2_w_sar <- w.sar$R2
c_s_aar <- s.aar$par[1]; z_s_aar <- s.aar$par[2]; p_s_aar <- s.aar$sigConf[8]; r2_s_aar <- s.aar$R2
c_w_aar <- w.aar$par[1]; z_w_aar <- w.aar$par[2]; p_w_aar <- w.aar$sigConf[8]; r2_w_aar <- w.aar$R2

# Predictions
make_pred <- function(c, z, logArea_seq) {
  data.frame(logArea = logArea_seq, logPred = log(c) + z * logArea_seq)
}

area_seq_s <- seq(min(s.sparea$logArea), max(s.sparea$logArea), length.out = 100)
area_seq_w <- seq(min(w.sparea$logArea), max(w.sparea$logArea), length.out = 100)

pred_s_sar <- make_pred(c_s_sar, z_s_sar, area_seq_s)
pred_w_sar <- make_pred(c_w_sar, z_w_sar, area_seq_w)
pred_s_aar <- make_pred(c_s_aar, z_s_aar, area_seq_s)
pred_w_aar <- make_pred(c_w_aar, z_w_aar, area_seq_w)

# Plot ranges
xlim_rich <- range(c(s.sparea$logArea, w.sparea$logArea), na.rm = TRUE)
ylim_rich <- range(c(s.sparea$logRich, w.sparea$logRich), na.rm = TRUE)

xlim_abund <- range(c(s.sparea$logArea, w.sparea$logArea), na.rm = TRUE)
ylim_abund <- range(c(s.sparea$logAbund, w.sparea$logAbund), na.rm = TRUE)

# Theme + colors
okabe_ito <- c(blue = "#0072B2", vermillion = "#D55E00")
my_theme <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  )

# Helper to place annotation consistently
ann_xy <- function(xlim, ylim) {
  list(
    x = min(xlim) + 0.02 * diff(xlim),
    y = max(ylim) - 0.02 * diff(ylim)
  )
}

# Plots
a1 <- ann_xy(xlim_rich, ylim_rich)
p1 <- ggplot(s.sparea, aes(x = logArea, y = logRich)) +
  geom_point(shape = 21, fill = okabe_ito["vermillion"], color = "black", size = 2, stroke = 0.4) +
  geom_line(data = pred_s_sar, aes(x = logArea, y = logPred), color = okabe_ito["vermillion"], size = 1) +
  labs(title = "Summer", x = "log(Wetland Area)", y = "log(Species Richness)") +
  scale_x_continuous(limits = xlim_rich) +
  scale_y_continuous(limits = ylim_rich) +
  annotate("text", x = a1$x, y = a1$y,
           label = paste0("z=", round(z_s_sar, 2),
                          ", p", ifelse(p_s_sar < 0.001, "<0.001",
                                        ifelse(p_s_sar < 0.01, "<0.01",
                                               ifelse(p_s_sar < 0.05, "<0.05", "\u2265 0.05"))),
                          "\nR\u00b2 = ", round(r2_s_sar, 2)),
           hjust = 0, vjust = 1, size = 3.5) +
  my_theme

a2 <- ann_xy(xlim_rich, ylim_rich)
p2 <- ggplot(w.sparea, aes(x = logArea, y = logRich)) +
  geom_point(shape = 21, fill = okabe_ito["blue"], color = "black", size = 2, stroke = 0.4) +
  geom_line(data = pred_w_sar, aes(x = logArea, y = logPred), color = okabe_ito["blue"], size = 1) +
  labs(title = "Winter", x = "log(Wetland Area)", y = NULL) +
  scale_x_continuous(limits = xlim_rich) +
  scale_y_continuous(limits = ylim_rich) +
  annotate("text", x = a2$x, y = a2$y,
           label = paste0("z=", round(z_w_sar, 2),
                          ", p", ifelse(p_w_sar < 0.001, "<0.001",
                                        ifelse(p_w_sar < 0.01, "<0.01",
                                               ifelse(p_w_sar < 0.05, "<0.05", "\u2265 0.05"))),
                          "\nR\u00b2 = ", round(r2_w_sar, 2)),
           hjust = 0, vjust = 1, size = 3.5) +
  my_theme

b1 <- ann_xy(xlim_abund, ylim_abund)
p3 <- ggplot(s.sparea, aes(x = logArea, y = logAbund)) +
  geom_point(shape = 21, fill = okabe_ito["vermillion"], color = "black", size = 2, stroke = 0.4) +
  geom_line(data = pred_s_aar, aes(x = logArea, y = logPred), color = okabe_ito["vermillion"], size = 1) +
  labs(title = "Summer", x = "log(Wetland Area)", y = "log(Abundance)") +
  scale_x_continuous(limits = xlim_abund) +
  scale_y_continuous(limits = ylim_abund) +
  annotate("text", x = b1$x, y = b1$y,
           label = paste0("z=", round(z_s_aar, 2),
                          ", p", ifelse(p_s_aar < 0.001, "<0.001",
                                        ifelse(p_s_aar < 0.01, "<0.01",
                                               ifelse(p_s_aar < 0.05, "<0.05", "\u2265 0.05"))),
                          "\nR\u00b2 = ", round(r2_s_aar, 2)),
           hjust = 0, vjust = 1, size = 3.5) +
  my_theme

b2 <- ann_xy(xlim_abund, ylim_abund)
p4 <- ggplot(w.sparea, aes(x = logArea, y = logAbund)) +
  geom_point(shape = 21, fill = okabe_ito["blue"], color = "black", size = 2, stroke = 0.4) +
  geom_line(data = pred_w_aar, aes(x = logArea, y = logPred), color = okabe_ito["blue"], size = 1) +
  labs(title = "Winter", x = "log(Wetland Area)", y = NULL) +
  scale_x_continuous(limits = xlim_abund) +
  scale_y_continuous(limits = ylim_abund) +
  annotate("text", x = b2$x, y = b2$y,
           label = paste0("z=", round(z_w_aar, 2),
                          ", p", ifelse(p_w_aar < 0.001, "<0.001",
                                        ifelse(p_w_aar < 0.01, "<0.01",
                                               ifelse(p_w_aar < 0.05, "<0.05", "\u2265 0.05"))),
                          "\nR\u00b2 = ", round(r2_w_aar, 2)),
           hjust = 0, vjust = 1, size = 3.5) +
  my_theme

final_plot <- (p1 + p2 + plot_spacer() + p3 + p4) +
  plot_layout(ncol = 5, widths = c(1, 1, 0.05, 1, 1))

print(final_plot)

# Save figure
fs::dir_create(here::here("Figures"))
ggsave(here::here("Figures", "SAR_AAR.png"), final_plot,
       width = 14, height = 5, dpi = 300)

# ---------------- Save results ----------------
OUT_ROOT <- here::here("results", "sar_aar")
STAMP    <- format(Sys.Date(), "%Y-%m-%d")
OUT_DIR  <- file.path(OUT_ROOT, STAMP)
fs::dir_create(OUT_DIR, recurse = TRUE)

coef_tbl <- dplyr::bind_rows(
  tibble::tibble(model = "summer_richness",  c = c_s_sar,  z = z_s_sar,  p = p_s_sar,  R2 = r2_s_sar),
  tibble::tibble(model = "winter_richness",  c = c_w_sar,  z = z_w_sar,  p = p_w_sar,  R2 = r2_w_sar),
  tibble::tibble(model = "summer_abundance", c = c_s_aar,  z = z_s_aar,  p = p_s_aar,  R2 = r2_s_aar),
  tibble::tibble(model = "winter_abundance", c = c_w_aar,  z = z_w_aar,  p = p_w_aar,  R2 = r2_w_aar)
)
readr::write_csv(coef_tbl, file.path(OUT_DIR, "sar_aar_coefficients.csv"))

readr::write_csv(pred_s_sar, file.path(OUT_DIR, "pred_summer_richness.csv"))
readr::write_csv(pred_w_sar, file.path(OUT_DIR, "pred_winter_richness.csv"))
readr::write_csv(pred_s_aar, file.path(OUT_DIR, "pred_summer_abundance.csv"))
readr::write_csv(pred_w_aar, file.path(OUT_DIR, "pred_winter_abundance.csv"))

ggsave(file.path(OUT_DIR, "SAR_AAR.png"), final_plot,
       width = 14, height = 5, dpi = 300)

saveRDS(
  list(s_sar = s.sar, s_aar = s.aar, w_sar = w.sar, w_aar = w.aar),
  file.path(OUT_DIR, "sar_aar_sar_power_objects.rds")
)

rm(s.sar, s.aar, w.sar, w.aar)

message("✅ SAR/AAR results written to: ", OUT_DIR)
