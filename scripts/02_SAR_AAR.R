# --- scripts/02_sar_aar.R ----------------------------------------------------

pacman::p_load(
  ggplot2, dplyr, readr, sars, cowplot, glmmTMB, patchwork,
  here, fs, tibble
)

######################################################
##-------Testing Species-Area Relationships---------##
######################################################

# ---- Fit SAR/AAR with sars::sar_power on RAW scale ----
w.sar <- sar_power(as.matrix(w.sparea[, c("aT", "s")]))  # SAR
w.aar <- sar_power(as.matrix(w.sparea[, c("aT", "a")]))  # AAR
s.sar <- sar_power(as.matrix(s.sparea[, c("aT", "s")]))
s.aar <- sar_power(as.matrix(s.sparea[, c("aT", "a")]))
# (Summaries available via summary(...))

# ---- Data prep for plotting (log₂ A on x; ln S on y) ----
s.sparea <- as.data.frame(s.sparea)
w.sparea <- as.data.frame(w.sparea)

safe_log <- function(x) { # plotting-only: drop zeros
  ifelse(x > 0, log(x), NA_real_)
}

s.sparea <- s.sparea %>%
  mutate(logArea  = log2(aT),
         logRich  = safe_log(s),
         logAbund = safe_log(a))

w.sparea <- w.sparea %>%
  mutate(logArea  = log2(aT),
         logRich  = safe_log(s),
         logAbund = safe_log(a))

# ---- Extract coefficients (c, z) and R2 / p (from sars object) ----
c_s_sar <- s.sar$par[1]; z_s_sar <- s.sar$par[2]; p_s_sar <- s.sar$sigConf[8]; r2_s_sar <- s.sar$R2
c_w_sar <- w.sar$par[1]; z_w_sar <- w.sar$par[2]; p_w_sar <- w.sar$sigConf[8]; r2_w_sar <- w.sar$R2
c_s_aar <- s.aar$par[1]; z_s_aar <- s.aar$par[2]; p_s_aar <- s.aar$sigConf[8]; r2_s_aar <- s.aar$R2
c_w_aar <- w.aar$par[1]; z_w_aar <- w.aar$par[2]; p_w_aar <- w.aar$sigConf[8]; r2_w_aar <- w.aar$R2

# ---- Predictions on log₂(A) axis with CORRECT slope scaling ----
# 🔧 ln(S) = ln(c) + z * ln(2) * log2(A)
ln2 <- log(2)
make_pred <- function(c_, z_, log2A_seq) {
  tibble(logArea = log2A_seq,
         logPred = log(c_) + z_ * ln2 * log2A_seq)
}

area_seq_s <- seq(min(s.sparea$logArea, na.rm = TRUE),
                  max(s.sparea$logArea, na.rm = TRUE), length.out = 200)
area_seq_w <- seq(min(w.sparea$logArea, na.rm = TRUE),
                  max(w.sparea$logArea, na.rm = TRUE), length.out = 200)

pred_s_sar <- make_pred(c_s_sar, z_s_sar, area_seq_s)
pred_w_sar <- make_pred(c_w_sar, z_w_sar, area_seq_w)
pred_s_aar <- make_pred(c_s_aar, z_s_aar, area_seq_s)
pred_w_aar <- make_pred(c_w_aar, z_w_aar, area_seq_w)

# ---- Optional: quick 95% CI ribbons via normal approx on (c, z) ----
# If you have SEs for c and z, set them here (from your NLS summaries);
# if not available, you can skip the ribbon or replace with bootstrap draws
# from your preferred fitting that returns vcov.
ci_ribbon <- function(c_, z_, log2A_seq, se_c = NA_real_, se_z = NA_real_, B = 1000L) {
  if (is.na(se_c) || is.na(se_z)) return(NULL)
  draws <- tibble(
    c_draw = rnorm(B, mean = c_, sd = se_c),
    z_draw = rnorm(B, mean = z_, sd = se_z)
  )
  grid <- tibble(logArea = log2A_seq)
  M <- sapply(1:B, function(b) log(draws$c_draw[b]) + draws$z_draw[b] * ln2 * log2A_seq)
  tibble(logArea = log2A_seq,
         lo = apply(M, 1, quantile, 0.025, na.rm = TRUE),
         hi = apply(M, 1, quantile, 0.975, na.rm = TRUE))
}
# Example: if you have SEs, plug them in; otherwise set ribbons <- NULL
ribbons <- list(
  s_sar = NULL, w_sar = NULL, s_aar = NULL, w_aar = NULL
)

# ---- Plot ranges ----
xlim_rich <- range(c(s.sparea$logArea, w.sparea$logArea), na.rm = TRUE)
ylim_rich <- range(c(s.sparea$logRich, w.sparea$logRich,
                     pred_s_sar$logPred, pred_w_sar$logPred), na.rm = TRUE)
xlim_abund <- range(c(s.sparea$logArea, w.sparea$logArea), na.rm = TRUE)
ylim_abund <- range(c(s.sparea$logAbund, w.sparea$logAbund,
                      pred_s_aar$logPred, pred_w_aar$logPred), na.rm = TRUE)

# ---- Theme + utils ----
okabe_ito <- c(blue = "#0072B2", vermillion = "#D55E00")
my_theme <- theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

fmt_p3 <- function(p, thresh = 0.001) {
  if (is.na(p)) return("NA")
  if (p < thresh) return(paste0("<", format(thresh, nsmall = 3)))
  paste0("=", sprintf("%.3f", p))
}
# lab_stats <- function(z, p, r2) {
#   paste0("z = ", round(z, 3),
#          ", p = ", fmt_p3(p),
#          "\nR\u00b2 = ", round(r2, 3))
# }
lab_stats <- function(p) {
  paste0("p", fmt_p3(p))
}
ann_xy <- function(xlim, ylim) {
  list(x = min(xlim) + 0.02 * diff(xlim),
       y = max(ylim) - 0.02 * diff(ylim))
}
# --- subplot tag styling (bold letter in top-left) ---
tag_theme <- theme(
  plot.tag = element_text(face = "bold", size = 12),
  plot.tag.position = c(0.02, 0.98)  # x,y in [0,1] npc
)


# ---- Plots (note axis labels changed) ----
a1 <- ann_xy(xlim_rich, ylim_rich)
p1 <- ggplot(s.sparea, aes(x = logArea, y = logRich)) +
  geom_point(shape = 21, fill = okabe_ito["vermillion"], color = "black", size = 2, stroke = 0.4, na.rm = TRUE) +
  { if (!is.null(ribbons$s_sar)) geom_ribbon(data = ribbons$s_sar, aes(x = logArea, ymin = lo, ymax = hi), inherit.aes = FALSE, alpha = 0.15) } +
  geom_line(data = pred_s_sar, aes(x = logArea, y = logPred), color = okabe_ito["vermillion"], linewidth = 1) +
  labs(title = "Summer", x = expression(log[2]("Wetland Area")), y = expression(log("Richness"))) +
  scale_x_continuous(limits = xlim_rich) +
  scale_y_continuous(limits = ylim_rich) +
  # annotate("text", x = a1$x, y = a1$y, label = lab_stats(z_s_sar, p_s_sar, r2_s_sar),
  #          hjust = 0, vjust = 1, size = 3.5) +
  annotate("text", x = a1$x, y = a1$y, label = lab_stats(p_s_sar),
           hjust = 0, vjust = 1, size = 3.5) +
  my_theme

a2 <- ann_xy(xlim_rich, ylim_rich)
p2 <- ggplot(w.sparea, aes(x = logArea, y = logRich)) +
  geom_point(shape = 21, fill = okabe_ito["blue"], color = "black", size = 2, stroke = 0.4, na.rm = TRUE) +
  { if (!is.null(ribbons$w_sar)) geom_ribbon(data = ribbons$w_sar, aes(x = logArea, ymin = lo, ymax = hi), inherit.aes = FALSE, alpha = 0.15) } +
  geom_line(data = pred_w_sar, aes(x = logArea, y = logPred), color = okabe_ito["blue"], linewidth = 1) +
  labs(title = "Winter", x = expression(log[2]("Wetland Area")), y = NULL) +
  scale_x_continuous(limits = xlim_rich) +
  scale_y_continuous(limits = ylim_rich) +
  # annotate("text", x = a2$x, y = a2$y, label = lab_stats(z_w_sar, p_w_sar, r2_w_sar),
  #          hjust = 0, vjust = 1, size = 3.5) +
  annotate("text", x = a2$x, y = a2$y, label = lab_stats(p_w_sar),
           hjust = 0, vjust = 1, size = 3.5) +
  my_theme

b1 <- ann_xy(xlim_abund, ylim_abund)
p3 <- ggplot(s.sparea, aes(x = logArea, y = logAbund)) +
  geom_point(shape = 21, fill = okabe_ito["vermillion"], color = "black", size = 2, stroke = 0.4, na.rm = TRUE) +
  { if (!is.null(ribbons$s_aar)) geom_ribbon(data = ribbons$s_aar, aes(x = logArea, ymin = lo, ymax = hi), inherit.aes = FALSE, alpha = 0.15) } +
  geom_line(data = pred_s_aar, aes(x = logArea, y = logPred), color = okabe_ito["vermillion"], linewidth = 1) +
  labs(title = "Summer", x = expression(log[2]("Wetland Area")), y = expression(log("Abundance"))) +
  scale_x_continuous(limits = xlim_abund) +
  scale_y_continuous(limits = ylim_abund) +
  # annotate("text", x = b1$x, y = b1$y, label = lab_stats(z_s_aar, p_s_aar, r2_s_aar),
  #          hjust = 0, vjust = 1, size = 3.5) +
  annotate("text", x = b1$x, y = b1$y, label = lab_stats(p_s_aar),
           hjust = 0, vjust = 1, size = 3.5) +
  my_theme

b2 <- ann_xy(xlim_abund, ylim_abund)
p4 <- ggplot(w.sparea, aes(x = logArea, y = logAbund)) +
  geom_point(shape = 21, fill = okabe_ito["blue"], color = "black", size = 2, stroke = 0.4, na.rm = TRUE) +
  { if (!is.null(ribbons$w_aar)) geom_ribbon(data = ribbons$w_aar, aes(x = logArea, ymin = lo, ymax = hi), inherit.aes = FALSE, alpha = 0.15) } +
  geom_line(data = pred_w_aar, aes(x = logArea, y = logPred), color = okabe_ito["blue"], linewidth = 1) +
  labs(title = "Winter", x = expression(log[2]("Wetland Area")), y = NULL) +
  scale_x_continuous(limits = xlim_abund) +
  scale_y_continuous(limits = ylim_abund) +
  # annotate("text", x = b2$x, y = b2$y, label = lab_stats(z_w_aar, p_w_aar, r2_w_aar),
  #          hjust = 0, vjust = 1, size = 3.5) +
  annotate("text", x = b2$x, y = b2$y, label = lab_stats(p_w_aar),
           hjust = 0, vjust = 1, size = 3.5) +
  my_theme

p1 <- p1 + labs(tag = "d") + tag_theme
p2 <- p2 + labs(tag = "e") + tag_theme
p3 <- p3 + labs(tag = "f") + tag_theme
p4 <- p4 + labs(tag = "g") + tag_theme


final_plot <- (p1 + p2 + plot_spacer() + p3 + p4) +
  plot_layout(ncol = 5, widths = c(1, 1, 0.05, 1, 1))

print(final_plot)

# ---- Save figure ----
fs::dir_create(here::here("Figures"))
ggsave(here::here("Figures", "SAR_AAR.png"), final_plot, width = 14, height = 5, dpi = 300)

# ---- Save results ----
OUT_ROOT <- here::here("results", "sar_aar")
STAMP    <- format(Sys.Date(), "%Y-%m-%d")
OUT_DIR  <- file.path(OUT_ROOT, STAMP)
fs::dir_create(OUT_DIR, recurse = TRUE)

coef_tbl <- dplyr::bind_rows(
  tibble(model = "summer_richness",  c = c_s_sar,  z = z_s_sar,  p = p_s_sar,  R2 = r2_s_sar),
  tibble(model = "winter_richness",  c = c_w_sar,  z = z_w_sar,  p = p_w_sar,  R2 = r2_w_sar),
  tibble(model = "summer_abundance", c = c_s_aar,  z = z_s_aar,  p = p_s_aar,  R2 = r2_s_aar),
  tibble(model = "winter_abundance", c = c_w_aar,  z = z_w_aar,  p = p_w_aar,  R2 = r2_w_aar)
)
readr::write_csv(coef_tbl, file.path(OUT_DIR, "sar_aar_coefficients.csv"))

readr::write_csv(pred_s_sar, file.path(OUT_DIR, "pred_summer_richness.csv"))
readr::write_csv(pred_w_sar, file.path(OUT_DIR, "pred_winter_richness.csv"))
readr::write_csv(pred_s_aar, file.path(OUT_DIR, "pred_summer_abundance.csv"))
readr::write_csv(pred_w_aar, file.path(OUT_DIR, "pred_winter_abundance.csv"))

ggsave(file.path(OUT_DIR, "SAR_AAR.png"), final_plot, width = 14, height = 5, dpi = 300)

saveRDS(list(s_sar = s.sar, s_aar = s.aar, w_sar = w.sar, w_aar = w.aar),
        file.path(OUT_DIR, "sar_aar_sar_power_objects.rds"))

rm(s.sar, s.aar, w.sar, w.aar)
message("✅ SAR/AAR results written to: ", OUT_DIR)
