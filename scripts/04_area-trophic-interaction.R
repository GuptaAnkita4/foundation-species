# scripts/04_trophic_SAR_by_level.R
suppressPackageStartupMessages({
  library(here)
  library(dplyr); library(tidyr); library(ggplot2); library(cowplot)
  library(car); library(emmeans); library(multcomp); library(patchwork)
  library(MASS); library(glmmTMB); library(ggeffects)
  library(readr); library(fs); library(purrr)
})

# ---------------------------------------------------------------------
# Preconditions (stop early if required columns missing)
# ---------------------------------------------------------------------
need_rich_cols <- c("grid","logArea","logWV","trophLevel","richness")
need_cnt_cols  <- c("grid","logArea","logWV","trophLevel","count")

stopifnot(
  exists("s.sparea.troph"), exists("w.sparea.troph"),
  exists("s.spcount.troph"), exists("w.spcount.troph")
)

stopifnot(all(need_rich_cols %in% names(s.sparea.troph)))
stopifnot(all(need_rich_cols %in% names(w.sparea.troph)))
stopifnot(all(need_cnt_cols  %in% names(s.spcount.troph)))
stopifnot(all(need_cnt_cols  %in% names(w.spcount.troph)))

dir.create(here("Figures"), showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------
# Labels for legends
# ---------------------------------------------------------------------
trophrich_labels <- c(
  "richness_2" = "Primary consumers",
  "richness_3" = "Secondary consumers",
  "richness_4" = "Top carnivores"
)
trophcount_labels <- c(
  "count_2" = "Primary consumers",
  "count_3" = "Secondary consumers",
  "count_4" = "Top carnivores"
)

# ---------------------------------------------------------------------
# Models (NB2; unchanged)
# ---------------------------------------------------------------------
troph_rich_A.s <- glmmTMB(
  richness ~ logArea * factor(trophLevel) + (1|grid),
  data = s.sparea.troph, family = nbinom2,
  control = glmmTMBControl(optimizer = optim, optArgs = list(method="BFGS"))
)
troph_rich_WV.s <- glmmTMB(
  richness ~ logWV * factor(trophLevel) + (1|grid),
  data = s.sparea.troph, family = nbinom2,
  control = glmmTMBControl(optimizer = optim, optArgs = list(method="BFGS"))
)
troph_count_A.s <- glmmTMB(
  count ~ logArea * factor(trophLevel) + (1|grid),
  data = s.spcount.troph, family = nbinom2,
  control = glmmTMBControl(optimizer = optim, optArgs = list(method="BFGS"))
)
troph_count_WV.s <- glmmTMB(
  count ~ logWV * factor(trophLevel) + (1|grid),
  data = s.spcount.troph, family = nbinom2,
  control = glmmTMBControl(optimizer = optim, optArgs = list(method="BFGS"))
)

troph_rich_A.w <- glmmTMB(
  richness ~ logArea * factor(trophLevel) + (1|grid),
  data = w.sparea.troph, family = nbinom2,
  control = glmmTMBControl(optimizer = optim, optArgs = list(method="BFGS"))
)
troph_rich_WV.w <- glmmTMB(
  richness ~ logWV * factor(trophLevel) + (1|grid),
  data = w.sparea.troph, family = nbinom2,
  control = glmmTMBControl(optimizer = optim, optArgs = list(method="BFGS"))
)
troph_count_A.w <- glmmTMB(
  count ~ logArea * factor(trophLevel) + (1|grid),
  data = w.spcount.troph, family = nbinom2,
  control = glmmTMBControl(optimizer = optim, optArgs = list(method="BFGS"))
)
troph_count_WV.w <- glmmTMB(
  count ~ logWV * factor(trophLevel) + (1|grid),
  data = w.spcount.troph, family = nbinom2,
  control = glmmTMBControl(optimizer = optim, optArgs = list(method="BFGS"))
)

# ---------------------------------------------------------------------
# Predictions (AREA terms for the 2x2 figure)
# ---------------------------------------------------------------------
pred_rich_A.s  <- ggpredict(troph_rich_A.s,  terms = c("logArea", "trophLevel"), bias_correction = TRUE)
pred_rich_A.w  <- ggpredict(troph_rich_A.w,  terms = c("logArea", "trophLevel"), bias_correction = TRUE)
pred_count_A.s <- ggpredict(troph_count_A.s, terms = c("logArea", "trophLevel"), bias_correction = TRUE)
pred_count_A.w <- ggpredict(troph_count_A.w, terms = c("logArea", "trophLevel"), bias_correction = TRUE)

# Keep prediction groups as expected
pred_rich_A.s$group  <- factor(pred_rich_A.s$group,  levels = c("richness_2","richness_3","richness_4"))
pred_rich_A.w$group  <- factor(pred_rich_A.w$group,  levels = c("richness_2","richness_3","richness_4"))
pred_count_A.s$group <- factor(pred_count_A.s$group, levels = c("count_2","count_3","count_4"))
pred_count_A.w$group <- factor(pred_count_A.w$group, levels = c("count_2","count_3","count_4"))

# ---------------------------------------------------------------------
# Plot transforms: pure log for everything (models unchanged)
# ---------------------------------------------------------------------
EPS_RICH <- 0.5   # half-count for richness points
EPS_CNT  <- 0.5   # half-count for abundance points

# points transform (avoid -Inf at zeros)
tlog_point_rich <- function(y) log(pmax(y, EPS_RICH))
tlog_point_cnt  <- function(y) log(pmax(y, EPS_CNT))

# lines transform: log of predicted mean (no epsilon)
tlog_line <- function(mu) log(mu)

y_lab_rich  <- "log(Species Richness)"
y_lab_abund <- "log(Abundance)"

# ---------------------------------------------------------------------
# Robust trophic mapping for points (fixes NA color issue)
# ---------------------------------------------------------------------
normalize_troph_code <- function(x) {
  v  <- as.character(x)
  n1 <- suppressWarnings(as.integer(v))
  digits <- suppressWarnings(as.integer(gsub("\\D+", "", v)))
  n2 <- ifelse(is.na(n1), digits, n1)
  vt <- tolower(trimws(v))
  n3 <- ifelse(is.na(n2) & grepl("primary",   vt), 2L, n2)
  n3 <- ifelse(is.na(n3) & grepl("second",    vt), 3L, n3)
  n3 <- ifelse(is.na(n3) & (grepl("top", vt) | grepl("carniv", vt) | grepl("pred", vt)), 4L, n3)
  ifelse(n3 %in% c(2L,3L,4L), n3, NA_integer_)
}

mk_points_rich <- function(df) {
  code <- normalize_troph_code(df$trophLevel)
  df |>
    mutate(
      group  = factor(paste0("richness_", code),
                      levels = c("richness_2","richness_3","richness_4")),
      y_plot = tlog_point_rich(richness)
    )
}
mk_points_cnt <- function(df) {
  code <- normalize_troph_code(df$trophLevel)
  df |>
    mutate(
      group  = factor(paste0("count_", code),
                      levels = c("count_2","count_3","count_4")),
      y_plot = tlog_point_cnt(count)
    )
}

s.sparea.troph  <- mk_points_rich(s.sparea.troph)
w.sparea.troph  <- mk_points_rich(w.sparea.troph)
s.spcount.troph <- mk_points_cnt(s.spcount.troph)
w.spcount.troph <- mk_points_cnt(w.spcount.troph)

# Optional sanity check
for (nm in c("s.sparea.troph","w.sparea.troph","s.spcount.troph","w.spcount.troph")) {
  df <- get(nm)
  if (any(is.na(df$group))) {
    message("⚠ NA trophic group in ", nm, " ; unmapped trophLevel values: ",
            paste(unique(df$trophLevel[is.na(df$group)]), collapse = ", "))
  }
}

# ---------------------------------------------------------------------
# Axes limits (include predictions + points under log transforms)
# ---------------------------------------------------------------------
xlim_rich  <- range(c(pred_rich_A.s$x,  pred_rich_A.w$x,
                      s.sparea.troph$logArea, w.sparea.troph$logArea), na.rm = TRUE)
xlim_abund <- range(c(pred_count_A.s$x, pred_count_A.w$x,
                      s.spcount.troph$logArea, w.spcount.troph$logArea), na.rm = TRUE)

ylim_rich  <- range(c(
  tlog_line(pred_rich_A.s$predicted),
  tlog_line(pred_rich_A.w$predicted),
  s.sparea.troph$y_plot, w.sparea.troph$y_plot
), na.rm = TRUE)

ylim_abund <- range(c(
  tlog_line(pred_count_A.s$predicted),
  tlog_line(pred_count_A.w$predicted),
  s.spcount.troph$y_plot, w.spcount.troph$y_plot
), na.rm = TRUE)

# ---------------------------------------------------------------------
# Plot styling
# ---------------------------------------------------------------------
okabe_ito <- c(
  "richness_2"="#009E73","richness_3"="#E69F00","richness_4"="#D55E00",
  "count_2"   ="#009E73","count_3"   ="#E69F00","count_4"   ="#D55E00"
)
custom_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11)
  )

# ---------------------------------------------------------------------
# Plots (AREA × trophic level): 2×2 (Summer/Winter × Richness/Abundance)
# ---------------------------------------------------------------------
richWVtroph.s <- ggplot(pred_rich_A.s, aes(x = x, y = tlog_line(predicted), color = group)) +
  annotate("text",
           x = min(pred_rich_A.s$x, na.rm = TRUE) + 0.1 * diff(range(pred_rich_A.s$x, na.rm = TRUE)),
           y = min(ylim_rich) + 0.95 * diff(range(ylim_rich)),
           label = "p>0.05", size = 3, hjust = 0) +
  # geom_point(data = s.sparea.troph, aes(x = logArea, y = y_plot, color = group),
  #            inherit.aes = FALSE, alpha = 0.45, size = 1.8, stroke = 0) +
  geom_line(linewidth = 1.1) +
  labs(title = "Summer", x = "log2(Wetland Area)", y = y_lab_rich, color = "Trophic Level") +
  scale_color_manual(values = okabe_ito, labels = trophrich_labels, drop = FALSE, na.translate = FALSE) +
  scale_x_continuous(limits = xlim_rich) + scale_y_continuous(limits = ylim_rich) +
  custom_theme

richWVtroph.w <- ggplot(pred_rich_A.w, aes(x = x, y = tlog_line(predicted), color = group)) +
  annotate("text",
           x = min(pred_rich_A.w$x, na.rm = TRUE) + 0.1 * diff(range(pred_rich_A.w$x, na.rm = TRUE)),
           y = min(ylim_rich) + 0.95 * diff(range(ylim_rich)),
           label = "p<0.01", size = 3, hjust = 0) +
  # geom_point(data = w.sparea.troph, aes(x = logArea, y = y_plot, color = group),
  #            inherit.aes = FALSE, alpha = 0.45, size = 1.8, stroke = 0) +
  geom_line(linewidth = 1.1) +
  labs(title = "Winter", x = "log2(Wetland Area)", y = NULL, color = "Trophic Level") +
  scale_color_manual(values = okabe_ito, labels = trophrich_labels, drop = FALSE, na.translate = FALSE) +
  scale_x_continuous(limits = xlim_rich) + scale_y_continuous(limits = ylim_rich) +
  custom_theme

countWVtroph.s <- ggplot(pred_count_A.s, aes(x = x, y = tlog_line(predicted), color = group)) +
  annotate("text",
           x = min(pred_count_A.s$x, na.rm = TRUE) + 0.1 * diff(range(pred_count_A.s$x, na.rm = TRUE)),
           y = min(ylim_abund) + 0.95 * diff(range(ylim_abund)),
           label = "p<0.001", size = 3, hjust = 0) +
  # geom_point(data = s.spcount.troph, aes(x = logArea, y = y_plot, color = group),
  #            inherit.aes = FALSE, alpha = 0.45, size = 1.8, stroke = 0) +
  geom_line(linewidth = 1.1) +
  labs(title = "Summer", x = "log2(Wetland Area)", y = y_lab_abund, color = "Trophic Level") +
  scale_color_manual(values = okabe_ito, labels = trophcount_labels, drop = FALSE, na.translate = FALSE) +
  scale_x_continuous(limits = xlim_abund) + scale_y_continuous(limits = ylim_abund) +
  custom_theme

countWVtroph.w <- ggplot(pred_count_A.w, aes(x = x, y = tlog_line(predicted), color = group)) +
  annotate("text",
           x = min(pred_count_A.w$x, na.rm = TRUE) + 0.1 * diff(range(pred_count_A.w$x, na.rm = TRUE)),
           y = min(ylim_abund) + 0.95 * diff(range(ylim_abund)),
           label = "p<0.01", size = 3, hjust = 0) +
  # geom_point(data = w.spcount.troph, aes(x = logArea, y = y_plot, color = group),
  #            inherit.aes = FALSE, alpha = 0.45, size = 1.8, stroke = 0) +
  geom_line(linewidth = 1.1) +
  labs(title = "Winter", x = "log2(Wetland Area)", y = NULL, color = "Trophic Level") +
  scale_color_manual(values = okabe_ito, labels = trophcount_labels, drop = FALSE, na.translate = FALSE) +
  scale_x_continuous(limits = xlim_abund) + scale_y_continuous(limits = ylim_abund) +
  custom_theme

combined_plot <- (
  richWVtroph.s + richWVtroph.w + plot_spacer() + countWVtroph.s + countWVtroph.w
) +
  plot_layout(ncol = 5, widths = c(1, 1, 0.05, 1, 1), guides = "collect") +
  plot_annotation(tag_levels = 'a') &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11)
  )

# Save figure to Figures/
ggsave(here("Figures", "troph_A_combined_color.png"),
       combined_plot, width = 1400/120, height = 500/120, dpi = 120, units = "in")

# ---------------------------------------------------------------------
# Save outputs: coefficients, CIs, ANOVA, predictions, figure, models
# ---------------------------------------------------------------------
extract_glmmTMB <- function(mod, label) {
  sm  <- summary(mod)
  ll  <- as.numeric(logLik(mod))
  aic <- AIC(mod); bic <- BIC(mod); n <- stats::nobs(mod)
  
  cf <- as.data.frame(sm$coefficients$cond)
  cf$term <- rownames(cf); rownames(cf) <- NULL
  names(cf) <- c("estimate","std_error","z_value","p_value","term")
  cf <- dplyr::select(cf, term, estimate, std_error, z_value, p_value); cf$model <- label
  
  ci <- tryCatch({
    ci_m <- suppressMessages(confint(mod, method = "wald"))
    out  <- as.data.frame(ci_m); out$term <- rownames(ci_m); rownames(ci_m) <- NULL
    names(out)[1:2] <- c("conf_low","conf_high"); out
  }, error = function(e) NULL)
  if (!is.null(ci)) cf <- dplyr::left_join(cf, ci, by = "term")
  
  av <- tryCatch({
    a <- car::Anova(mod, type = "III")
    at <- as.data.frame(a); at$term <- rownames(at); rownames(at) <- NULL
    names(at) <- sub("Pr\\(>Chisq\\)", "p_value", names(at)); at$model <- label; at
  }, error = function(e) NULL)
  
  fit <- tibble::tibble(
    model  = label,
    family = as.character(family(mod)$family),
    theta  = sm$sigma,
    logLik = ll, AIC = aic, BIC = bic, nobs = n
  )
  list(coefs = cf, anova = av, fit = fit)
}

OUT_ROOT <- here("results", "trophic_sar")
STAMP    <- format(Sys.Date(), "%Y-%m-%d")
OUT_DIR  <- file.path(OUT_ROOT, STAMP)
fs::dir_create(OUT_DIR, recurse = TRUE)

# Extract & write for all 8 models
mouts <- list(
  s_rich_A = extract_glmmTMB(troph_rich_A.s,  "summer_richness_logArea_by_troph"),
  s_rich_W = extract_glmmTMB(troph_rich_WV.s, "summer_richness_logWV_by_troph"),
  s_cnt_A  = extract_glmmTMB(troph_count_A.s, "summer_count_logArea_by_troph"),
  s_cnt_W  = extract_glmmTMB(troph_count_WV.s,"summer_count_logWV_by_troph"),
  w_rich_A = extract_glmmTMB(troph_rich_A.w,  "winter_richness_logArea_by_troph"),
  w_rich_W = extract_glmmTMB(troph_rich_WV.w, "winter_richness_logWV_by_troph"),
  w_cnt_A  = extract_glmmTMB(troph_count_A.w, "winter_count_logArea_by_troph"),
  w_cnt_W  = extract_glmmTMB(troph_count_WV.w,"winter_count_logWV_by_troph")
)

purrr::iwalk(mouts, function(x, nm) {
  write_csv(x$coefs, file.path(OUT_DIR, paste0(nm, "_coefficients.csv")))
  if (!is.null(x$anova)) write_csv(x$anova, file.path(OUT_DIR, paste0(nm, "_anova_typeIII.csv")))
  write_csv(x$fit,  file.path(OUT_DIR, paste0(nm, "_fitstats.csv")))
})

# Compact index for AREA×trophic models (interaction p-values)
interaction_p <- function(av) {
  if (is.null(av)) return(NA_real_)
  r <- av[grepl("logArea:.*troph", av$term), , drop = FALSE]
  if (!nrow(r)) NA_real_ else suppressWarnings(as.numeric(r$p_value[1]))
}
model_index <- dplyr::bind_rows(
  mouts$s_rich_A$fit, mouts$s_cnt_A$fit, mouts$w_rich_A$fit, mouts$w_cnt_A$fit
) |>
  dplyr::mutate(p_interaction = c(
    interaction_p(mouts$s_rich_A$anova),
    interaction_p(mouts$s_cnt_A$anova),
    interaction_p(mouts$w_rich_A$anova),
    interaction_p(mouts$w_cnt_A$anova)
  ))
write_csv(model_index, file.path(OUT_DIR, "model_index_area_by_trophic.csv"))

# Save prediction grids used in figure (original scale)
write_csv(pred_rich_A.s,  file.path(OUT_DIR, "pred_rich_A_summer.csv"))
write_csv(pred_rich_A.w,  file.path(OUT_DIR, "pred_rich_A_winter.csv"))
write_csv(pred_count_A.s, file.path(OUT_DIR, "pred_count_A_summer.csv"))
write_csv(pred_count_A.w, file.path(OUT_DIR, "pred_count_A_winter.csv"))

# Save figure here too
ggsave(file.path(OUT_DIR, "troph_A_combined_color.png"),
       combined_plot, width = 1400/120, height = 500/120, dpi = 120, units = "in")

# Save model objects
saveRDS(
  list(
    troph_rich_A_s = troph_rich_A.s,  troph_rich_WV_s = troph_rich_WV.s,
    troph_count_A_s = troph_count_A.s, troph_count_WV_s = troph_count_WV.s,
    troph_rich_A_w = troph_rich_A.w,  troph_rich_WV_w = troph_rich_WV.w,
    troph_count_A_w = troph_count_A.w, troph_count_WV_w = troph_count_WV.w
  ),
  file.path(OUT_DIR, "models_glmmTMB_trophicSAR.rds")
)

message("✅ Trophic-by-level SAR: results written to ", OUT_DIR)
