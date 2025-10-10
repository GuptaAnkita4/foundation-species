# scripts/10_analysis_wv_interactions.R
suppressPackageStartupMessages({
  library(here)
  library(pacman)
})
pacman::p_load(
  dplyr, ggplot2, readr, glmmTMB, ggeffects, emmeans, car, DHARMa,
  patchwork, GGally, cowplot, sars
)

dir.create(here("Figures"), showWarnings = FALSE)

# ---- global epsilons for safe log on observed counts -----------------
EPS_RICH <- 0.5     # half-count for richness
EPS_ABUN <- 0.5     # half-count for abundance
tlog_point_rich <- function(y) log(pmax(y, EPS_RICH))
tlog_point_abun <- function(y) log(pmax(y, EPS_ABUN))
tlog_line       <- function(mu) log(mu)   # predictions (no epsilon)

# Ensure WVclass factor order and build log columns for plotting (pure log)
normalize_glm_dat <- function(df) {
  stopifnot(all(c("grid","logArea","perWV","WVclass","richness","abundance") %in% names(df)))
  df %>%
    mutate(
      WVclass  = factor(as.character(WVclass), levels = c("0","1","2")),
      log_rich = tlog_point_rich(richness),
      log_abun = tlog_point_abun(abundance)
    )
}
s.glm.dat <- normalize_glm_dat(s.glm.dat)
w.glm.dat <- normalize_glm_dat(w.glm.dat)

# quick pairplot (optional)
pairplot <- GGally::ggpairs(s.glm.dat[, c("logArea", "perWV", "richness", "abundance")])
ggsave(here("Figures", "pairplot_summer.png"), pairplot, width = 8, height = 6, dpi = 300)

# ---------------------- Modeling helpers ------------------------------------
fit_nb <- function(df, response) {
  form_full <- as.formula(paste0(response, " ~ logArea*factor(WVclass) + (1|grid)"))
  glmmTMB::glmmTMB(
    form_full, data = df, family = nbinom2,
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")),
    na.action = na.fail
  )
}

predict_grid <- function(mod, df, term_name = "logArea") {
  stopifnot(term_name %in% names(df), "WVclass" %in% names(df), "grid" %in% names(df))
  xs     <- seq(min(df[[term_name]], na.rm = TRUE),
                max(df[[term_name]], na.rm = TRUE),
                length.out = 100)
  groups <- levels(factor(df$WVclass, levels = c("0","1","2")))
  newdata <- expand.grid(logArea = xs, WVclass = groups, KEEP.OUT.ATTRS = FALSE)
  newdata$grid <- df$grid[1]  # dummy for excluding REs
  pr <- predict(mod, newdata = newdata, type = "response", se.fit = TRUE, re.form = NA)
  data.frame(
    x         = newdata[[term_name]],
    predicted = as.numeric(pr$fit),
    se        = as.numeric(pr$se.fit),
    group     = factor(newdata$WVclass, levels = c("0","1","2"))
  )
}

anova_sig_label <- function(mod) {
  a <- car::Anova(mod, type = "III")
  rn <- rownames(a); term <- which(grepl("logArea:factor\\(WVclass\\)", rn))
  p <- if (length(term)) a$`Pr(>Chisq)`[term] else NA_real_
  if (is.na(p)) "p = NA" else if (p < 0.001) "p < 0.001" else if (p < 0.01) "p < 0.01" else if (p < 0.05) "p < 0.05" else "p > 0.05"
}

diag_dharma <- function(mod, tag) {
  sim <- DHARMa::simulateResiduals(mod)
  plot(sim)
  ggplot2::ggsave(here("Figures", paste0("dharma_", tag, ".png")), width = 7, height = 5, dpi = 300)
  invisible(sim)
}

# ---------------------- Fit models ------------------------------------------
# Summer
s_nb_rich  <- fit_nb(s.glm.dat, "richness");  print(summary(s_nb_rich));  print(Anova(s_nb_rich, type = "III"))
s_nb_count <- fit_nb(s.glm.dat, "abundance"); print(summary(s_nb_count)); print(Anova(s_nb_count, type = "III"))

diag_dharma(s_nb_rich,  "summer_richness")
diag_dharma(s_nb_count, "summer_abundance")

pred_rich_s  <- predict_grid(s_nb_rich,  s.glm.dat)
pred_count_s <- predict_grid(s_nb_count, s.glm.dat)

lab_rich_s  <- anova_sig_label(s_nb_rich)
lab_count_s <- anova_sig_label(s_nb_count)

# Winter
w_nb_rich  <- fit_nb(w.glm.dat, "richness");  print(summary(w_nb_rich));  print(Anova(w_nb_rich, type = "III"))
w_nb_count <- fit_nb(w.glm.dat, "abundance"); print(summary(w_nb_count)); print(Anova(w_nb_count, type = "III"))

diag_dharma(w_nb_rich,  "winter_richness")
diag_dharma(w_nb_count, "winter_abundance")

pred_rich_w  <- predict_grid(w_nb_rich,  w.glm.dat)
pred_count_w <- predict_grid(w_nb_count, w.glm.dat)

lab_rich_w  <- anova_sig_label(w_nb_rich)
lab_count_w <- anova_sig_label(w_nb_count)

# ---------------------- Plotting (pure log everywhere) -----------------------
make_wv_plot <- function(pred_data, obs_data, season_label, yvar_logcol,
                         annotation_text, ylab_text, xlims = NULL, ylims = NULL) {
  
  pred_data <- as.data.frame(pred_data)
  pred_data$line_group <- factor(pred_data$group, levels = c("0","1","2"))
  obs_data$point_group <- factor(obs_data$WVclass, levels = c("0","1","2"))
  
  fill_vals  <- c("0" = "#009E73", "1" = "#E69F00", "2" = "#CC79A7")
  shape_vals <- c("0" = 21, "1" = 22, "2" = 24)
  line_vals  <- c("0" = "#007256", "1" = "#B66F00", "2" = "#994C8D")
  wv_labels  <- c("Low (<30%)", "Medium (30–70%)", "High (>70%)")
  
  if (is.null(xlims)) xlims <- range(c(obs_data$logArea, pred_data$x), na.rm = TRUE)
  if (is.null(ylims)) ylims <- range(obs_data[[yvar_logcol]], na.rm = TRUE)
  
  ggplot() +
    # POINTS (observed; already on log scale)
    # geom_jitter(
    #   data = obs_data,
    #   aes(x = logArea, y = .data[[yvar_logcol]], shape = point_group, fill = point_group),
    #   width = 0.1, alpha = 0.4, size = 2, stroke = 0.4, color = "black"
    # ) +
    # # LINES (predicted; log of the mean)
    geom_line(
      data = pred_data,
      aes(x = x, y = tlog_line(predicted), color = line_group, group = line_group),
      linewidth = 1, inherit.aes = FALSE
    ) +
    annotate("text",
             x = min(obs_data$logArea, na.rm = TRUE) + 0.5,
             y = min(ylims, na.rm = TRUE) + 0.95 * diff(range(ylims, na.rm = TRUE)),
             label = annotation_text, size = 4, hjust = 0
    ) +
    scale_shape_manual(values = shape_vals, labels = wv_labels, name = "Points: % Wetland Vegetation") +
    scale_fill_manual(values = fill_vals,  labels = wv_labels, name = "Points: % Wetland Vegetation") +
    scale_color_manual(values = line_vals, labels = wv_labels, name = "% Wetland Vegetation") +
    labs(title = season_label, x = "log2(Wetland Area)", y = ylab_text) +
    scale_x_continuous(limits = xlims) +
    scale_y_continuous(limits = ylims) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text  = element_text(size = 14),
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 12)
    ) +
    guides(
      color = guide_legend(override.aes = list(linewidth = 2), order = 1),
      fill  = guide_legend(override.aes = list(size = 2), order = 2),
      shape = guide_legend(override.aes = list(size = 2), order = 2)
    )
}

# axis limits using log(y) (safe via +epsilon in normalize_glm_dat)
xlim_rich  <- range(c(s.glm.dat$logArea, w.glm.dat$logArea), na.rm = TRUE)
ylim_rich  <- range(c(s.glm.dat$log_rich, w.glm.dat$log_rich), na.rm = TRUE)
xlim_abund <- xlim_rich
ylim_abund <- range(c(s.glm.dat$log_abun, w.glm.dat$log_abun), na.rm = TRUE)

p_a <- make_wv_plot(pred_rich_s,  s.glm.dat, "Summer", "log_rich", lab_rich_s,  "log(Richness)",  xlim_rich,  ylim_rich)
p_b <- make_wv_plot(pred_rich_w,  w.glm.dat, "Winter", "log_rich", lab_rich_w,  "log(Richness)",  xlim_rich,  ylim_rich) + labs(y = NULL)
p_c <- make_wv_plot(pred_count_s, s.glm.dat, "Summer", "log_abun", lab_count_s, "log(Abundance)", xlim_abund, ylim_abund)
p_d <- make_wv_plot(pred_count_w, w.glm.dat, "Winter", "log_abun", lab_count_w, "log(Abundance)", xlim_abund, ylim_abund) + labs(y = NULL)

combined_plot <- (
  (p_a + p_b + plot_spacer() + p_c + p_d) +
    plot_layout(ncol = 5, widths = c(1, 1, 0.05, 1, 1), guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 12)
    ) &
    guides(
      color = guide_legend(nrow = 1),
      fill  = guide_legend(nrow = 1),
      shape = guide_legend(nrow = 1)
    )
)

# Save outputs
ggsave(here("Figures", "SAR_WV_combined_color.png"), combined_plot, width = 14, height = 5, dpi = 600)
ggsave(here("Figures", "SAR_WV_summer.png"), p_a, width = 6, height = 5, dpi = 600)
ggsave(here("Figures", "SAR_WV_winter.png"), p_b, width = 6, height = 5, dpi = 600)

print(combined_plot)

# ---- SAVE MODEL OUTPUTS -----------------------------------------------------
suppressPackageStartupMessages({ library(fs); library(readr); library(dplyr) })

OUT_ROOT <- here::here("results", "models")
STAMP    <- format(Sys.Date(), "%Y-%m-%d")
OUT_DIR  <- file.path(OUT_ROOT, STAMP)
fs::dir_create(OUT_DIR, recurse = TRUE)

extract_glmmTMB <- function(mod, label) {
  sm  <- summary(mod)
  ll  <- as.numeric(logLik(mod))
  aic <- AIC(mod); bic <- BIC(mod); n <- stats::nobs(mod)
  
  fam_obj  <- family(mod)
  fam_name <- if (is.list(fam_obj) && !is.null(fam_obj$family)) fam_obj$family else as.character(fam_obj)
  
  cf <- as.data.frame(sm$coefficients$cond, stringsAsFactors = FALSE)
  cf$term <- rownames(cf); rownames(cf) <- NULL
  names(cf) <- c("estimate","std_error","z_value","p_value","term")
  cf <- dplyr::select(cf, term, estimate, std_error, z_value, p_value)
  cf$component <- "cond"; cf$model <- label
  
  ci <- tryCatch({
    ci_m <- suppressMessages(confint(mod, parm = rownames(sm$coefficients$cond), method = "wald"))
    ci_df <- as.data.frame(ci_m); ci_df$term <- rownames(ci_m); rownames(ci_m) <- NULL
    names(ci_df)[1:2] <- c("conf_low","conf_high"); ci_df
  }, error = function(e) NULL)
  if (!is.null(ci)) cf <- dplyr::left_join(cf, ci, by = "term")
  
  av <- tryCatch({
    a <- car::Anova(mod, type = "III")
    at <- as.data.frame(a); at$term <- rownames(at); rownames(at) <- NULL
    names(at) <- sub("Pr\\(>Chisq\\)", "p_value", names(at)); at$model <- label; at
  }, error = function(e) NULL)
  
  fit <- dplyr::tibble(
    model   = label,
    family  = fam_name,
    theta   = sm$sigma,
    logLik  = ll, AIC = aic, BIC = bic, nobs = n
  )
  list(coefs = cf, anova = av, fit = fit)
}

mouts <- list(
  s_rich  = extract_glmmTMB(s_nb_rich,  "summer_richness"),
  s_abun  = extract_glmmTMB(s_nb_count, "summer_abundance"),
  w_rich  = extract_glmmTMB(w_nb_rich,  "winter_richness"),
  w_abun  = extract_glmmTMB(w_nb_count, "winter_abundance")
)

purrr::iwalk(mouts, function(x, nm) {
  readr::write_csv(x$coefs, file.path(OUT_DIR, paste0(nm, "_coefficients.csv")))
  if (!is.null(x$anova)) readr::write_csv(x$anova, file.path(OUT_DIR, paste0(nm, "_anova_typeIII.csv")))
  readr::write_csv(x$fit,  file.path(OUT_DIR, paste0(nm, "_fitstats.csv")))
})

combine_interaction_p <- function(av, pattern = "logArea:factor\\(WVclass\\)") {
  if (is.null(av)) return(NA_real_)
  row <- av[grepl(pattern, av$term), , drop = FALSE]
  if (!nrow(row)) NA_real_ else as.numeric(row$p_value[1])
}
model_index <- dplyr::bind_rows(
  mouts$s_rich$fit, mouts$s_abun$fit, mouts$w_rich$fit, mouts$w_abun$fit
) |>
  dplyr::mutate(
    p_interaction = c(
      combine_interaction_p(mouts$s_rich$anova),
      combine_interaction_p(mouts$s_abun$anova),
      combine_interaction_p(mouts$w_rich$anova),
      combine_interaction_p(mouts$w_abun$anova)
    )
  )
readr::write_csv(model_index, file.path(OUT_DIR, "models_index.csv"))

preds <- list(
  pred_rich_s  = pred_rich_s,
  pred_count_s = pred_count_s,
  pred_rich_w  = pred_rich_w,
  pred_count_w = pred_count_w
)
purrr::iwalk(preds, ~ readr::write_csv(.x, file.path(OUT_DIR, paste0(.y, ".csv"))))

saveRDS(
  list(
    s_nb_rich  = s_nb_rich,
    s_nb_count = s_nb_count,
    w_nb_rich  = w_nb_rich,
    w_nb_count = w_nb_count
  ),
  file.path(OUT_DIR, "models_glmmTMB.rds")
)

message("✅ Model outputs written to: ", OUT_DIR)
