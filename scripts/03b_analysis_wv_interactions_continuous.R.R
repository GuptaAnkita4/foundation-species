# scripts/10_analysis_wv_interactions_continuous.R
suppressPackageStartupMessages({
  library(here); library(pacman)
})
pacman::p_load(
  dplyr, ggplot2, readr, glmmTMB, ggeffects, emmeans, car, DHARMa,
  patchwork, GGally, cowplot, fs, purrr
  # mgcv  # <- uncomment if you run the GAM sensitivity at bottom
)

dir.create(here("Figures"), showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Expect s.glm.dat and w.glm.dat already in memory with columns:
# grid, logArea, perWV (0-100), WVclass, richness, abundance
# This script uses perWV as CONTINUOUS. We scale it to 0-1 for stability.
# -----------------------------------------------------------------------------

EPS_RICH <- 0.5
EPS_ABUN <- 0.5
tlog_point_rich <- function(y) log(pmax(y, EPS_RICH))
tlog_point_abun <- function(y) log(pmax(y, EPS_ABUN))
tlog_line       <- function(mu) log(mu)

normalize_glm_dat_cont <- function(df) {
  stopifnot(all(c("grid","logArea","perWV","richness","abundance") %in% names(df)))
  df %>%
    mutate(
      perWV_sc = pmin(pmax(perWV/100, 0), 1),   # 0–1
      log_rich = tlog_point_rich(richness),
      log_abun = tlog_point_abun(abundance)
    )
}
s.glm.dat <- normalize_glm_dat_cont(s.glm.dat)
w.glm.dat <- normalize_glm_dat_cont(w.glm.dat)

# Optional quick pairs
pairplot <- GGally::ggpairs(s.glm.dat[, c("logArea", "perWV", "richness", "abundance")])
ggsave(here("Figures", "pairplot_summer_continuous.png"), pairplot, width = 8, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# Modeling helpers (CONTINUOUS cover)
# -----------------------------------------------------------------------------
fit_nb_cont <- function(df, response, quadratic = TRUE) {
  if (quadratic) {
    form_full <- as.formula(paste0(
      response, " ~ logArea * perWV_sc + logArea * I(perWV_sc^2) + (1|grid)"
    ))
  } else {
    form_full <- as.formula(paste0(
      response, " ~ logArea * perWV_sc + (1|grid)"
    ))
  }
  glmmTMB::glmmTMB(
    form_full, data = df, family = nbinom2,
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")),
    na.action = na.fail
  )
}

predict_grid_cont <- function(mod, df, term_name = "logArea",
                              cover_vals_pct = c(10, 30, 50, 70, 90)) {
  xs <- seq(min(df[[term_name]], na.rm = TRUE),
            max(df[[term_name]], na.rm = TRUE), length.out = 120)
  covers <- cover_vals_pct/100
  newdata <- expand.grid(logArea = xs, perWV_sc = covers, KEEP.OUT.ATTRS = FALSE)
  newdata$grid <- df$grid[1]  # dummy to exclude REs
  pr <- predict(mod, newdata = newdata, type = "response", se.fit = TRUE, re.form = NA)
  out <- data.frame(
    x         = newdata[[term_name]],
    predicted = as.numeric(pr$fit),
    se        = as.numeric(pr$se.fit),
    cover_pct = factor(paste0(round(newdata$perWV_sc*100), "%"),
                       levels = paste0(cover_vals_pct, "%"))
  )
  out
}

anova_sig_label_cont <- function(mod, quadratic = TRUE) {
  a <- car::Anova(mod, type = "III")
  rn <- rownames(a)
  if (quadratic) {
    p1 <- a$`Pr(>Chisq)`[grep("logArea:perWV_sc$", rn)]
    p2 <- a$`Pr(>Chisq)`[grep("logArea:I\\(perWV_sc\\^2\\)", rn)]
    fmt <- function(p) if (is.na(p)) NA_character_ else if (p<0.001) "p < 0.001" else
      if (p<0.01) "p < 0.01" else if (p<0.05) "p < 0.05" else "p > 0.05"
    paste0("Area×Cover: ", fmt(p1), "; Area×Cover^2: ", fmt(p2))
  } else {
    p1 <- a$`Pr(>Chisq)`[grep("logArea:perWV_sc", rn)]
    if (is.na(p1)) "p = NA" else if (p1<0.001) "p < 0.001" else
      if (p1<0.01) "p < 0.01" else if (p1<0.05) "p < 0.05" else "p > 0.05"
  }
}

diag_dharma <- function(mod, tag) {
  sim <- DHARMa::simulateResiduals(mod)
  plot(sim)
  ggplot2::ggsave(here("Figures", paste0("dharma_", tag, "_continuous.png")),
                  width = 7, height = 5, dpi = 300)
  invisible(sim)
}

# -----------------------------------------------------------------------------
# Fit models (continuous; quadratic enabled by default)
# -----------------------------------------------------------------------------
QUAD <- TRUE  # set FALSE to drop the quadratic term

# Summer
s_nb_rich_c  <- fit_nb_cont(s.glm.dat, "richness",  quadratic = QUAD)
s_nb_abun_c  <- fit_nb_cont(s.glm.dat, "abundance", quadratic = QUAD)
print(summary(s_nb_rich_c));  print(Anova(s_nb_rich_c, type = "III"))
print(summary(s_nb_abun_c));  print(Anova(s_nb_abun_c, type = "III"))

diag_dharma(s_nb_rich_c, "summer_richness_cont")
diag_dharma(s_nb_abun_c, "summer_abundance_cont")

pred_rich_s_c  <- predict_grid_cont(s_nb_rich_c,  s.glm.dat)
pred_abun_s_c  <- predict_grid_cont(s_nb_abun_c,  s.glm.dat)

lab_rich_s_c <- anova_sig_label_cont(s_nb_rich_c, quadratic = QUAD)
lab_abun_s_c <- anova_sig_label_cont(s_nb_abun_c, quadratic = QUAD)

# Winter
w_nb_rich_c  <- fit_nb_cont(w.glm.dat, "richness",  quadratic = QUAD)
w_nb_abun_c  <- fit_nb_cont(w.glm.dat, "abundance", quadratic = QUAD)
print(summary(w_nb_rich_c));  print(Anova(w_nb_rich_c, type = "III"))
print(summary(w_nb_abun_c));  print(Anova(w_nb_abun_c, type = "III"))

diag_dharma(w_nb_rich_c, "winter_richness_cont")
diag_dharma(w_nb_abun_c, "winter_abundance_cont")

pred_rich_w_c  <- predict_grid_cont(w_nb_rich_c,  w.glm.dat)
pred_abun_w_c  <- predict_grid_cont(w_nb_abun_c,  w.glm.dat)

lab_rich_w_c <- anova_sig_label_cont(w_nb_rich_c, quadratic = QUAD)
lab_abun_w_c <- anova_sig_label_cont(w_nb_abun_c, quadratic = QUAD)

# -----------------------------------------------------------------------------
# Plotting (lines = chosen cover percentages)
# -----------------------------------------------------------------------------
make_cont_plot <- function(pred_data, obs_data, season_label, yvar_logcol, annotation_text,
                           ylab_text, xlims = NULL, ylims = NULL) {
  if (is.null(xlims)) xlims <- range(c(obs_data$logArea, pred_data$x), na.rm = TRUE)
  if (is.null(ylims)) ylims <- range(obs_data[[yvar_logcol]], na.rm = TRUE)
  
  ggplot() +
    geom_line(data = pred_data,
              aes(x = x, y = tlog_line(predicted), color = cover_pct, group = cover_pct),
              linewidth = 1) +
    annotate("text",
             x = min(obs_data$logArea, na.rm = TRUE) + 0.5,
             y = min(ylims, na.rm = TRUE) + 0.95*diff(range(ylims, na.rm = TRUE)),
             label = annotation_text, size = 4, hjust = 0) +
    scale_color_brewer(palette = "Set1", name = "Vegetation cover") +
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
    )
}

xlim_rich  <- range(c(s.glm.dat$logArea, w.glm.dat$logArea), na.rm = TRUE)
ylim_rich  <- range(c(s.glm.dat$log_rich, w.glm.dat$log_rich), na.rm = TRUE)
xlim_abun  <- xlim_rich
ylim_abun  <- range(c(s.glm.dat$log_abun, w.glm.dat$log_abun), na.rm = TRUE)

p_a_c <- make_cont_plot(pred_rich_s_c, s.glm.dat, "Summer", "log_rich", lab_rich_s_c, "log(Richness)",  xlim_rich,  ylim_rich)
p_b_c <- make_cont_plot(pred_rich_w_c, w.glm.dat, "Winter", "log_rich", lab_rich_w_c, "log(Richness)",  xlim_rich,  ylim_rich) + labs(y = NULL)
p_c_c <- make_cont_plot(pred_abun_s_c, s.glm.dat, "Summer", "log_abun", lab_abun_s_c, "log(Abundance)", xlim_abun, ylim_abun)
p_d_c <- make_cont_plot(pred_abun_w_c, w.glm.dat, "Winter", "log_abun", lab_abun_w_c, "log(Abundance)", xlim_abun, ylim_abun) + labs(y = NULL)

combined_plot_c <- (
  (p_a_c + p_b_c + plot_spacer() + p_c_c + p_d_c) +
    plot_layout(ncol = 5, widths = c(1, 1, 0.05, 1, 1), guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 12)
    ) &
    guides(color = guide_legend(nrow = 1))
)

# Save figures
ggsave(here("Figures", "SAR_WV_continuous_combined.png"), combined_plot_c, width = 14, height = 5, dpi = 600)
ggsave(here("Figures", "SAR_WV_continuous_summer.png"), p_a_c, width = 6, height = 5, dpi = 600)
ggsave(here("Figures", "SAR_WV_continuous_winter.png"), p_b_c, width = 6, height = 5, dpi = 600)

print(combined_plot_c)

# -----------------------------------------------------------------------------
# SAVE MODEL OUTPUTS
# -----------------------------------------------------------------------------
OUT_ROOT <- here::here("results", "models_continuous")
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
    at <- as.data.frame(a); at$term <- rownames(a); rownames(a) <- NULL
    names(at) <- sub("Pr\\(>Chisq\\)", "p_value", names(at)); at$model <- label; at
  }, error = function(e) NULL)
  
  fit <- dplyr::tibble(
    model   = label, family = fam_name, theta = sm$sigma,
    logLik  = ll, AIC = aic, BIC = bic, nobs = n
  )
  list(coefs = cf, anova = av, fit = fit)
}

mouts <- list(
  s_rich_c  = extract_glmmTMB(s_nb_rich_c,  "summer_richness_cont"),
  s_abun_c  = extract_glmmTMB(s_nb_abun_c,  "summer_abundance_cont"),
  w_rich_c  = extract_glmmTMB(w_nb_rich_c,  "winter_richness_cont"),
  w_abun_c  = extract_glmmTMB(w_nb_abun_c,  "winter_abundance_cont")
)

purrr::iwalk(mouts, function(x, nm) {
  readr::write_csv(x$coefs, file.path(OUT_DIR, paste0(nm, "_coefficients.csv")))
  if (!is.null(x$anova)) readr::write_csv(x$anova, file.path(OUT_DIR, paste0(nm, "_anova_typeIII.csv")))
  readr::write_csv(x$fit,  file.path(OUT_DIR, paste0(nm, "_fitstats.csv")))
})

# Index with p-values for interactions
get_p <- function(av, pattern) {
  if (is.null(av)) return(NA_real_)
  row <- av[grepl(pattern, av$term), , drop = FALSE]
  if (!nrow(row)) NA_real_ else suppressWarnings(as.numeric(row$p_value[1]))
}
model_index <- dplyr::bind_rows(
  mouts$s_rich_c$fit, mouts$s_abun_c$fit, mouts$w_rich_c$fit, mouts$w_abun_c$fit
) |>
  dplyr::mutate(
    p_logA_x_cover   = c(
      get_p(mouts$s_rich_c$anova,  "logArea:perWV_sc$"),
      get_p(mouts$s_abun_c$anova,  "logArea:perWV_sc$"),
      get_p(mouts$w_rich_c$anova,  "logArea:perWV_sc$"),
      get_p(mouts$w_abun_c$anova,  "logArea:perWV_sc$")
    ),
    p_logA_x_cover2  = c(
      get_p(mouts$s_rich_c$anova,  "logArea:I\\(perWV_sc\\^2\\)"),
      get_p(mouts$s_abun_c$anova,  "logArea:I\\(perWV_sc\\^2\\)"),
      get_p(mouts$w_rich_c$anova,  "logArea:I\\(perWV_sc\\^2\\)"),
      get_p(mouts$w_abun_c$anova,  "logArea:I\\(perWV_sc\\^2\\)")
    )
  )
readr::write_csv(model_index, file.path(OUT_DIR, "models_index_continuous.csv"))

preds <- list(
  pred_rich_s_c  = pred_rich_s_c,
  pred_abun_s_c  = pred_abun_s_c,
  pred_rich_w_c  = pred_rich_w_c,
  pred_abun_w_c  = pred_abun_w_c
)
purrr::iwalk(preds, ~ readr::write_csv(.x, file.path(OUT_DIR, paste0(.y, ".csv"))))

saveRDS(
  list(
    s_nb_rich_c  = s_nb_rich_c,
    s_nb_abun_c  = s_nb_abun_c,
    w_nb_rich_c  = w_nb_rich_c,
    w_nb_abun_c  = w_nb_abun_c
  ),
  file.path(OUT_DIR, "models_glmmTMB_continuous.rds")
)

message("✅ Continuous-cover model outputs written to: ", OUT_DIR)

# -----------------------------------------------------------------------------
# (Optional) GAM sensitivity with smooth interaction surface
# -----------------------------------------------------------------------------
pacman::p_load(mgcv)
gam_rich_s <- mgcv::gam(richness ~ s(logArea, k=4) + s(perWV_sc, k=4) +
                          ti(logArea, perWV_sc, k=c(4,4)) + s(grid, bs="re"),
                        family = nb(), data = s.glm.dat)
summary(gam_rich_s); plot(gam_rich_s, pages=1)
ggeffects::ggpredict(gam_rich_s, terms = c("logArea","perWV_sc[0.1,0.3,0.5,0.7,0.9]"))
