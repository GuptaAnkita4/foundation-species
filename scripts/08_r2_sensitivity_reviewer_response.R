#!/usr/bin/env Rscript
# ==============================================================================
# Reviewer-requested analysis: R^2 reporting + grid-cell random-effect
# sensitivity + leverage check on smallest wetlands.
#
# Addresses Reviewer 1: General Comment 4 (site-selection / random-effect /
# leverage) and Specific Comment 9 (R^2 for GLMs/GLMMs, with vs without the
# grid-cell random factor).

# ==============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(glmmTMB)
  library(car); library(performance); library(insight); library(tibble)
})

options(width = 130)
source("src/load_data.R")
source("src/build_datasets.R")

dat <- load_all_data()
s.habitat <- dat$s.habitat; s.indices <- dat$s.indices
w.habitat <- dat$w.habitat; w.indices <- dat$w.indices

s.glm.dat <- make_glm_dat(s.habitat, s.indices, size_breaks = c(0, 2, Inf))
w.glm.dat <- make_glm_dat(w.habitat, w.indices, size_breaks = c(0, 2, Inf))

s.troph <- make_trophic_dist(s.habitat, s.indices)
w.troph <- make_trophic_dist(w.habitat, w.indices)
s.sparea.troph <- s.troph$richness;  s.spcount.troph <- s.troph$count
w.sparea.troph <- w.troph$richness;  w.spcount.troph <- w.troph$count

# ------------------------------------------------------------------------
# McFadden pseudo-R2: works for ANY model with a logLik method, with or
# without random effects, any glmmTMB family -- this is what lets us make a
# genuine apples-to-apples "with grid vs without grid" comparison, since
# performance::r2_nakagawa() only works when the RE is well identified and
# performance::r2() refuses nbinom2-without-RE glmmTMB objects outright.
# ------------------------------------------------------------------------
mcfadden_r2 <- function(mod, data, response, is_glmmTMB_re = FALSE) {
  null_form <- if (is_glmmTMB_re) as.formula(paste0(response, " ~ 1 + (1|grid)"))
               else as.formula(paste0(response, " ~ 1"))
  fam <- family(mod)
  if (is_glmmTMB_re || inherits(mod, "glmmTMB")) {
    null_mod <- glmmTMB(null_form, data = data, family = fam)
  } else {
    null_mod <- glm(null_form, data = data, family = fam)
  }
  1 - as.numeric(logLik(mod)) / as.numeric(logLik(null_mod))
}

get_interaction_p <- function(aov_tab, pattern) {
  rn <- rownames(aov_tab); idx <- grep(pattern, rn)
  if (!length(idx)) return(NA_real_)
  aov_tab$`Pr(>Chisq)`[idx[1]]
}

# ==========================================================================
# PART 1: FS-cover models (script 03, Fig 2) -- ONE ROW PER WETLAND
# ==========================================================================
cat("\n================ PART 1: FS-COVER MODELS (script 03 / Fig 2) ================\n")
cat("These tables have exactly one row per wetland (grid), confirmed below:\n")
cat("  summer: n =", nrow(s.glm.dat), ", n_distinct(grid) =", n_distinct(s.glm.dat$grid), "\n")
cat("  winter: n =", nrow(w.glm.dat), ", n_distinct(grid) =", n_distinct(w.glm.dat$grid), "\n")

run_fscover <- function(df, response, label) {
  form_re   <- as.formula(paste0(response, " ~ logArea*factor(WVclass) + (1|grid)"))
  form_nore <- as.formula(paste0(response, " ~ logArea*factor(WVclass)"))
  m_re   <- glmmTMB(form_re,   data = df, family = nbinom2,
                     control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
  m_nore <- glmmTMB(form_nore, data = df, family = nbinom2,
                     control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))

  aov_re   <- car::Anova(m_re,   type = "III")
  aov_nore <- car::Anova(m_nore, type = "III")

  theta_re   <- sigma(m_re)
  gridvar_re <- as.numeric(VarCorr(m_re)$cond$grid[1])

  tibble(
    model = label, n = nobs(m_re),
    logArea_slope_with_grid    = fixef(m_re)$cond["logArea"],
    logArea_slope_without_grid = fixef(m_nore)$cond["logArea"],
    AIC_with_grid = AIC(m_re), AIC_without_grid = AIC(m_nore),
    McFadden_R2_with_grid    = mcfadden_r2(m_re,   df, response, is_glmmTMB_re = TRUE),
    McFadden_R2_without_grid = mcfadden_r2(m_nore, df, response, is_glmmTMB_re = FALSE),
    grid_intercept_variance = gridvar_re,
    NB_theta_with_grid    = theta_re,
    NB_theta_without_grid = sigma(m_nore),
    p_interaction_with_grid    = get_interaction_p(aov_re,   "logArea:factor\\(WVclass\\)"),
    p_interaction_without_grid = get_interaction_p(aov_nore, "logArea:factor\\(WVclass\\)"),
    identifiable_grid_RE = n_distinct(df$grid) < nrow(df)  # FALSE => 1 row/grid => non-identifiable
  )
}

fscover_tab <- bind_rows(
  run_fscover(s.glm.dat, "richness",  "summer_richness_FScover"),
  run_fscover(s.glm.dat, "abundance", "summer_abundance_FScover"),
  run_fscover(w.glm.dat, "richness",  "winter_richness_FScover"),
  run_fscover(w.glm.dat, "abundance", "winter_abundance_FScover")
)
print(as.data.frame(fscover_tab))
cat("\nInterpretation: NB_theta_with_grid should blow up (>>1000) relative to NB_theta_without_grid\n")
cat("if the random intercept is non-identifiable (confounded with residual dispersion) --\n")
cat("i.e., one random-effect level per observation gives the model no repeated measures to pool.\n")

# ==========================================================================
# PART 2: Trophic-level models (script 04, Fig 4) -- THREE ROWS PER WETLAND
# ==========================================================================
cat("\n================ PART 2: TROPHIC-LEVEL MODELS (script 04 / Fig 4) ================\n")
cat("These long-format tables have THREE rows per wetland (one per trophic level):\n")
cat("  summer richness table: n =", nrow(s.sparea.troph), ", n_distinct(grid) =", n_distinct(s.sparea.troph$grid), "\n")
cat("  winter richness table: n =", nrow(w.sparea.troph), ", n_distinct(grid) =", n_distinct(w.sparea.troph$grid), "\n")

run_trophic <- function(df, response, label) {
  form_re   <- as.formula(paste0(response, " ~ logArea*factor(trophLevel) + (1|grid)"))
  form_nore <- as.formula(paste0(response, " ~ logArea*factor(trophLevel)"))
  m_re   <- glmmTMB(form_re,   data = df, family = nbinom2,
                     control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
  m_nore <- glmmTMB(form_nore, data = df, family = nbinom2,
                     control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
  r2n <- tryCatch(performance::r2_nakagawa(m_re), error = function(e) NULL)
  aov_re <- car::Anova(m_re, type = "III")

  tibble(
    model = label, n = nobs(m_re),
    AIC_with_grid = AIC(m_re), AIC_without_grid = AIC(m_nore),
    R2_marginal_with_grid    = if (!is.null(r2n)) as.numeric(r2n$R2_marginal) else NA,
    R2_conditional_with_grid = if (!is.null(r2n)) as.numeric(r2n$R2_conditional) else NA,
    McFadden_R2_without_grid = mcfadden_r2(m_nore, df, response, is_glmmTMB_re = FALSE),
    grid_intercept_variance = as.numeric(VarCorr(m_re)$cond$grid[1]),
    p_interaction_with_grid    = get_interaction_p(aov_re, "logArea:factor\\(trophLevel\\)"),
    p_interaction_without_grid = get_interaction_p(car::Anova(m_nore, type = "III"), "logArea:factor\\(trophLevel\\)")
  )
}

trophic_tab <- bind_rows(
  run_trophic(s.sparea.troph,  "richness", "summer_richness_trophic"),
  run_trophic(s.spcount.troph, "count",    "summer_count_trophic"),
  run_trophic(w.sparea.troph,  "richness", "winter_richness_trophic"),
  run_trophic(w.spcount.troph, "count",    "winter_count_trophic")
)
print(as.data.frame(trophic_tab))

# ==========================================================================
# PART 3: Herbivore/omnivore/carnivore presence + abundance (script 06, Fig S3)
# Original analysis NEVER includes a grid random effect. Confirm whether one
# row per grid also applies here (it should, same wide table structure).
# ==========================================================================
cat("\n================ PART 3: PRESENCE/ABUNDANCE-BY-COVER MODELS (script 06 / Fig S3) ================\n")
recode_wv <- function(x) {
  x <- as.character(x)
  look <- c("0"="Low","1"="Medium","2"="High","3"="High")
  factor(dplyr::recode(x, !!!look), levels = c("Low","Medium","High"))
}
compute_trophic_metrics <- function(df) {
  df %>%
    rowwise() %>%
    mutate(
      Herbivore_present = as.integer(count_2 > 0),
      Omnivore_present  = as.integer(count_3 > 0),
      Carnivore_present = as.integer(count_4 > 0)
    ) %>%
    ungroup() %>%
    mutate(WVclass = recode_wv(WVclass))
}
s.troph_metrics <- compute_trophic_metrics(s.troph$wide)
w.troph_metrics <- compute_trophic_metrics(w.troph$wide)
cat("summer: n =", nrow(s.troph_metrics), ", n_distinct(grid) =", n_distinct(s.troph_metrics$grid), "\n")
cat("winter: n =", nrow(w.troph_metrics), ", n_distinct(grid) =", n_distinct(w.troph_metrics$grid), "\n")

run_presence_abundance <- function(df, presence_col, count_col, label) {
  form_p_re   <- as.formula(paste0(presence_col, " ~ WVclass + logArea + (1|grid)"))
  form_p_nore <- as.formula(paste0(presence_col, " ~ WVclass + logArea"))
  glm_re   <- tryCatch(glmmTMB(form_p_re, data = df, family = binomial), error = function(e) NULL)
  glm_nore <- glm(form_p_nore, data = df, family = binomial)

  form_a_re   <- as.formula(paste0(count_col, " ~ WVclass + logArea + (1|grid)"))
  form_a_nore <- as.formula(paste0(count_col, " ~ WVclass + logArea"))
  nb_re   <- tryCatch(glmmTMB(form_a_re, data = df, family = nbinom2()), error = function(e) NULL)
  nb_nore <- glmmTMB(form_a_nore, data = df, family = nbinom2())

  tibble(
    model = label, n = nrow(df),
    presence_RE_fit_status = if (is.null(glm_re)) "failed to converge / non-identifiable" else "converged",
    presence_grid_var_if_fit = if (!is.null(glm_re)) tryCatch(as.numeric(VarCorr(glm_re)$cond$grid[1]), error=function(e) NA) else NA,
    presence_McFadden_R2_no_grid = 1 - as.numeric(logLik(glm_nore)) / as.numeric(logLik(glm(as.formula(paste0(presence_col," ~ 1")), data=df, family=binomial))),
    abundance_RE_fit_status = if (is.null(nb_re)) "failed to converge / non-identifiable" else "converged",
    abundance_grid_var_if_fit = if (!is.null(nb_re)) tryCatch(as.numeric(VarCorr(nb_re)$cond$grid[1]), error=function(e) NA) else NA,
    abundance_theta_no_grid = sigma(nb_nore),
    abundance_theta_with_grid_if_fit = if (!is.null(nb_re)) sigma(nb_re) else NA,
    abundance_McFadden_R2_no_grid = mcfadden_r2(nb_nore, df, count_col, is_glmmTMB_re = FALSE)
  )
}

presence_tab <- bind_rows(
  run_presence_abundance(s.troph_metrics, "Herbivore_present", "count_2", "summer_herbivore"),
  run_presence_abundance(s.troph_metrics, "Omnivore_present",  "count_3", "summer_omnivore"),
  run_presence_abundance(s.troph_metrics, "Carnivore_present", "count_4", "summer_carnivore"),
  run_presence_abundance(w.troph_metrics, "Herbivore_present", "count_2", "winter_herbivore"),
  run_presence_abundance(w.troph_metrics, "Omnivore_present",  "count_3", "winter_omnivore"),
  run_presence_abundance(w.troph_metrics, "Carnivore_present", "count_4", "winter_carnivore")
)
print(as.data.frame(presence_tab))

# ==========================================================================
# PART 4: Leverage / robustness check on the smallest wetlands (Fig 2 / Fig 1)
# Drop the k smallest-area wetlands one at a time (and cumulatively) and see
# how much the logArea x WVclass slope moves.
# ==========================================================================
cat("\n================ PART 4: LEVERAGE OF SMALLEST WETLANDS (FS-cover models) ================\n")
leverage_check <- function(df, response, label, k = 5) {
  ord <- order(df$logArea)
  smallest_ids <- df$grid[ord[1:k]]
  base_form <- as.formula(paste0(response, " ~ logArea*factor(WVclass) + (1|grid)"))
  m_full <- glmmTMB(base_form, data = df, family = nbinom2,
                     control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
  base_slope <- fixef(m_full)$cond["logArea"]

  drop_one <- purrr::map_dfr(seq_len(k), function(i) {
    df_sub <- df %>% filter(grid != smallest_ids[i])
    m_sub <- tryCatch(glmmTMB(base_form, data = df_sub, family = nbinom2,
                               control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))),
                       error = function(e) NULL)
    if (is.null(m_sub)) return(tibble(dropped_grid = smallest_ids[i], area_ha = NA, new_slope = NA, pct_change = NA))
    new_slope <- fixef(m_sub)$cond["logArea"]
    tibble(dropped_grid = smallest_ids[i],
           new_slope = new_slope,
           pct_change = 100 * (new_slope - base_slope) / base_slope)
  })

  df_sub_k <- df %>% filter(!grid %in% smallest_ids)
  m_dropk <- glmmTMB(base_form, data = df_sub_k, family = nbinom2,
                      control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
  slope_dropk <- fixef(m_dropk)$cond["logArea"]

  list(
    label = label, base_slope = base_slope, drop_one = drop_one,
    slope_dropping_k_smallest = slope_dropk,
    pct_change_dropping_k_smallest = 100 * (slope_dropk - base_slope) / base_slope
  )
}

lev_s_rich <- leverage_check(s.glm.dat, "richness",  "summer_richness")
lev_s_abun <- leverage_check(s.glm.dat, "abundance", "summer_abundance")
lev_w_rich <- leverage_check(w.glm.dat, "richness",  "winter_richness")
lev_w_abun <- leverage_check(w.glm.dat, "abundance", "winter_abundance")

for (lv in list(lev_s_rich, lev_s_abun, lev_w_rich, lev_w_abun)) {
  cat("\n--", lv$label, "-- base logArea slope =", round(lv$base_slope, 4), "\n")
  print(as.data.frame(lv$drop_one))
  cat("Dropping all 5 smallest wetlands at once: new slope =", round(lv$slope_dropping_k_smallest, 4),
      " (", round(lv$pct_change_dropping_k_smallest, 1), "% change)\n")
}

# ==========================================================================
# SAVE EVERYTHING
# ==========================================================================
dir.create("results/r2_sensitivity", recursive = TRUE, showWarnings = FALSE)
write_csv(fscover_tab,  "results/r2_sensitivity/fscover_models_r2.csv")
write_csv(trophic_tab,  "results/r2_sensitivity/trophic_models_r2.csv")
write_csv(presence_tab, "results/r2_sensitivity/presence_abundance_models_r2.csv")

leverage_out <- bind_rows(
  lev_s_rich$drop_one %>% mutate(model = lev_s_rich$label, base_slope = lev_s_rich$base_slope),
  lev_s_abun$drop_one %>% mutate(model = lev_s_abun$label, base_slope = lev_s_abun$base_slope),
  lev_w_rich$drop_one %>% mutate(model = lev_w_rich$label, base_slope = lev_w_rich$base_slope),
  lev_w_abun$drop_one %>% mutate(model = lev_w_abun$label, base_slope = lev_w_abun$base_slope)
)
write_csv(leverage_out, "results/r2_sensitivity/leverage_smallest_wetlands.csv")

cat("\n✅ All results written to results/r2_sensitivity/\n")
