# =============================================================================
# run_all.R -- Master pipeline for the foundation-species manuscript
# =============================================================================
#
# Runs every analysis/figure script in the dependency order their in-memory
# object handoffs require
#
# USAGE
#   Rscript run_all.R
#   Rscript run_all.R --include-exploratory   # also runs 03b and 07 (see below)
#
#   Or from within an R/RStudio session, just: source("run_all.R")
#   (set INCLUDE_EXPLORATORY <- TRUE before sourcing to include 03b/07)
#
# WORKING DIRECTORY
#   Run this from the project root (the folder containing this file,
#   scripts/, src/, and data/). It calls here::i_am() below to confirm.
#
# WHAT THIS PRODUCES
#   Figures/*.png            -- all main-text and supplementary figures
#                                except Fig. 1 panel a (see "NOT INCLUDED")
#   results/<name>/<date>/   -- coefficients, ANOVA tables, prediction grids,
#                                and saved model objects behind every figure
#   Data/derived/<date>/     -- the two derived trophic-metrics CSVs
#

#
# EXPLORATORY SCRIPTS (skipped by default; --include-exploratory to run)
#   03b_analysis_wv_interactions_continuous.R
#   07_biomass_fseffects_supp.R
#
#
# ENVIRONMENT NOTE
#   `sars` (needed by 02 and 02b) and `DHARMa` (needed by 03 and 03b) may
#   not be installed in every R environment -- both require a working
#   compiler toolchain for their dependencies. If a step needs one of these
#   and it isn't installed, install manually first:
#     install.packages(c("sars", "DHARMa"))
# =============================================================================

suppressPackageStartupMessages(library(here))
here::i_am("run_all.R")

INCLUDE_EXPLORATORY <- isTRUE(get0("INCLUDE_EXPLORATORY", ifnotfound = FALSE)) ||
  ("--include-exploratory" %in% commandArgs(trailingOnly = TRUE))

cat(sprintf(
  "\n=== run_all.R starting | project root: %s | exploratory scripts: %s ===\n\n",
  here::here(), if (INCLUDE_EXPLORATORY) "INCLUDED" else "skipped (default)"
))

# ---- bookkeeping -----------------------------------------------------------

.log_rows <- list()

run_step <- function(step, script_path, produces, required = TRUE) {
  cat(sprintf("--- [%s] %s ---\n", step, script_path))
  t0 <- Sys.time()
  note <- ""
  result <- tryCatch({
    source(here::here(script_path), echo = FALSE, chdir = FALSE)
    TRUE
  }, error = function(e) {
    note <<- conditionMessage(e)
    FALSE
  })
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  status <- if (result) "OK" else if (required) "FAILED" else "FAILED (non-blocking)"
  .log_rows[[length(.log_rows) + 1]] <<- list(
    step = step, script = script_path, produces = produces,
    status = status, seconds = round(secs, 1), note = note
  )
  if (result) {
    cat(sprintf("    OK (%.1fs)\n\n", secs))
  } else {
    cat(sprintf("    FAILED (%.1fs): %s\n\n", secs, note))
    if (required) {
      cat("    This step's outputs are required in-memory by later steps in\n",
          "    this run -- stopping the pipeline here. Fix the error above and\n",
          "    re-run (or `source(\"run_all.R\")` again after fixing it interactively).\n\n", sep = "")
      print_summary()
      stop(sprintf("run_all.R halted at step [%s]: %s", step, note), call. = FALSE)
    }
  }
  invisible(result)
}

print_summary <- function() {
  .log <- as.data.frame(dplyr::bind_rows(.log_rows), stringsAsFactors = FALSE)
  cat("\n=== run_all.R summary ===\n")
  print(.log[, c("step", "produces", "status", "seconds")], row.names = FALSE)
  n_fail <- sum(grepl("FAILED", .log$status))
  cat(sprintf("\n%d step(s) run, %d failed.\n", nrow(.log), n_fail))
  if (n_fail > 0) {
    cat("\nFailure detail:\n")
    for (i in which(grepl("FAILED", .log$status))) {
      cat(sprintf("  [%s] %s\n      %s\n", .log$step[i], .log$script[i], .log$note[i]))
    }
  }
}

# =============================================================================
# Any step can be skipped, reordered, or run standalone via `Rscript
# scripts/NN_*.R` without affecting any other step.
# =============================================================================

# ---- 00: core data build + writes the derived trophic-metrics CSVs --------

run_step("00", "scripts/00_create_datasets.R", "Data/derived/<date>/*_trophic_metrics.csv", required = TRUE)

# ---- 01: EDA (independent of every figure/table; safe to skip) -------------
run_step("01", "scripts/01_eda_waterbirds.R", "reports/eda/*.png", required = FALSE)

# ---- 02: SAR/AAR power-law fits -> Fig. 1 panels d-g ------------------------
run_step("02", "scripts/02_SAR_AAR.R", "Figures/SAR_AAR.png (Fig. 1 d-g)", required = FALSE)

# ---- 03: FS-cover x area GLMs -> Fig. 2 -------------------------------------
run_step("03", "scripts/03_species-area-foundationspecies-models.R",
          "Figures/SAR_WV_combined_color.png (Fig. 2), results/models/<date>/", required = FALSE)

# ---- 03b: exploratory continuous-cover models (opt-in) ---------------------
if (INCLUDE_EXPLORATORY) {
  run_step("03b", "scripts/03b_analysis_wv_interactions_continuous.R.R",
            "Figures/SAR_WV_continuous_*.png (exploratory, not in manuscript)", required = FALSE)
} else {
  cat("--- [03b] skipped (exploratory, not wired to any manuscript figure/table; pass --include-exploratory to run) ---\n\n")
}

# ---- 04: trophic-level x area/veg GLMMs -> Fig. 4 ---------------------------
run_step("04", "scripts/04_area-trophic-interaction.R",
          "Figures/troph_A_combined_color.png (Fig. 4), results/trophic_sar/<date>/", required = FALSE)


# ---- 06: herbivore/omnivore/carnivore presence-abundance -> Fig. S3 --------
run_step("06", "scripts/06_trophic-reponse-to-fs.R",
          "Figures/Herbivore_Presence_Abundance.png (Fig. S3), results/trophic/<date>/", required = FALSE)

# ---- 07: exploratory Tweedie biomass models (opt-in) ------------------------
if (INCLUDE_EXPLORATORY) {
  run_step("07", "scripts/07_biomass_fseffects_supp.R",
            "Figures/biomass_tweedie_*.png (exploratory, not in manuscript)", required = FALSE)
} else {
  cat("--- [07] skipped (exploratory, not wired to any manuscript figure/table; pass --include-exploratory to run) ---\n\n")
}

# ---- 05: area-class x vegetation-cover hotspots -> Fig. 3 -------------------
run_step("05", "scripts/05_fseffects_areaclass.R",
          "Figures/hotspots_ci_row.png (Fig. 3), results/hotspots_area_veg/<date>/", required = FALSE)

# ---- 08: size x vegetation-class summary -> Fig. S2 -------------------------
run_step("08", "scripts/08_figureS2_size_veg_diversity.R",
          "Figures/FigureS2_size_veg_diversity.png (Fig. S2), results/figureS2_size_veg/<date>/", required = FALSE)

# ---- 02b: standalone SAR/AAR render -> Fig. S1 ----
run_step("02b", "scripts/02b_figureS1_render.R",
          "Figures/FigureS1_SAR_AAR.png (Fig. S1)", required = FALSE)

# =============================================================================
print_summary()

cat("\nNot run by this script (see header comments):\n")
cat("  - map_inset_bysize.R          (Fig. 1 panel a; needs spatial/shapefile data)\n")
if (!INCLUDE_EXPLORATORY) {
  cat("  - 03b (exploratory, continuous-cover models; --include-exploratory to run)\n")
  cat("  - 07  (exploratory, Tweedie biomass models; --include-exploratory to run)\n")
}
cat("\n=== run_all.R finished ===\n")
