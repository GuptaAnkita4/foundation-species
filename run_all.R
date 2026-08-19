# =============================================================================
# run_all.R -- Master pipeline for the foundation-species manuscript
# =============================================================================
#
# Runs every analysis/figure script in the dependency order their in-memory
# object handoffs require, in ONE shared R session (not subprocesses -- see
# "WHY ONE SESSION" below), and prints a pass/fail summary at the end.
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
# WHY ONE SESSION, NOT ONE-SUBPROCESS-PER-SCRIPT
#   Scripts 00, 02, 03, 03b, and 04 do NOT reload their own data -- they
#   assume s.glm.dat / w.glm.dat / s.sparea / w.sparea / s.sparea.troph /
#   s.troph / etc. are already sitting in memory from 00_create_datasets.R,
#   in ITS specific object shape. Running each script as an isolated
#   `Rscript` subprocess would lose that shared state and every one of
#   these would fail immediately. So this script source()s everything into
#   the SAME R session, in the order below, and relies on that ordering to
#   avoid a real object-collision hazard (next paragraph).
#
# THE OBJECT-COLLISION HAZARD THIS ORDERING AVOIDS
#   05 and 08 are each self-contained: they call load_all_data() themselves
#   and build their OWN s.glm.dat / w.glm.dat -- in two different shapes,
#   neither identical to 00's (05 skips size_class entirely; 08 uses
#   build_datasets.R's make_glm_dat() but with its *default* 2-class
#   small/large size_breaks, not 00's explicit size_breaks = c(0, 2, Inf)).
#   If either ran BEFORE 03 or 03b, it would silently overwrite the
#   s.glm.dat/w.glm.dat that 03/03b require (which need columns
#   grid/logArea/perWV/WVclass/richness/abundance in 00's shape) and both
#   would fail with a confusing "column not found" error instead of the
#   real cause. So 05/08 are ordered strictly AFTER 03 and 03b below, and
#   03/03b are never re-run afterward in the same session.
#   [UPDATE, run_all.R code-quality pass: 07 used to also rebuild its own
#   s.glm.dat/w.glm.dat here (a third shape) -- that was dead code (the
#   objects it built were never read again) and has been removed from
#   scripts/07_biomass_fseffects_supp.R, so 07 no longer poses any part of
#   this collision hazard. Left its position unchanged since nothing
#   required moving it.]
#
#   06 only needs s.troph/w.troph (list objects) or s.troph.dist/
#   w.troph.dist, which nothing else in this pipeline touches -- its
#   position relative to 05/08 doesn't matter, but it's kept in the
#   "Group A" block below for clarity since it does share 00's objects
#   directly (unlike 05/08, it has no source()/load_all_data() of its own).
#
# EXPLORATORY SCRIPTS (skipped by default; --include-exploratory to run)
#   03b_analysis_wv_interactions_continuous.R.R -- continuous+quadratic
#     vegetation-cover models; not referenced by any current figure/table.
#     Still keeps a (1|grid) random intercept that 03's models deliberately
#     dropped (see that script's own header note) -- if this is ever
#     promoted into a live figure/table, that needs revisiting first.
#   07_biomass_fseffects_supp.R -- Tweedie biomass models; confirmed (in
#     the script's own header) not referenced by any figure in the
#     manuscript or supplement.
#
# NOT INCLUDED -- run separately if needed
#   map_inset_bysize.R (Fig. 1 panel a, the area-coded site map) needs
#   shapefile/spatial habitat data that isn't part of this data pipeline.
#   It also isn't yet in this project's script/ or src/ folders -- run it
#   from wherever it currently lives, with its own inputs.
#
# FIXED, formerly "KNOWN FRAGILITY" -- run_all.R code-quality pass
#   claude/02b (Fig. S1) used to read cached sar_aar outputs from a path
#   hardcoded to the date they were generated (results/sar_aar/2026-07-30/),
#   not from whatever results/sar_aar/<date>/ folder THIS run's
#   02_SAR_AAR.R step had just written. scripts/02b_figureS1_render.R now
#   resolves the most recent results/sar_aar/<date>/ folder at runtime
#   instead (find_latest_sar_aar_dir()) -- confirmed this doesn't change
#   what Figure S1 shows (a live refit and the old frozen snapshot agree to
#   the precision the figure displays; see
#   claude/code-quality-pass-2026-08-18.md for the numeric check), and it
#   keeps Fig. S1 permanently in sync with Fig. 1 panels d-g instead of the
#   two silently drifting apart after every re-run.
#
# ENVIRONMENT NOTE
#   `sars` (needed by 02 and 02b) and `DHARMa` (needed by 03 and 03b) may
#   not be installed in every R environment -- both require a working
#   compiler toolchain for their dependencies. If a step needs one of these
#   and it isn't installed, install manually first:
#     install.packages(c("sars", "DHARMa"))
#   Every script here now loads packages via a plain library() call (the
#   pacman::p_load() this note used to describe -- which silently
#   auto-installed missing packages -- was removed from every script in
#   this pass; see claude/code-quality-pass-2026-08-18.md). A missing
#   package now fails immediately with R's normal "there is no package
#   called X" error, exactly like any other error a step can throw -- that
#   step is reported as FAILED below and the rest of the pipeline still
#   runs.
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
# NOTE (run_all.R code-quality pass): `.log` used to be a zero-row
# data.frame that every run_step() call grew via rbind(.log, data.frame(...))
# -- the classic "growing a data.frame in a loop" R anti-pattern. Harmless
# at this pipeline's scale (11 steps), but easy to fix cleanly: accumulate
# each step's result as a list in `.log_rows` instead (O(1) append) and
# only build the data.frame once, in print_summary(), via
# dplyr::bind_rows(). Also dropped `ok`, a variable run_step() set
# (TRUE/FALSE, mirroring `result`) but never actually read anywhere --
# `status` below is computed from `result` directly, not `ok`.
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
  # as.data.frame(): dplyr::bind_rows() on a list of lists returns a tibble,
  # whose print method truncates wide character columns (the "produces"
  # column, in particular) with "~" -- exactly the column this summary
  # exists to show in full. Converting back to a plain data.frame keeps the
  # untruncated print() behavior this script always had.
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
# GROUP A -- must run first, in this order, sharing one R session
# =============================================================================

# ---- 00: core data build (MUST run before anything else) -------------------
# Creates: dat, s.habitat/s.indices/w.habitat/w.indices, s.sparea/w.sparea,
# s.troph/w.troph (lists), s.sparea.troph/w.sparea.troph, s.spcount.troph/
# w.spcount.troph, s.spbio.troph/w.spbio.troph, s.glm.dat/w.glm.dat (shape:
# grid/logArea/perWV/WVclass/size_class/richness/abundance/biomass),
# s.f.dat/w.f.dat. Also writes Data/derived/<date>/{summer,winter}_trophic_metrics.csv.
run_step("00", "scripts/00_create_datasets.R", "Data/derived/<date>/*_trophic_metrics.csv", required = TRUE)

# ---- 01: EDA (independent of every figure/table; safe to skip) -------------
# Reloads its own data; harmlessly overwrites `dat`/`s.habitat`/etc. with
# identical values. Nothing downstream reads those names again.
run_step("01", "scripts/01_eda_waterbirds.R", "reports/eda/*.png", required = FALSE)

# ---- 02: SAR/AAR power-law fits -> Fig. 1 panels d-g ------------------------
# Needs s.sparea/w.sparea from 00. Requires the `sars` package (sar_power()).
run_step("02", "scripts/02_SAR_AAR.R", "Figures/SAR_AAR.png (Fig. 1 d-g)", required = FALSE)

# ---- 03: FS-cover x area GLMs -> Fig. 2 -------------------------------------
# Needs s.glm.dat/w.glm.dat in 00's shape (grid/logArea/perWV/WVclass/
# richness/abundance). Requires the `DHARMa` package for residual diagnostics.
run_step("03", "scripts/03_species-area-foundationspecies-models.R",
         "Figures/SAR_WV_combined_color.png (Fig. 2), results/models/<date>/", required = FALSE)

# ---- 03b: exploratory continuous-cover models (opt-in) ---------------------
# Not wired to any current figure/table. Must run before 05/08 rebuild
# s.glm.dat/w.glm.dat in a different shape, since it also needs 00's shape.
if (INCLUDE_EXPLORATORY) {
  run_step("03b", "scripts/03b_analysis_wv_interactions_continuous.R.R",
           "Figures/SAR_WV_continuous_*.png (exploratory, not in manuscript)", required = FALSE)
} else {
  cat("--- [03b] skipped (exploratory, not wired to any manuscript figure/table; pass --include-exploratory to run) ---\n\n")
}

# ---- 04: trophic-level x area/veg GLMMs -> Fig. 4 ---------------------------
# Needs s.sparea.troph/w.sparea.troph/s.spcount.troph/w.spcount.troph from 00.
run_step("04", "scripts/04_area-trophic-interaction.R",
         "Figures/troph_A_combined_color.png (Fig. 4), results/trophic_sar/<date>/", required = FALSE)

# ---- MASS/dplyr::select() masking guard -- REMOVED, run_all.R code-quality
# pass (verified safe; see claude/code-quality-pass-2026-08-18.md) --------
# A `select <- dplyr::select` .GlobalEnv rebind used to sit here. History,
# for anyone re-adding a script that loads MASS (directly or transitively):
# MASS exports its own S3 generic called `select()` (stepwise model
# selection -- unrelated to picking data-frame columns), and once MASS is
# attached ahead of dplyr in the search path, every later BARE `select(...)`
# call anywhere in this shared session resolves to MASS::select() instead
# of dplyr::select() -- re-calling library(dplyr) later does NOT fix this,
# because library() only changes search-path position on a package's FIRST
# attachment. This was originally found because scripts/04 called
# library(MASS) directly, and src/build_datasets.R's make_glm_dat()/
# make_sparea() both called bare (unnamespaced) select() internally, so 08
# and 02b -- which call those helpers after 04 had run -- both failed with
# a cryptic "unused arguments (grid, logArea, ...)" error.
#
# Both root causes are gone now, checked exhaustively rather than assumed:
# (1) scripts/04 no longer calls library(MASS) at all -- confirmed unused
#     and removed. Note this does NOT mean MASS never gets attached in this
#     session anymore: scripts/06 still needs `multcomp` for cld(), and
#     multcomp attaches MASS transitively via TH.data (confirmed via
#     search() after a fresh library(multcomp) call) -- so the masking
#     MECHANISM is still technically live during/after step 06.
# (2) What actually made that matter is gone, though: every bare select()
#     call anywhere in scripts/ and src/ was namespaced to dplyr::select()
#     during this pass (build_datasets.R, eda_utils.R, and every script
#     that called select() at all). Grepped scripts/ and src/ end to end
#     for any remaining bare select( -- excluding dplyr::/tidyr::-prefixed
#     calls and comment text -- and found none. With no bare select() call
#     left anywhere in the pipeline, MASS being on the search path no
#     longer has anything to shadow, so the .GlobalEnv rebind has nothing
#     left to protect. Verified by removing it and re-running the full
#     pipeline (`--include-exploratory`): 11/11 steps still pass, including
#     08 and 02b, the two steps this guard was originally written to fix.
# If a future edit reintroduces a bare select() call anywhere downstream of
# a script that (directly or transitively) attaches MASS, it will resurface
# as the same "unused arguments" error 08/02b hit originally -- the fix is
# to namespace that call as dplyr::select(), matching the rest of this
# codebase, not to re-add this rebind.

# ---- 06: herbivore/omnivore/carnivore presence-abundance -> Fig. S3 --------
# Needs s.troph/w.troph (or s.troph.dist/w.troph.dist) from 00.
run_step("06", "scripts/06_trophic-reponse-to-fs.R",
         "Figures/Herbivore_Presence_Abundance.png (Fig. S3), results/trophic/<date>/", required = FALSE)

# ---- 07: exploratory Tweedie biomass models (opt-in) ------------------------
# Not referenced by any figure/table (confirmed in the script's own header).
# Needs s.spbio.troph/w.spbio.troph from 00. [UPDATE, run_all.R code-quality
# pass: 07 used to also rebuild its own s.glm.dat/w.glm.dat here as a side
# effect (a third shape) -- that was dead code (never read again) and has
# been removed from the script, so 07 no longer touches s.glm.dat/w.glm.dat
# at all. Left in this position regardless; nothing required moving it.]
if (INCLUDE_EXPLORATORY) {
  run_step("07", "scripts/07_biomass_fseffects_supp.R",
           "Figures/biomass_tweedie_*.png (exploratory, not in manuscript)", required = FALSE)
} else {
  cat("--- [07] skipped (exploratory, not wired to any manuscript figure/table; pass --include-exploratory to run) ---\n\n")
}

# =============================================================================
# GROUP B -- self-contained (each reloads its own data); order-safe from here
# =============================================================================

# ---- 05: area-class x vegetation-cover hotspots -> Fig. 3 -------------------
run_step("05", "scripts/05_fseffects_areaclass.R",
         "Figures/hotspots_ci_row.png (Fig. 3), results/hotspots_area_veg/<date>/", required = FALSE)

# ---- 08: size x vegetation-class summary -> Fig. S2 -------------------------
run_step("08", "scripts/08_figureS2_size_veg_diversity.R",
         "Figures/FigureS2_size_veg_diversity.png (Fig. S2), results/figureS2_size_veg/<date>/", required = FALSE)

# ---- 02b: standalone SAR/AAR render -> Fig. S1 (see "FIXED" note above) ----
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