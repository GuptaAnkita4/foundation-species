source(here::here("src/build_datasets.R"))
source(here::here("src/load_data.R"))

dat <- load_all_data()

s.habitat <- dat$s.habitat
s.indices <- dat$s.indices

w.habitat <- dat$w.habitat
w.indices <- dat$w.indices

# 1) species–area
s.sparea <- make_sparea(s.habitat, s.indices, area_mode = "x10_plus1")
w.sparea <- make_sparea(w.habitat, w.indices, area_mode = "x10_plus1")

# 2) trophic distributions
s.troph <- make_trophic_dist(s.habitat, s.indices)
w.troph <- make_trophic_dist(w.habitat, w.indices)

# long tables (names match your usage)
s.sparea.troph <- s.troph$richness
s.spcount.troph <- s.troph$count
s.spbio.troph  <- s.troph$biomass

w.sparea.troph <- w.troph$richness
w.spcount.troph <- w.troph$count
w.spbio.troph  <- w.troph$biomass

# 3) glm data with size classes (bin at 2 ha; adjust as needed)
s.glm.dat <- make_glm_dat(s.habitat, s.indices, size_breaks = c(0, 2, Inf))
w.glm.dat <- make_glm_dat(w.habitat, w.indices, size_breaks = c(0, 2, Inf))

# 4) functional diversity tables
s.f.dat <- make_functional_dat(s.habitat, s.indices)
w.f.dat <- make_functional_dat(w.habitat, w.indices)

# optional: relabel for plots on the fly
# s.glm.dat$WVclass <- label_wvclass(s.glm.dat$WVclass)
# w.sparea.troph$trophLevel <- label_troph_rich(w.sparea.troph$trophLevel)
# w.spcount.troph$trophLevel <- label_troph_count(w.spcount.troph$trophLevel)


# save these datasets for future use
suppressPackageStartupMessages({
  library(readr)
  library(fs)
  library(jsonlite)
})

# 1) choose an output folder (dated for versioning)
OUT_ROOT <- "Data/derived"
STAMP    <- format(Sys.Date(), "%Y-%m-%d")  # e.g., 2025-09-18
OUT_DIR  <- file.path(OUT_ROOT, STAMP)
dir_create(OUT_DIR, recurse = TRUE)

# tiny helpers
write_csv_safe <- function(x, path, ...) {
  # ensure data.frame/tibble + no rownames
  if (!inherits(x, c("data.frame", "tbl_df"))) x <- as.data.frame(x)
  write_csv(x, path, na = "", ...)
}

write_matrix_wide <- function(mat_or_df, path, rowname_as = "grid") {
  df <- mat_or_df
  # if object has rownames but no explicit grid column, promote rownames
  rn <- rownames(df)
  if (!is.null(rn) && !(rowname_as %in% names(df))) {
    df <- cbind(!!rowname_as := rn, df)
  }
  # ensure plain data.frame for readr
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  # put id column first
  nm <- names(df); nm <- c(rowname_as, setdiff(nm, rowname_as)); df <- df[, nm, drop = FALSE]
  write_csv_safe(df, path)
}

# 2) files to write (adjust to taste)
files <- list(
  s_sparea_troph = list(obj = s.sparea.troph,  file = file.path(OUT_DIR, "summer_trophic_metrics.csv")),
  w_sparea_troph = list(obj = w.sparea.troph,  file = file.path(OUT_DIR, "winter_trophic_metrics.csv"))
 )

# 3) write CSVs
purrr::walk(files, function(x) {
  obj <- x$obj; fp <- x$file
  if (grepl("birds_.*_csv$", names(files)[which(vapply(files, identical, logical(1), x))])) {
    # for community matrices keep first column as 'grid'
    write_matrix_wide(obj, fp, rowname_as = "grid")
  } else {
    write_csv_safe(obj, fp)
  }
})


# 5) simple manifest (what was written + basic dims)
manifest <- tibble::tibble(
  file = basename(vapply(files, `[[`, "", "file")),
  path = normalizePath(vapply(files, `[[`, "", "file"), winslash = "/"),
  nrow = vapply(files, function(x) nrow(as.data.frame(x$obj)), integer(1)),
  ncol = vapply(files, function(x) ncol(as.data.frame(x$obj)), integer(1))
)
write_json(manifest, file.path(OUT_DIR, "manifest.json"), pretty = TRUE, auto_unbox = TRUE)

message("✅ wrote derived datasets to: ", OUT_DIR)



