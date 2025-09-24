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
