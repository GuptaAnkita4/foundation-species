# src/R/build_datasets.R
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
})

# -----------------------------
# 0) label maps (yours, kept)
# -----------------------------
wv_labels <- c("0" = "Low", "1" = "Medium", "2" = "High")
trophcount_labels <- c(
  "count_2" = "Primary consumers",
  "count_3" = "Secondary consumers",
  "count_4" = "Top carnivores"
)
trophrich_labels <- c(
  "richness_2" = "Primary consumers",
  "richness_3" = "Secondary consumers",
  "richness_4" = "Top carnivores"
)

# -----------------------------
# helpers
# -----------------------------
num <- function(x) suppressWarnings(as.numeric(x))
safe_log <- function(x) ifelse(is.na(x) | x <= 0, NA_real_, log(x))

area_scale <- function(x, mode = c("x10_plus1", "none")) {
  mode <- match.arg(mode)
  if (mode == "x10_plus1") return(x * 10 + 1)
  x
}

warn_missing_keys <- function(lhs_keys, rhs_keys, lhs_name = "left", rhs_name = "right") {
  miss_rhs <- setdiff(lhs_keys, rhs_keys)
  if (length(miss_rhs)) {
    message("⚠️ Keys in ", lhs_name, " not found in ", rhs_name, ": ",
            paste(utils::head(miss_rhs, 10), collapse = ", "),
            if (length(miss_rhs) > 10) " …")
  }
}

# -----------------------------
# 1) Basic species–area table
#    (replicates your s.sparea / w.sparea)
# -----------------------------
make_sparea <- function(habitat, indices, area_mode = c("x10_plus1","none")) {
  area_mode <- match.arg(area_mode)
  
  warn_missing_keys(habitat$grid, indices$grid, "habitat", "indices")
  
  out <- habitat %>%
    select(grid, wetlandTot, wetland, wetlandVegetation, perWV) %>%
    mutate(
      aT  = area_scale(num(wetlandTot), area_mode),
      aOW = area_scale(num(wetland), area_mode),
      aWV = area_scale(num(wetlandVegetation), area_mode),
      perWV = num(perWV)
    ) %>%
    select(grid, aT, aOW, aWV, perWV) %>%
    left_join(
      indices %>%
        select(grid, s = richness, a = count, d = shannonDiv),
      by = "grid"
    ) %>%
    mutate(
      s = num(s), a = num(a), d = num(d)
    ) %>%
    filter(!if_any(everything(), ~ is.na(.x)))
  
  # columns/ordering exactly like your original
  out %>% select(grid, aT, aOW, aWV, perWV, s, a, d)
}

# -----------------------------
# 2) Trophic distribution wide + long
#    (replicates w.troph.dist / s.troph.dist and long pivots)
# -----------------------------
make_trophic_dist <- function(habitat, indices) {
  warn_missing_keys(habitat$grid, indices$grid, "habitat", "indices")
  
  wide <- habitat %>%
    select(grid, area = wetlandTot, logArea, logWV, perWV, WVclass, change, changeSign) %>%
    mutate(
      area = num(area),
      logArea = num(logArea),
      logWV = num(logWV),
      perWV = num(perWV),
      WVclass = as.factor(WVclass),
      change = num(change),
      changeSign = as.factor(changeSign)
    ) %>%
    left_join(
      indices %>%
        select(
          grid,
          richness_2 = species_richness_2,
          richness_3 = species_richness_3,
          richness_4 = species_richness_4,
          count_2, count_3, count_4,
          MaxTrophicPosition, MeanTrophicLevel,
          biomass_2 = biomass2, biomass_3 = biomass3, biomass_4 = biomass4
        ),
      by = "grid"
    ) %>%
    mutate(
      across(starts_with("richness_"), num),
      across(starts_with("count_"), num),
      across(starts_with("biomass_"), num),
      MeanTrophicLevel = num(MeanTrophicLevel)
    ) %>%
    filter(!if_any(everything(), ~ is.na(.x)))
  
  long_rich <- wide %>%
    pivot_longer(starts_with("richness_"),
                 names_to = "trophLevel", values_to = "richness")
  long_count <- wide %>%
    pivot_longer(starts_with("count_"),
                 names_to = "trophLevel", values_to = "count")
  long_biomass <- wide %>%
    pivot_longer(starts_with("biomass_"),
                 names_to = "trophLevel", values_to = "biomass")
  
  list(wide = wide, richness = long_rich, count = long_count, biomass = long_biomass)
}

# -----------------------------
# 3) Area × vegetation interaction GLM data
#    (replicates w.glm.dat / s.glm.dat)
# -----------------------------
make_glm_dat <- function(habitat, indices, size_breaks = c(0, 2, Inf),
                         size_labels = c("small", "large")) {
  warn_missing_keys(habitat$grid, indices$grid, "habitat", "indices")
  
  habitat2 <- habitat %>%
    mutate(
      size_class = cut(num(wetlandTot), breaks = size_breaks,
                       labels = size_labels, right = TRUE, include.lowest = TRUE)
    )
  
  out <- habitat2 %>%
    select(grid, logArea, perWV, WVclass, size_class) %>%
    mutate(
      logArea = num(logArea),
      perWV = num(perWV),
      WVclass = as.factor(WVclass),
      size_class = as.factor(size_class)
    ) %>%
    left_join(
      indices %>%
        select(grid,
               richness,
               abundance = count,
               biomass = totalBiomass),
      by = "grid"
    ) %>%
    mutate(across(c(richness, abundance, biomass), num)) %>%
    filter(!if_any(everything(), ~ is.na(.x)))
  
  out %>% select(grid, logArea, perWV, WVclass, size_class, richness, abundance, biomass)
}

# -----------------------------
# 4) Functional diversity tables
#    (replicates w.f.dat / s.f.dat)
# -----------------------------
make_functional_dat <- function(habitat, indices) {
  warn_missing_keys(habitat$grid, indices$grid, "habitat", "indices")
  
  out <- habitat %>%
    select(grid, logArea, logOW, logWV, perWV, WVclass) %>%
    mutate(
      logArea = num(logArea),
      logOW = num(logOW),
      logWV = num(logWV),
      perWV = num(perWV),
      WVclass = as.factor(WVclass)
    ) %>%
    left_join(
      indices %>% select(grid, fdiv = FDiv, fric = FRic, feve = FEve),
      by = "grid"
    ) %>%
    mutate(
      across(c(fdiv, fric, feve), num),
      logFric = safe_log(fric)
    ) %>%
    filter(!if_any(everything(), ~ is.na(.x)))
  
  out %>% select(grid, logArea, logOW, logWV, perWV, WVclass, fdiv, fric, feve, logFric)
}

# -----------------------------
# 5) Factor label helpers for plotting
# -----------------------------
label_wvclass <- function(x) forcats::fct_recode(as.factor(x), !!!wv_labels)
label_troph_rich <- function(x) forcats::fct_recode(as.factor(x), !!!trophrich_labels)
label_troph_count <- function(x) forcats::fct_recode(as.factor(x), !!!trophcount_labels)
