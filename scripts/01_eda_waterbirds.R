# notebooks/EDA_waterbirds.R
suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
  source(here("src/R/eda_utils.R"))
})

# ---- params / output dir ----
OUT_DIR <- here("reports", "eda")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- 0) quick base plots (optional; keep if you find them useful) ----
# Prefer ggplot, but if you want to keep these, at least label axes:
# plot(w.habitat$wetland ~ w.habitat$wetlandVegetation,
#      xlab = "Winter wetland vegetation", ylab = "Winter wetland area")
# boxplot(s.habitat$wetlandTot ~ s.habitat$WVclass,
#         xlab = "WV class (summer)", ylab = "Wetland area (ha)")

# ---- 1) season data frames (safe joins, clean names) ----
s.temp <- make_season_df(s.habitat, s.indices)
w.temp <- make_season_df(w.habitat, w.indices)

# ---- 2) SAR plots (use log1p to handle zeros) ----
p1 <- plot_scatter_loess(
  s.temp, area, richness,
  "Summer: Richness vs Area", "Wetland area (ha)", "Waterbird richness"
)

p2 <- plot_scatter_loess(
  s.temp, logArea, log1p(richness),
  "Summer: Log–Log SAR (safe)", "Log wetland area", "log1p(richness)"
)

p3 <- plot_scatter_loess(
  w.temp, area, richness,
  "Winter: Richness vs Area", "Wetland area (ha)", "Waterbird richness"
)

p4 <- plot_scatter_loess(
  w.temp, logArea, log1p(richness),
  "Winter: Log–Log SAR (safe)", "Log wetland area", "log1p(richness)"
)

sar_grid <- (p1 | p2) / (p3 | p4)
ggsave(filename = file.path(OUT_DIR, "SAR_plots_combined.png"),
       plot = sar_grid, width = 10, height = 8, dpi = 300)

# ---- 3) species-level metadata (traits + residency) ----
# Expect: waterbirds (species table with Code, residentStat)
#         func.dat.diet (with Code, trophLevel)
species_info <- waterbirds %>%
  dplyr::select(Code, residentStat) %>%
  left_join(func.dat.diet %>% dplyr::select(Code, trophLevel), by = "Code") %>%
  mutate(
    trophLevel = as.factor(trophLevel),
    residentStat = as.factor(residentStat)
  )

# Sanity: ensure matrices’ species exist in species_info
summer_species <- colnames(s.birds.ab)[-1]  # drop 'grid'
winter_species <- colnames(w.birds.ab)[-1]
missing_traits <- setdiff(union(summer_species, winter_species), species_info$Code)
if (length(missing_traits)) {
  message("⚠️ Species missing in species_info: ", paste(missing_traits, collapse = ", "))
}

# ---- 4) counts of species by trophic level × residency × season ----
combined_info <- bind_rows(
  species_info %>% filter(Code %in% summer_species) %>% mutate(season = "Summer"),
  species_info %>% filter(Code %in% winter_species) %>% mutate(season = "Winter")
)

summary_counts <- combined_info %>%
  group_by(season, trophLevel, residentStat, .drop = FALSE) %>%
  summarise(n_species = dplyr::n(), .groups = "drop")

print(summary_counts)

p_counts <- ggplot(summary_counts,
                   aes(x = trophLevel, y = n_species, fill = residentStat)) +
  geom_col() +
  facet_wrap(~ season) +
  labs(x = "Trophic level", y = "Number of species", fill = "Status",
       title = "Species counts by trophic level and residency") +
  theme_min_clean()
ggsave(file.path(OUT_DIR, "species_counts_by_traits.png"),
       p_counts, width = 9, height = 5, dpi = 300)

# ---- 5) total abundance by trophic level × residency × season ----
summer_counts <- sum_abundance_by_traits(s.birds.ab, species_info, "Summer")
winter_counts <- sum_abundance_by_traits(w.birds.ab, species_info, "Winter")
combined_counts <- bind_rows(summer_counts, winter_counts)

print(combined_counts)

p_abund <- ggplot(combined_counts,
                  aes(x = trophLevel, y = total_abundance, fill = residentStat)) +
  geom_col() +
  facet_wrap(~ season) +
  labs(x = "Trophic level", y = "Total abundance (individuals)", fill = "Status",
       title = "Total abundance by trophic level and residency") +
  theme_min_clean()
ggsave(file.path(OUT_DIR, "abundance_by_traits.png"),
       p_abund, width = 9, height = 5, dpi = 300)

# ---- 6) (optional) Proportions of nonbreeders per trophic level within season ----
props <- combined_info %>%
  count(season, trophLevel, residentStat, name = "n_species") %>%
  group_by(season, trophLevel) %>%
  mutate(total = sum(n_species),
         prop_nonbreed = ifelse(residentStat == "Native nonbreeding",
                                n_species / total, NA_real_)) %>%
  filter(!is.na(prop_nonbreed)) %>%
  ungroup()

print(props)

p_props <- ggplot(props, aes(x = trophLevel, y = prop_nonbreed)) +
  geom_col() +
  facet_wrap(~ season) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Trophic level", y = "Nonbreeders (% of species)",
       title = "Share of nonbreeding species by trophic level") +
  theme_min_clean()
ggsave(file.path(OUT_DIR, "prop_nonbreeders_by_trophic.png"),
       p_props, width = 9, height = 5, dpi = 300)

message("✅ EDA complete. Figures in: ", OUT_DIR)
