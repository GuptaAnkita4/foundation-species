
suppressPackageStartupMessages({
  library(glmmTMB)
  library(car)
  library(ggeffects)
  library(ggplot2)
  library(here)
  library(fs)
  library(readr)
})

source(here::here("src/load_data.R"))
source(here::here("src/build_datasets.R"))

dat <- load_all_data()
s.spbio.troph <- make_trophic_dist(dat$s.habitat, dat$s.indices)$biomass
w.spbio.troph <- make_trophic_dist(dat$w.habitat, dat$w.indices)$biomass

# Output folders
fs::dir_create(here::here("Figures"))
STAMP   <- format(Sys.Date(), "%Y-%m-%d")
OUT_DIR <- here::here("results", "biomass_area_wv", STAMP)
fs::dir_create(OUT_DIR, recurse = TRUE)


# SUMMER
m_biomass_s_tw <- glmmTMB(
  biomass ~ logArea * trophLevel + (1|grid),
  data   = s.spbio.troph,
  family = tweedie(link = "log")   # estimates variance power p (1<p<2)
)
aov_biomass_s_tw <- car::Anova(m_biomass_s_tw, type = "III")

pred_biomass_s_tw <- ggeffects::ggpredict(
  m_biomass_s_tw, terms = c("logArea","trophLevel")
)

# WINTER
m_biomass_w_tw <- glmmTMB(
  biomass ~ logArea * trophLevel + (1|grid),
  data   = w.spbio.troph,
  family = tweedie(link = "log")
)
aov_biomass_w_tw <- car::Anova(m_biomass_w_tw, type = "III")

pred_biomass_w_tw <- ggeffects::ggpredict(
  m_biomass_w_tw, terms = c("logArea","trophLevel")
)

# ---- Plots------------------------
p_bio_s <- ggplot(as.data.frame(pred_biomass_s_tw),
                  aes(x = x, y = log(predicted), color = group)) +
  geom_line(linewidth = 1) +
  labs(x = "log(Area)", y = "log(Predicted biomass)", color = "Trophic level") +
  theme_classic()

p_bio_w <- ggplot(as.data.frame(pred_biomass_w_tw),
                  aes(x = x, y = log(predicted), color = group)) +
  geom_line(linewidth = 1) +
  labs(x = "log(Area)", y = "log(Predicted biomass)", color = "Trophic level") +
  theme_classic()

# ---- Save outputs -----------------------------------------------------------
readr::write_csv(
  tibble::rownames_to_column(as.data.frame(summary(m_biomass_s_tw)$coefficients$cond), "term"),
  file.path(OUT_DIR, "coef_biomass_tweedie_summer.csv")
)
readr::write_csv(
  tibble::rownames_to_column(as.data.frame(summary(m_biomass_w_tw)$coefficients$cond), "term"),
  file.path(OUT_DIR, "coef_biomass_tweedie_winter.csv")
)
readr::write_csv(
  tibble::rownames_to_column(as.data.frame(aov_biomass_s_tw), "term"),
  file.path(OUT_DIR, "anova_biomass_tweedie_summer.csv")
)
readr::write_csv(
  tibble::rownames_to_column(as.data.frame(aov_biomass_w_tw), "term"),
  file.path(OUT_DIR, "anova_biomass_tweedie_winter.csv")
)

ggsave(here::here("Figures", "biomass_tweedie_summer.png"), p_bio_s, width = 5, height = 4, dpi = 300)
ggsave(here::here("Figures", "biomass_tweedie_winter.png"), p_bio_w, width = 5, height = 4, dpi = 300)

saveRDS(
  list(tweedie = list(summer = m_biomass_s_tw, winter = m_biomass_w_tw)),
  file.path(OUT_DIR, "biomass_models_tweedie.rds")
)
