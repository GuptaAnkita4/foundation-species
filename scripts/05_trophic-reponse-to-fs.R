# scripts/20_trophic_analysis.R
suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(emmeans)
  library(glmmTMB)
  library(multcomp)
  library(multcompView)
  library(patchwork)
  library(vegan)   # for diversity()
  library(readr)
  library(fs)
})

dir.create(here("Figures"), showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Inputs: expect wide trophic tables:
#   - s.troph.dist, w.troph.dist  (columns: grid, logArea, WVclass, richness_2/3/4, count_2/3/4, ...)
# If you used our earlier helpers (make_trophic_dist), they live at s.troph$wide / w.troph$wide.
# This makes it compatible with either.
# -----------------------------------------------------------------------------
if (exists("s.troph") && is.list(s.troph) && "wide" %in% names(s.troph)) s.troph.dist <- s.troph$wide
if (exists("w.troph") && is.list(w.troph) && "wide" %in% names(w.troph)) w.troph.dist <- w.troph$wide

stopifnot(all(c("grid","logArea","WVclass","richness_2","richness_3","richness_4",
                "count_2","count_3","count_4") %in% names(s.troph.dist)))
stopifnot(all(c("grid","logArea","WVclass","richness_2","richness_3","richness_4",
                "count_2","count_3","count_4") %in% names(w.troph.dist)))

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
recode_wv <- function(x) {
  # map 0/1/2 or 1/2/3 -> Low/Medium/High
  x <- as.character(x)
  look <- c("0"="Low","1"="Medium","2"="High","3"="High")
  factor(dplyr::recode(x, !!!look), levels = c("Low","Medium","High"))
}

compute_trophic_metrics <- function(df) {
  df %>%
    rowwise() %>%
    mutate(
      TotalRichness       = richness_2 + richness_3 + richness_4,
      TotalAbundance      = count_2 + count_3 + count_4,
      trophicRichness     = sum(c_across(c(richness_2, richness_3, richness_4)) > 0),
      shannonRichnessTroph= vegan::diversity(c(richness_2, richness_3, richness_4), index = "shannon"),
      shannonAbundanceTroph= vegan::diversity(c(count_2, count_3, count_4), index = "shannon"),
      PropRich_2          = ifelse(TotalRichness > 0, richness_2/TotalRichness, NA_real_),
      PropRich_3          = ifelse(TotalRichness > 0, richness_3/TotalRichness, NA_real_),
      PropRich_4          = ifelse(TotalRichness > 0, richness_4/TotalRichness, NA_real_),
      PropAbund_2         = ifelse(TotalAbundance > 0, count_2/TotalAbundance, NA_real_),
      PropAbund_3         = ifelse(TotalAbundance > 0, count_3/TotalAbundance, NA_real_),
      PropAbund_4         = ifelse(TotalAbundance > 0, count_4/TotalAbundance, NA_real_),
      MeanTrophicLevel    = ifelse(TotalAbundance > 0, (2*count_2 + 3*count_3 + 4*count_4)/TotalAbundance, NA_real_),
      Herbivore_present   = as.integer(count_2 > 0),
      Omnivore_present    = as.integer(count_3 > 0),
      Carnivore_present   = as.integer(count_4 > 0)
    ) %>%
    ungroup() %>%
    mutate(WVclass = recode_wv(WVclass))
}

tidy_glm <- function(mod, label) {
  sm <- summary(mod)
  cf <- as.data.frame(sm$coefficients)
  cf$term <- rownames(cf); rownames(cf) <- NULL
  names(cf) <- c("estimate","std_error","z_value","p_value","term")
  ci <- suppressMessages(as.data.frame(confint(mod)))
  if (nrow(ci)) { ci$term <- rownames(ci); rownames(ci) <- NULL; names(ci)[1:2] <- c("conf_low","conf_high") }
  out <- dplyr::left_join(cf, ci, by = "term")
  out$model <- label
  out
}

tidy_tmb <- function(mod, label) {
  sm <- summary(mod)
  cf <- as.data.frame(sm$coefficients$cond)
  cf$term <- rownames(cf); rownames(cf) <- NULL
  names(cf) <- c("estimate","std_error","z_value","p_value","term")
  ci <- tryCatch({
    ci_m <- suppressMessages(confint(mod, method = "wald"))
    ci_df <- as.data.frame(ci_m); ci_df$term <- rownames(ci_m); rownames(ci_m) <- NULL
    names(ci_df)[1:2] <- c("conf_low","conf_high"); ci_df
  }, error = function(e) data.frame())
  out <- dplyr::left_join(cf, ci, by = "term")
  out$model <- label
  out
}

# Colors
veg_colors <- c("Low"="#A8E6A3","Medium"="#4CAF50","High"="#2E7D32")

# -----------------------------------------------------------------------------
# Build metrics
# -----------------------------------------------------------------------------
s.troph_metrics <- compute_trophic_metrics(s.troph.dist)
w.troph_metrics <- compute_trophic_metrics(w.troph.dist)

# -----------------------------------------------------------------------------
# Models
# -----------------------------------------------------------------------------
# Presence (binomial GLMs)
glm_herb_s <- glm(Herbivore_present ~ WVclass + logArea, family = binomial, data = s.troph_metrics)
glm_omn_s  <- glm(Omnivore_present  ~ WVclass + logArea, family = binomial, data = s.troph_metrics)
glm_carn_s <- glm(Carnivore_present ~ WVclass + logArea, family = binomial, data = s.troph_metrics)

glm_herb_w <- glm(Herbivore_present ~ WVclass + logArea, family = binomial, data = w.troph_metrics)
glm_omn_w  <- glm(Omnivore_present  ~ WVclass + logArea, family = binomial, data = w.troph_metrics)
glm_carn_w <- glm(Carnivore_present ~ WVclass + logArea, family = binomial, data = w.troph_metrics)

# Abundance (negative binomial)
nb_herb_s <- glmmTMB(count_2 ~ WVclass + logArea, family = nbinom2(), data = s.troph_metrics)
nb_omn_s  <- glmmTMB(count_3 ~ WVclass + logArea, family = nbinom2(), data = s.troph_metrics)
nb_carn_s <- glmmTMB(count_4 ~ WVclass + logArea, family = nbinom2(), data = s.troph_metrics)

nb_herb_w <- glmmTMB(count_2 ~ WVclass + logArea, family = nbinom2(), data = w.troph_metrics)
nb_omn_w  <- glmmTMB(count_3 ~ WVclass + logArea, family = nbinom2(), data = w.troph_metrics)
nb_carn_w <- glmmTMB(count_4 ~ WVclass + logArea, family = nbinom2(), data = w.troph_metrics)

# -----------------------------------------------------------------------------
# Predictions & CLDs
# -----------------------------------------------------------------------------
median_logArea_s <- median(s.troph_metrics$logArea, na.rm = TRUE)
median_logArea_w <- median(w.troph_metrics$logArea, na.rm = TRUE)

newdata_s <- data.frame(WVclass = factor(c("Low","Medium","High"), levels = c("Low","Medium","High")),
                        logArea = median_logArea_s)
newdata_w <- data.frame(WVclass = factor(c("Low","Medium","High"), levels = c("Low","Medium","High")),
                        logArea = median_logArea_w)

# Presence preds + CLDs (herbivores shown; replicate for omn/carn if desired)
newdata_s$predicted_prob <- predict(glm_herb_s, newdata = newdata_s, type = "response")
newdata_w$predicted_prob <- predict(glm_herb_w, newdata = newdata_w, type = "response")
cld_herb_s <- multcomp::cld(emmeans(glm_herb_s, ~ WVclass), adjust = "tukey", Letters = letters)
cld_herb_w <- multcomp::cld(emmeans(glm_herb_w, ~ WVclass), adjust = "tukey", Letters = letters)
cld_herb_s$WVclass <- as.character(cld_herb_s$WVclass)
cld_herb_w$WVclass <- as.character(cld_herb_w$WVclass)
annot_s <- merge(newdata_s, cld_herb_s[, c("WVclass",".group")], by = "WVclass")
annot_w <- merge(newdata_w, cld_herb_w[, c("WVclass",".group")], by = "WVclass")

# Abundance preds + CLDs (herbivores shown)
pred_s <- predict(nb_herb_s, newdata = newdata_s, type = "response", se.fit = TRUE)
pred_w <- predict(nb_herb_w, newdata = newdata_w, type = "response", se.fit = TRUE)
newdata_s$predicted_abundance <- pred_s$fit
newdata_s$lower_CI <- pmax(pred_s$fit - 1.96*pred_s$se.fit, 0)
newdata_s$upper_CI <- pred_s$fit + 1.96*pred_s$se.fit
newdata_w$predicted_abundance <- pred_w$fit
newdata_w$lower_CI <- pmax(pred_w$fit - 1.96*pred_w$se.fit, 0)
newdata_w$upper_CI <- pred_w$fit + 1.96*pred_w$se.fit

cld_abund_s <- multcomp::cld(emmeans(nb_herb_s, ~ WVclass), adjust = "tukey", Letters = letters)
cld_abund_w <- multcomp::cld(emmeans(nb_herb_w, ~ WVclass), adjust = "tukey", Letters = letters)
cld_abund_s$WVclass <- as.character(cld_abund_s$WVclass)
cld_abund_w$WVclass <- as.character(cld_abund_w$WVclass)
annot_abund_s <- merge(newdata_s, cld_abund_s[, c("WVclass",".group")], by = "WVclass")
annot_abund_w <- merge(newdata_w, cld_abund_w[, c("WVclass",".group")], by = "WVclass")

# -----------------------------------------------------------------------------
# Plots (herbivores; same as your original layout)
# -----------------------------------------------------------------------------
p1 <- ggplot(annot_s, aes(x = WVclass, y = predicted_prob, fill = WVclass)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = .group, y = pmax(predicted_prob - 0.05, 0.02)), size = 4) +
  geom_jitter(data = s.troph_metrics, aes(x = WVclass, y = Herbivore_present),
              inherit.aes = FALSE, width = 0.15, alpha = 0.5, size = 1.2, height = 0.02,
              shape = 21, stroke = 0.6, fill = NA, color = "black") +
  scale_fill_manual(values = veg_colors) +
  labs(x = "Vegetation Class", y = "Predicted Probability") +
  theme_minimal(base_size = 12) +
  ylim(0, 1) +
  theme(legend.position = "none")

p2 <- ggplot(annot_w, aes(x = WVclass, y = predicted_prob, fill = WVclass)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = .group, y = pmax(predicted_prob - 0.05, 0.02)), size = 4) +
  geom_jitter(data = w.troph_metrics, aes(x = WVclass, y = Herbivore_present),
              inherit.aes = FALSE, width = 0.15, alpha = 0.5, size = 1.2, height = 0.02,
              shape = 21, stroke = 0.6, fill = NA, color = "black") +
  scale_fill_manual(values = veg_colors) +
  labs(x = "Vegetation Class", y = "Predicted Probability") +
  theme_minimal(base_size = 12) +
  ylim(0, 1) +
  theme(legend.position = "none")

common_ymax <- max(c(w.troph_metrics$count_2, s.troph_metrics$count_2), na.rm = TRUE) + 1

p3 <- ggplot(annot_abund_s, aes(x = WVclass, y = predicted_abundance, fill = WVclass)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.08) +
  geom_text(aes(label = .group, y = upper_CI + 0.05*common_ymax), size = 4) +
  geom_jitter(data = s.troph_metrics, aes(x = WVclass, y = count_2),
              inherit.aes = FALSE, width = 0.15, alpha = 0.5, size = 1.2,
              shape = 21, stroke = 0.6, fill = NA, color = "black") +
  scale_fill_manual(values = veg_colors) +
  labs(x = "Vegetation Class", y = "Abundance") +
  theme_minimal(base_size = 12) +
  ylim(0, common_ymax) +
  theme(legend.position = "none")

p4 <- ggplot(annot_abund_w, aes(x = WVclass, y = predicted_abundance, fill = WVclass)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.08) +
  geom_text(aes(label = .group, y = upper_CI + 0.05*common_ymax), size = 4) +
  geom_jitter(data = w.troph_metrics, aes(x = WVclass, y = count_2),
              inherit.aes = FALSE, width = 0.15, alpha = 0.5, size = 1.2,
              shape = 21, stroke = 0.6, fill = NA, color = "black") +
  scale_fill_manual(values = veg_colors) +
  labs(x = "Vegetation Class", y = "Abundance") +
  theme_minimal(base_size = 12) +
  ylim(0, common_ymax) +
  theme(legend.position = "none")

final_plot <- (p1 + p2) / (p3 + p4)
ggsave(here("Figures", "Herbivore_Presence_Abundance.png"),
       final_plot, width = 7, height = 6, dpi = 600)

# -----------------------------------------------------------------------------
# Nonparametric examples (leave as examples; extend as needed)
# -----------------------------------------------------------------------------
kr_s <- kruskal.test(PropRich_2 ~ WVclass, data = s.troph_metrics)
pw_s <- pairwise.wilcox.test(s.troph_metrics$PropRich_2, s.troph_metrics$WVclass,
                             p.adjust.method = "holm", exact = FALSE)

kr_w <- kruskal.test(PropRich_2 ~ WVclass, data = w.troph_metrics)
pw_w <- pairwise.wilcox.test(w.troph_metrics$PropRich_2, w.troph_metrics$WVclass,
                             p.adjust.method = "holm", exact = FALSE)

# ============================================================
# OMNIVORES: presence (GLM) + abundance (glmmTMB) 2×2 panel
# ============================================================

# --- presence preds + CLDs (already fit glm_omn_s / glm_omn_w) ---
newdata_s_omn <- data.frame(WVclass = factor(c("Low","Medium","High"), levels = c("Low","Medium","High")),
                            logArea = median_logArea_s)
newdata_w_omn <- data.frame(WVclass = factor(c("Low","Medium","High"), levels = c("Low","Medium","High")),
                            logArea = median_logArea_w)
newdata_s_omn$predicted_prob <- predict(glm_omn_s, newdata = newdata_s_omn, type = "response")
newdata_w_omn$predicted_prob <- predict(glm_omn_w, newdata = newdata_w_omn, type = "response")

cld_omn_s <- cld(emmeans(glm_omn_s, ~ WVclass), adjust = "tukey", Letters = letters)
cld_omn_w <- cld(emmeans(glm_omn_w, ~ WVclass), adjust = "tukey", Letters = letters)
cld_omn_s$WVclass <- as.character(cld_omn_s$WVclass)
cld_omn_w$WVclass <- as.character(cld_omn_w$WVclass)

annot_omn_s <- merge(newdata_s_omn, cld_omn_s[, c("WVclass", ".group")], by = "WVclass")
annot_omn_w <- merge(newdata_w_omn, cld_omn_w[, c("WVclass", ".group")], by = "WVclass")

# --- abundance preds + CLDs (nb_omn_s / nb_omn_w already fit) ---
pred_s_omn <- predict(nb_omn_s, newdata = newdata_s_omn, type = "response", se.fit = TRUE)
pred_w_omn <- predict(nb_omn_w, newdata = newdata_w_omn, type = "response", se.fit = TRUE)

newdata_s_omn$predicted_abundance <- pred_s_omn$fit
newdata_s_omn$lower_CI <- pmax(pred_s_omn$fit - 1.96 * pred_s_omn$se.fit, 0)
newdata_s_omn$upper_CI <- pred_s_omn$fit + 1.96 * pred_s_omn$se.fit
newdata_w_omn$predicted_abundance <- pred_w_omn$fit
newdata_w_omn$lower_CI <- pmax(pred_w_omn$fit - 1.96 * pred_w_omn$se.fit, 0)
newdata_w_omn$upper_CI <- pred_w_omn$fit + 1.96 * pred_w_omn$se.fit

cld_abund_omn_s <- cld(emmeans(nb_omn_s, ~ WVclass), adjust = "tukey", Letters = letters)
cld_abund_omn_w <- cld(emmeans(nb_omn_w, ~ WVclass), adjust = "tukey", Letters = letters)
cld_abund_omn_s$WVclass <- as.character(cld_abund_omn_s$WVclass)
cld_abund_omn_w$WVclass <- as.character(cld_abund_omn_w$WVclass)

annot_abund_omn_s <- merge(newdata_s_omn, cld_abund_omn_s[, c("WVclass", ".group")], by = "WVclass")
annot_abund_omn_w <- merge(newdata_w_omn, cld_abund_omn_w[, c("WVclass", ".group")], by = "WVclass")

# --- y-limits for abundance (use count_3) ---
common_ymax_omn <- max(c(w.troph_metrics$count_3, s.troph_metrics$count_3), na.rm = TRUE) + 1

# --- plots (presence top row) ---
p1_omn <- ggplot(annot_omn_s, aes(x = WVclass, y = predicted_prob, fill = WVclass)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = .group, y = pmax(predicted_prob - 0.05, 0.02)), size = 4) +
  geom_jitter(data = s.troph_metrics, aes(x = WVclass, y = Omnivore_present),
              inherit.aes = FALSE, width = 0.15, alpha = 0.5, size = 1.2, height = 0.02,
              shape = 21, stroke = 0.6, fill = NA, color = "black") +
  scale_fill_manual(values = veg_colors) +
  labs(x = "Vegetation Class", y = "Predicted Probability", title = "Summer") +
  theme_minimal(base_size = 12) + ylim(0, 1) + theme(legend.position = "none")

p2_omn <- ggplot(annot_omn_w, aes(x = WVclass, y = predicted_prob, fill = WVclass)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = .group, y = pmax(predicted_prob - 0.05, 0.02)), size = 4) +
  geom_jitter(data = w.troph_metrics, aes(x = WVclass, y = Omnivore_present),
              inherit.aes = FALSE, width = 0.15, alpha = 0.5, size = 1.2, height = 0.02,
              shape = 21, stroke = 0.6, fill = NA, color = "black") +
  scale_fill_manual(values = veg_colors) +
  labs(x = "Vegetation Class", y = "Predicted Probability", title = "Winter") +
  theme_minimal(base_size = 12) + ylim(0, 1) + theme(legend.position = "none")

# --- plots (abundance bottom row) ---
p3_omn <- ggplot(annot_abund_omn_s, aes(x = WVclass, y = predicted_abundance, fill = WVclass)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.08) +
  geom_text(aes(label = .group, y = upper_CI + 0.05 * common_ymax_omn), size = 4) +
  geom_jitter(data = s.troph_metrics, aes(x = WVclass, y = count_3),
              inherit.aes = FALSE, width = 0.15, alpha = 0.5, size = 1.2,
              shape = 21, stroke = 0.6, fill = NA, color = "black") +
  scale_fill_manual(values = veg_colors) +
  labs(x = "Vegetation Class", y = "Abundance") +
  theme_minimal(base_size = 12) + ylim(0, common_ymax_omn) + theme(legend.position = "none")

p4_omn <- ggplot(annot_abund_omn_w, aes(x = WVclass, y = predicted_abundance, fill = WVclass)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.08) +
  geom_text(aes(label = .group, y = upper_CI + 0.05 * common_ymax_omn), size = 4) +
  geom_jitter(data = w.troph_metrics, aes(x = WVclass, y = count_3),
              inherit.aes = FALSE, width = 0.15, alpha = 0.5, size = 1.2,
              shape = 21, stroke = 0.6, fill = NA, color = "black") +
  scale_fill_manual(values = veg_colors) +
  labs(x = "Vegetation Class", y = "Abundance") +
  theme_minimal(base_size = 12) + ylim(0, common_ymax_omn) + theme(legend.position = "none")

final_plot_omn_combined <- (p1_omn + p2_omn) / (p3_omn + p4_omn)
ggsave("Figures/Omnivore_Presence_Abundance.png", final_plot_omn_combined,
       width = 7, height = 6, dpi = 600)

# ============================================================
# CARNIVORES: presence (GLM) + abundance (glmmTMB) 2×2 panel
# ============================================================

# --- presence preds + CLDs ---
newdata_s_car <- data.frame(WVclass = factor(c("Low","Medium","High"), levels = c("Low","Medium","High")),
                            logArea = median_logArea_s)
newdata_w_car <- data.frame(WVclass = factor(c("Low","Medium","High"), levels = c("Low","Medium","High")),
                            logArea = median_logArea_w)
newdata_s_car$predicted_prob <- predict(glm_carn_s, newdata = newdata_s_car, type = "response")
newdata_w_car$predicted_prob <- predict(glm_carn_w, newdata = newdata_w_car, type = "response")

cld_carn_s <- cld(emmeans(glm_carn_s, ~ WVclass), adjust = "tukey", Letters = letters)
cld_carn_w <- cld(emmeans(glm_carn_w, ~ WVclass), adjust = "tukey", Letters = letters)
cld_carn_s$WVclass <- as.character(cld_carn_s$WVclass)
cld_carn_w$WVclass <- as.character(cld_carn_w$WVclass)

annot_carn_s <- merge(newdata_s_car, cld_carn_s[, c("WVclass", ".group")], by = "WVclass")
annot_carn_w <- merge(newdata_w_car, cld_carn_w[, c("WVclass", ".group")], by = "WVclass")

# --- abundance preds + CLDs ---
pred_s_car <- predict(nb_carn_s, newdata = newdata_s_car, type = "response", se.fit = TRUE)
pred_w_car <- predict(nb_carn_w, newdata = newdata_w_car, type = "response", se.fit = TRUE)

newdata_s_car$predicted_abundance <- pred_s_car$fit
newdata_s_car$lower_CI <- pmax(pred_s_car$fit - 1.96 * pred_s_car$se.fit, 0)
newdata_s_car$upper_CI <- pred_s_car$fit + 1.96 * pred_s_car$se.fit
newdata_w_car$predicted_abundance <- pred_w_car$fit
newdata_w_car$lower_CI <- pmax(pred_w_car$fit - 1.96 * pred_w_car$se.fit, 0)
newdata_w_car$upper_CI <- pred_w_car$fit + 1.96 * pred_w_car$se.fit

cld_abund_carn_s <- cld(emmeans(nb_carn_s, ~ WVclass), adjust = "tukey", Letters = letters)
cld_abund_carn_w <- cld(emmeans(nb_carn_w, ~ WVclass), adjust = "tukey", Letters = letters)
cld_abund_carn_s$WVclass <- as.character(cld_abund_carn_s$WVclass)
cld_abund_carn_w$WVclass <- as.character(cld_abund_carn_w$WVclass)

annot_abund_carn_s <- merge(newdata_s_car, cld_abund_carn_s[, c("WVclass", ".group")], by = "WVclass")
annot_abund_carn_w <- merge(newdata_w_car, cld_abund_carn_w[, c("WVclass", ".group")], by = "WVclass")

# --- y-limits for abundance (use count_4) ---
common_ymax_carn <- max(c(w.troph_metrics$count_4, s.troph_metrics$count_4), na.rm = TRUE) + 1

# --- plots (presence top row) ---
p1_carn <- ggplot(annot_carn_s, aes(x = WVclass, y = predicted_prob, fill = WVclass)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = .group, y = pmax(predicted_prob - 0.05, 0.02)), size = 4) +
  geom_jitter(data = s.troph_metrics, aes(x = WVclass, y = Carnivore_present),
              inherit.aes = FALSE, width = 0.15, alpha = 0.5, size = 1.2, height = 0.02,
              shape = 21, stroke = 0.6, fill = NA, color = "black") +
  scale_fill_manual(values = veg_colors) +
  labs(x = "Vegetation Class", y = "Predicted Probability", title = "Summer") +
  theme_minimal(base_size = 12) + ylim(0, 1) + theme(legend.position = "none")

p2_carn <- ggplot(annot_carn_w, aes(x = WVclass, y = predicted_prob, fill = WVclass)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = .group, y = pmax(predicted_prob - 0.05, 0.02)), size = 4) +
  geom_jitter(data = w.troph_metrics, aes(x = WVclass, y = Carnivore_present),
              inherit.aes = FALSE, width = 0.15, alpha = 0.5, size = 1.2, height = 0.02,
              shape = 21, stroke = 0.6, fill = NA, color = "black") +
  scale_fill_manual(values = veg_colors) +
  labs(x = "Vegetation Class", y = "Predicted Probability", title = "Winter") +
  theme_minimal(base_size = 12) + ylim(0, 1) + theme(legend.position = "none")

# --- plots (abundance bottom row) ---
p3_carn <- ggplot(annot_abund_carn_s, aes(x = WVclass, y = predicted_abundance, fill = WVclass)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.08) +
  geom_text(aes(label = .group, y = upper_CI + 0.05 * common_ymax_carn), size = 4) +
  geom_jitter(data = s.troph_metrics, aes(x = WVclass, y = count_4),
              inherit.aes = FALSE, width = 0.15, alpha = 0.5, size = 1.2,
              shape = 21, stroke = 0.6, fill = NA, color = "black") +
  scale_fill_manual(values = veg_colors) +
  labs(x = "Vegetation Class", y = "Abundance") +
  theme_minimal(base_size = 12) + ylim(0, common_ymax_carn) + theme(legend.position = "none")

p4_carn <- ggplot(annot_abund_carn_w, aes(x = WVclass, y = predicted_abundance, fill = WVclass)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.08) +
  geom_text(aes(label = .group, y = upper_CI + 0.05 * common_ymax_carn), size = 4) +
  geom_jitter(data = w.troph_metrics, aes(x = WVclass, y = count_4),
              inherit.aes = FALSE, width = 0.15, alpha = 0.5, size = 1.2,
              shape = 21, stroke = 0.6, fill = NA, color = "black") +
  scale_fill_manual(values = veg_colors) +
  labs(x = "Vegetation Class", y = "Abundance") +
  theme_minimal(base_size = 12) + ylim(0, common_ymax_carn) + theme(legend.position = "none")

final_plot_carn_combined <- (p1_carn + p2_carn) / (p3_carn + p4_carn)
ggsave("Figures/Carnivore_Presence_Abundance.png", final_plot_carn_combined,
       width = 7, height = 6, dpi = 600)

# -----------------------------------------------------------------------------
# SAVE RESULTS (metrics, models, tests, CLDs, preds, figure)
# -----------------------------------------------------------------------------
OUT_ROOT <- here("results", "trophic")
STAMP    <- format(Sys.Date(), "%Y-%m-%d")
OUT_DIR  <- file.path(OUT_ROOT, STAMP)
fs::dir_create(OUT_DIR, recurse = TRUE)

# metrics
write_csv(s.troph_metrics, file.path(OUT_DIR, "summer_trophic_metrics.csv"))
write_csv(w.troph_metrics, file.path(OUT_DIR, "winter_trophic_metrics.csv"))

# model tables
glm_tables <- bind_rows(
  tidy_glm(glm_herb_s, "summer_herbivore_presence"),
  tidy_glm(glm_omn_s,  "summer_omnivore_presence"),
  tidy_glm(glm_carn_s, "summer_carnivore_presence"),
  tidy_glm(glm_herb_w, "winter_herbivore_presence"),
  tidy_glm(glm_omn_w,  "winter_omnivore_presence"),
  tidy_glm(glm_carn_w, "winter_carnivore_presence")
)
write_csv(glm_tables, file.path(OUT_DIR, "glm_presence_coefficients.csv"))

tmb_tables <- bind_rows(
  tidy_tmb(nb_herb_s, "summer_herbivore_abundance"),
  tidy_tmb(nb_omn_s,  "summer_omnivore_abundance"),
  tidy_tmb(nb_carn_s, "summer_carnivore_abundance"),
  tidy_tmb(nb_herb_w, "winter_herbivore_abundance"),
  tidy_tmb(nb_omn_w,  "winter_omnivore_abundance"),
  tidy_tmb(nb_carn_w, "winter_carnivore_abundance")
)
write_csv(tmb_tables, file.path(OUT_DIR, "glmmTMB_abundance_coefficients.csv"))

# anova (GLM, Type II by default here)
glm_anova <- bind_rows(
  {a<-anova(glm_herb_s, test="Chisq"); tibble::as_tibble(cbind(term=rownames(a), a), .name_repair="minimal") %>% mutate(model="summer_herbivore_presence")},
  {a<-anova(glm_omn_s,  test="Chisq"); tibble::as_tibble(cbind(term=rownames(a), a), .name_repair="minimal") %>% mutate(model="summer_omnivore_presence")},
  {a<-anova(glm_carn_s, test="Chisq"); tibble::as_tibble(cbind(term=rownames(a), a), .name_repair="minimal") %>% mutate(model="summer_carnivore_presence")},
  {a<-anova(glm_herb_w, test="Chisq"); tibble::as_tibble(cbind(term=rownames(a), a), .name_repair="minimal") %>% mutate(model="winter_herbivore_presence")},
  {a<-anova(glm_omn_w,  test="Chisq"); tibble::as_tibble(cbind(term=rownames(a), a), .name_repair="minimal") %>% mutate(model="winter_omnivore_presence")},
  {a<-anova(glm_carn_w, test="Chisq"); tibble::as_tibble(cbind(term=rownames(a), a), .name_repair="minimal") %>% mutate(model="winter_carnivore_presence")}
)
write_csv(glm_anova, file.path(OUT_DIR, "glm_presence_anova.csv"))

# CLDs & prediction newdata (herbivores shown)
write_csv(newdata_s,      file.path(OUT_DIR, "presence_pred_summer_newdata.csv"))
write_csv(newdata_w,      file.path(OUT_DIR, "presence_pred_winter_newdata.csv"))
write_csv(cld_herb_s,     file.path(OUT_DIR, "presence_emmeans_cld_summer.csv"))
write_csv(cld_herb_w,     file.path(OUT_DIR, "presence_emmeans_cld_winter.csv"))
write_csv(annot_abund_s,  file.path(OUT_DIR, "abundance_pred_summer_newdata.csv"))
write_csv(annot_abund_w,  file.path(OUT_DIR, "abundance_pred_winter_newdata.csv"))
write_csv(cld_abund_s,    file.path(OUT_DIR, "abundance_emmeans_cld_summer.csv"))
write_csv(cld_abund_w,    file.path(OUT_DIR, "abundance_emmeans_cld_winter.csv"))

# Nonparametric outputs
sink(file.path(OUT_DIR, "nonparametric_tests.txt"))
cat("== SUMMER ==\n"); print(kr_s); print(pw_s)
cat("\n== WINTER ==\n"); print(kr_w); print(pw_w)
sink()

# Save figure here too
ggsave(file.path(OUT_DIR, "Herbivore_Presence_Abundance.png"),
       final_plot, width = 7, height = 6, dpi = 600)

# Save fitted models bundle
saveRDS(
  list(
    glm_herb_s = glm_herb_s, glm_omn_s = glm_omn_s, glm_carn_s = glm_carn_s,
    glm_herb_w = glm_herb_w, glm_omn_w = glm_omn_w, glm_carn_w = glm_carn_w,
    nb_herb_s = nb_herb_s,   nb_omn_s = nb_omn_s,   nb_carn_s = nb_carn_s,
    nb_herb_w = nb_herb_w,   nb_omn_w = nb_omn_w,   nb_carn_w = nb_carn_w
  ),
  file.path(OUT_DIR, "trophic_models.rds")
)

# omnivore presence/abundance predictions & CLDs
readr::write_csv(newdata_s_omn,      file.path(OUT_DIR, "presence_pred_summer_omnivore.csv"))
readr::write_csv(newdata_w_omn,      file.path(OUT_DIR, "presence_pred_winter_omnivore.csv"))
readr::write_csv(cld_omn_s,          file.path(OUT_DIR, "presence_emmeans_cld_summer_omnivore.csv"))
readr::write_csv(cld_omn_w,          file.path(OUT_DIR, "presence_emmeans_cld_winter_omnivore.csv"))
readr::write_csv(annot_abund_omn_s,  file.path(OUT_DIR, "abundance_pred_summer_omnivore.csv"))
readr::write_csv(annot_abund_omn_w,  file.path(OUT_DIR, "abundance_pred_winter_omnivore.csv"))
readr::write_csv(cld_abund_omn_s,    file.path(OUT_DIR, "abundance_emmeans_cld_summer_omnivore.csv"))
readr::write_csv(cld_abund_omn_w,    file.path(OUT_DIR, "abundance_emmeans_cld_winter_omnivore.csv"))

# carnivore presence/abundance predictions & CLDs
readr::write_csv(newdata_s_car,      file.path(OUT_DIR, "presence_pred_summer_carnivore.csv"))
readr::write_csv(newdata_w_car,      file.path(OUT_DIR, "presence_pred_winter_carnivore.csv"))
readr::write_csv(cld_carn_s,         file.path(OUT_DIR, "presence_emmeans_cld_summer_carnivore.csv"))
readr::write_csv(cld_carn_w,         file.path(OUT_DIR, "presence_emmeans_cld_winter_carnivore.csv"))
readr::write_csv(annot_abund_carn_s, file.path(OUT_DIR, "abundance_pred_summer_carnivore.csv"))
readr::write_csv(annot_abund_carn_w, file.path(OUT_DIR, "abundance_pred_winter_carnivore.csv"))
readr::write_csv(cld_abund_carn_s,   file.path(OUT_DIR, "abundance_emmeans_cld_summer_carnivore.csv"))
readr::write_csv(cld_abund_carn_w,   file.path(OUT_DIR, "abundance_emmeans_cld_winter_carnivore.csv"))

# save the new combined panels
ggsave(file.path(OUT_DIR, "Omnivore_Presence_Abundance.png"),
       final_plot_omn_combined, width = 7, height = 6, dpi = 600)
ggsave(file.path(OUT_DIR, "Carnivore_Presence_Abundance.png"),
       final_plot_carn_combined, width = 7, height = 6, dpi = 600)

message("✅ Trophic analysis complete. Results in: ", OUT_DIR)

