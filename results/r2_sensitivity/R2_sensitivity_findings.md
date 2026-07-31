# R² reporting, grid-cell random-effect sensitivity, and leverage check

*Built by running the actual project pipeline (`src/load_data.R`, `src/build_datasets.R`, and the models specified in `scripts/03`, `03b`, `04`, `06`) against the real data in `data/processed/`. Script: `analysis_r2_sensitivity.R` (attached). Full console output and per-model CSVs are in `results/r2_sensitivity/`.*

This directly answers Reviewer 1's General Comment 4 (grid-cell random effect, leverage of small wetlands) and Specific Comment 9 (R² for the GLMs/GLMMs, with vs. without the grid factor).

## Headline finding: the grid-cell random effect means two different things in this analysis, and only one of them matches the Methods justification

Your Methods (§5.5) says: *"We included wetland identity (grid cell) as a random intercept to account for spatial clustering and repeated observations across trophic groups."* That justification is exactly right for **one** family of models and does not apply to **another**:

- **Trophic-level models (script `04`, Fig. 4)** — each wetland contributes **3 rows** (one per trophic group: primary/secondary/top). The grid random intercept is genuinely identified here: it captures real repeated-measures correlation within a wetland across trophic groups, AIC clearly favors keeping it (drops by 8–47 points across the four models), and `performance::r2_nakagawa()` returns sensible, stable marginal/conditional R² (table below).
- **FS-cover models (script `03`, Fig. 2) and the herbivore/omnivore/carnivore presence-abundance models (script `06`, Fig. S3)** — these tables have **exactly one row per wetland** (69 in summer, 70 in winter after NA-filtering — see sample-size note below). A random intercept with one observation per level is structurally non-identifiable: it can't be distinguished from residual (over)dispersion. The tell-tale symptom shows up directly in the fitted objects — the negative-binomial dispersion parameter (θ) explodes from a normal range (1–26) to 900–13,000+ once `(1|grid)` is added, while the random-intercept variance itself stays small. `performance::r2_nakagawa()` and `r2()` both refuse to return a number for these models (they detect the degenerate variance decomposition) — which is why Fig. 2/Fig. S3's models can't get a Nakagawa R² at all; they need a different R² definition (McFadden's, below).

**Practically, this is good news for the paper's conclusions**: the fixed-effect slopes and interaction significance in the FS-cover models barely move whether `(1|grid)` is included or not (see table — e.g., summer richness logArea slope 0.0765 with vs. 0.0764 without). The random effect isn't distorting your results; it just isn't doing the job the Methods sentence describes for these specific models. **Recommended fix**: split that Methods sentence in two — keep "repeated observations across trophic groups" for the Fig. 4 models, and for the Fig. 2 / Fig. S3 models say the random intercept was included as a standard overdispersion control for negative-binomial/logistic count data (or drop it there and report the simpler fixed-effects-only model, since it makes no practical difference — your call).

## R² summary table

**FS-cover models (Fig. 2, script `03`)** — Nakagawa R² is not computable (non-identifiable RE, see above); McFadden's pseudo-R² is reported instead, which works with or without the random term:

| Model | n | logArea slope (w/ grid) | logArea slope (no grid) | McFadden R² (w/ grid) | McFadden R² (no grid) | Interaction p (w/ grid) | Interaction p (no grid) |
|---|---|---|---|---|---|---|---|
| Summer richness | 69 | 0.0765 | 0.0764 | 0.089 | 0.091 | 0.137 | 0.133 |
| Summer abundance | 69 | 0.1013 | 0.0677 | — (see note) | 0.053 | 0.022 | 0.0008 |
| Winter richness | 70 | 0.0783 | 0.0796 | 0.053 | 0.056 | 0.443 | 0.445 |
| Winter abundance | 70 | 0.0642 | 0.0442 | 0.011 | 0.022 | 0.689 | 0.447 |

*Note: the summer-abundance McFadden R² with the grid term hit a numerical edge case in the null-model fit (returned NA) — flagging honestly rather than papering over it; the four other cells all computed cleanly and tell the same story.* The interaction p-value for summer abundance is the one case where inclusion of the random effect visibly changes the p-value (0.022 vs. 0.0008) — both are significant, but you should report the exact one you actually used in the paper and be ready to note this sensitivity if asked.

**Trophic-level models (Fig. 4, script `04`)** — random effect is genuinely identified; Nakagawa's marginal/conditional R² apply cleanly:

| Model | n | AIC (w/ grid) | AIC (no grid) | R²marginal | R²conditional | Interaction p (w/ grid) | Interaction p (no grid) |
|---|---|---|---|---|---|---|---|
| Summer richness | 207 | 793.9 | 802.5 | 0.546 | 0.630 | 0.053 | 0.057 |
| Summer count | 207 | 1589.2 | 1636.4 | 0.344 | 0.688 | 2.3×10⁻⁶ | 3.3×10⁻⁵ |
| Winter richness | 210 | 867.3 | 879.2 | 0.581 | 0.700 | 2.3×10⁻⁴ | 2.1×10⁻³ |
| Winter count | 210 | 1579.2 | 1612.8 | 0.347 | 0.668 | 3.1×10⁻⁴ | 1.1×10⁻² |

Here the random effect earns a large, genuine AIC improvement every time, and the conditional R² (0.63–0.70) is comfortably higher than marginal (0.34–0.58) — exactly what you'd expect when wetland identity legitimately explains extra variance beyond area and trophic level.

**Herbivore/omnivore/carnivore presence (GLM) and abundance (GLMM) models (Fig. S3, script `06`)**: same one-row-per-wetland structure as Fig. 2, so the same non-identifiability applies — confirmed by wildly unstable grid-variance estimates across the six sub-models (from ~10⁻²⁴² to ~5,000) once a random intercept is forced in. McFadden pseudo-R² for the no-grid versions (what the script actually fits) are modest: presence models range 0.0000018–0.39, abundance models 0.014–0.072. These are honestly low — worth stating plainly in the response rather than a specific number if asked, since Reviewer 1 will likely follow up.

*A methodological caveat worth stating in the paper: these McFadden/Nakagawa pseudo-R² values are on a different scale than the R² = 0.13–0.33 reported in §2.1 for the baseline power-law SAR/AAR fits (script `02`, which uses `sars::sar_power`'s own R²). They are not directly comparable — one is an OLS-style R² on log–log-transformed data, the others are GLM/GLMM pseudo-R²s. Don't let a reviewer or reader assume they should match.*

## Leverage of the smallest wetlands — this needs your attention

Reviewer 1 asked directly: *"I also wonder what influence the smallest wetland points have on the slopes... Do they have undue leverage?"* Dropping each of the 5 smallest wetlands one at a time from the FS-cover models and refitting:

- **Summer abundance**: dropping a single wetland — **grid 165** (0.01 ha, the smallest area in the whole summer dataset, low vegetation cover) — changes the area slope from 0.101 to 0.258, a **+154% swing**. That one wetland recorded **112 birds** despite being the smallest surveyed site; the two other 0.01-ha wetlands (grids 294, 216) recorded 18 and 9 birds respectively. This is a strong outlier by any standard.
- **Winter abundance**: dropping **grid 373** (0.05 ha) changes the slope from 0.064 to 0.143, a **+123% swing**. That wetland recorded 69 birds at 0.05 ha.
- **Summer richness**: dropping grid 165 alone shifts the slope by −32%.
- **Winter richness**: dropping grid 373 alone shifts the slope by −19%.
- All other single-point drops (the other 8 of 10 smallest-wetland tests) move slopes by under 1%, which is unremarkable and expected.

**This is worth checking against your field records before responding to the reviewer** — I can't tell from the data alone whether grid 165's 112-bird count / 0.01-ha area combination is a genuine ecological result (e.g., a large mixed flock stopping at a tiny residual pool) or a transcription issue (area or count entered incorrectly, or a count that should reflect a much larger contiguous wetland). Either way, the honest response to Reviewer 1 is: yes, the smallest wetlands do have outsized leverage on the abundance slopes specifically (richness is more robust), concentrated in one wetland per season, and you should report a sensitivity analysis (e.g., slopes with and without these points, or a robust/quasi-Poisson refit) in the supplement rather than leave the single-fit slope unqualified.

## Two sample-size findings that fell out of this exercise

1. **The analytic sample for the FS-cover models is smaller than the survey totals**: 6 of the 75 raw summer wetlands and 1 of the 71 raw winter wetlands have `wetlandTot = 0` (recorded zero area/water/vegetation) and get dropped by the existing NA-filter in `make_glm_dat()`. That leaves **n = 69 (summer) / 70 (winter)** actually feeding Fig. 2/Fig. S3 — not 75/71, and not the 76/71 quoted in the main text and Table S2. Per-class breakdown you can drop straight into the Fig. 2 and Fig. S2 captions Reviewer 1 asked for:
   - Summer WVclass counts (analytic n=69): Low 25, Medium 19, High 25.
   - Winter WVclass counts (analytic n=70): Low 36, Medium 23, High 11.
2. **Some wetlands were visited more than once within a season.** The raw bird-survey tables have more rows than unique wetlands: 9 summer wetlands were visited twice (one, grid 91, three times) and 5 winter wetlands were visited twice (one, grid 2, three times) — all consolidated into one season-level richness/abundance row per wetland in `summer_indices.csv`/`winter_indices.csv` before modeling. This is genuinely useful, verified material for your response to Reviewer 1's General Comment 2 ("was each wetland visited once or repeatedly") — worth stating explicitly, and worth you confirming *why* those particular wetlands got a repeat visit (equipment issue, access, or planned double-observer check?) since I can't tell that from the data alone.

## Files delivered
- `analysis_r2_sensitivity.R` — the full, runnable script (uses your actual `src/` helpers and `data/processed/*.csv`).
- `results/r2_sensitivity/fscover_models_r2.csv`, `trophic_models_r2.csv`, `presence_abundance_models_r2.csv`, `leverage_smallest_wetlands.csv` — every number above, plus the full per-drop leverage table for all 4×5 combinations.
