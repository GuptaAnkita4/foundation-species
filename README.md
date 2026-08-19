# Foundation species cover reshapes biodiversity-area relationships and trophic structure in waterbird communities

**R version:** 4.3.1\
**Repository DOI:** (to be added upon acceptance)

------------------------------------------------------------------------

## Overview

This repository contains analysis and figure-generation code for the manuscript:

> *Foundation species cover reshapes biodiversity-area relationships and trophic structure in waterbird communities.*

The project investigates how foundation species (emergent vegetation) and trophic-level structure shape biodiversity–area and abundance–area relationships (SARs and AARs) across small, seasonal wetlands in an irrigated agricultural landscape.

------------------------------------------------------------------------

## Repository structure

    foundation-species/
    ├── data/                    # Input datasets (data/processed/ — the 10 processed CSVs the pipeline reads)
    ├── data/                    # Pipeline-derived outputs (data/derived/<date>/*_trophic_metrics.csv + manifest.json)
    ├── scripts/                 # Analysis scripts (listed below)
    ├── src/                     # Shared helper functions used by every script
    ├── reports/                 # Manuscript text, tables, supplements, and EDA output (reports/eda/)
    ├── Figures/                 # Auto-generated figures
    ├── results/                 # Model outputs and diagnostics, one dated subfolder per script
    ├── run_all.R                # Master pipeline script — runs every step below in order
    ├── foundation-species.Rproj
    ├── README.md
    └── .gitignore, .gitattributes

> All file paths use `{here}` so scripts can be run directly from the repository root, either individually or via `run_all.R`.

### `src/` — shared helper functions

Every script in `scripts/` sources these directly rather than depending on objects left behind by another script (see "Reproduce results" below):

| File | Provides |
|:---|:---|
| `load_data.R` | `load_all_data()` — reads all 10 CSVs from `data/processed/` |
| `build_datasets.R` | `make_sparea()`, `make_trophic_dist()`, `make_glm_dat()`, `make_functional_dat()` — the functions that build every derived dataset used downstream |
| `eda_utils.R` | Plotting/summary helpers used by `01_eda_waterbirds.R` |

------------------------------------------------------------------------

## Software & dependencies

- **R version:** 4.3.1
- **Key packages:**\
  `glmmTMB`, `car`, `emmeans`, `ggeffects`, `DHARMa`, `patchwork`,\
  `dplyr`, `tidyr`, `ggplot2`, `readr`, `fs`, `here`, `purrr`,\
  `GGally`, `sars`, `vegan`, `mgcv`, `multcomp`, `multcompView`,\
  `stringr`, `forcats`, `tibble`, `jsonlite`, `performance`

  `sars` and `DHARMa` in particular require a working compiler toolchain for their own dependencies — install these manually first if a fresh environment doesn't already have them:

  ``` r
  install.packages(c("sars", "DHARMa"))
  ```

### Recommended reproducible setup

``` r
install.packages("renv")
renv::init()
renv::restore()
```

------------------------------------------------------------------------

## Scripts

| Script | Description | Main output |
|:---|:---|:---|
| `00_create_datasets.R` | Builds every derived dataset (species–area tables, trophic distributions, GLM data, functional-diversity tables) from `data/processed/`. | `Data/derived/<date>/{summer,winter}_trophic_metrics.csv` + `manifest.json` |
| `01_eda_waterbirds.R` | Exploratory analyses of richness, abundance, and foundation-species cover distributions. | `reports/eda/*.png` |
| `02_SAR_AAR.R` | Fits baseline Species–Area and Abundance–Area Relationships (`sars::sar_power()`). | Fig. 1, panels d–g (`Figures/SAR_AAR.png`) |
| `02b_figureS1_render.R` | Standalone render of the SAR/AAR fits for the supplement, reading the most recent `results/sar_aar/<date>/` output. | Fig. S1 (`Figures/FigureS1_SAR_AAR.png`) |
| `03_species-area-foundationspecies-models.R` | Tests categorical foundation-species (FS) cover effects on richness and abundance using GLMMs. | Fig. 2 (`Figures/SAR_WV_combined_color.png`), `results/models/<date>/` |
| `03b_analysis_wv_interactions_continuous.R.R` | *(Exploratory — not run by default; see below.)* Continuous/quadratic vegetation-cover interaction models. Not referenced by any current figure or table. | `Figures/SAR_WV_continuous_*.png`, `results/models_continuous/<date>/` |
| `04_area-trophic-interaction.R` | Fits trophic-level SAR and AAR models (herbivores, omnivores, carnivores). | Fig. 4 (`Figures/troph_A_combined_color.png`), `results/trophic_sar/<date>/` |
| `05_fseffects_areaclass.R` | Tests size × FS interactions using discrete wetland size classes. | Fig. 3 (`Figures/hotspots_ci_row.png`), `results/hotspots_area_veg/<date>/` |
| `06_trophic-reponse-to-fs.R` | Models trophic-level (herbivore/omnivore/carnivore) presence and abundance responses to FS cover. | Fig. S3 (`Figures/Herbivore_Presence_Abundance.png` + Omnivore/Carnivore equivalents), `results/trophic/<date>/` |
| `07_biomass_fseffects_supp.R` | *(Exploratory — not run by default; see below.)* Tweedie biomass–area–FS models. Not referenced by any current figure or table. | `Figures/biomass_tweedie_*.png`, `results/biomass_area_wv/<date>/` |
| `08_figureS2_size_veg_diversity.R` | Size × vegetation-class diversity summary. | Fig. S2 (`Figures/FigureS2_size_veg_diversity.png`), `results/figureS2_size_veg/<date>/` |

------------------------------------------------------------------------

## How to reproduce results

### Prepare data

Input data files go under `data/processed/` (10 CSVs: habitat, biodiversity indices, bird occupancy/abundance, species list, and functional-trait data, for summer and winter). Raw data are available upon request.

### Run the full pipeline

The simplest way to reproduce everything is the master script, run from the repository root:

``` bash
Rscript run_all.R                        # core pipeline (skips the two exploratory scripts, 03b and 07)
Rscript run_all.R --include-exploratory   # also runs 03b and 07
```

`run_all.R` runs every step in one shared R session, in manuscript order, and prints a pass/fail summary table at the end. Each script's outputs land in `Figures/` and `results/<name>/<date>/`.

### Run scripts individually

Each can also be run on its own, in any order:

``` bash
Rscript scripts/00_create_datasets.R
Rscript scripts/01_eda_waterbirds.R
Rscript scripts/02_SAR_AAR.R
Rscript scripts/02b_figureS1_render.R
Rscript scripts/03_species-area-foundationspecies-models.R
Rscript scripts/04_area-trophic-interaction.R
Rscript scripts/05_fseffects_areaclass.R
Rscript scripts/06_trophic-reponse-to-fs.R
Rscript scripts/08_figureS2_size_veg_diversity.R

# exploratory, not wired to any manuscript figure/table:
Rscript scripts/03b_analysis_wv_interactions_continuous.R.R
Rscript scripts/07_biomass_fseffects_supp.R
```

------------------------------------------------------------------------

## Outputs

### `Figures/`

Manuscript and supplement figures:

    SAR_AAR.png                        Fig. 1, panels d-g
    SAR_WV_combined_color.png          Fig. 2
    hotspots_ci_row.png                Fig. 3
    troph_A_combined_color.png         Fig. 4
    FigureS1_SAR_AAR.png               Fig. S4.1
    FigureS2_size_veg_diversity.png    Fig. S4.2
    Herbivore_Presence_Abundance.png   Fig. S4.3 (herbivores)

### `results/`

One dated subfolder per script (`results/<name>/<date>/`), each containing:

- Model coefficient tables (`*_coefficients.csv`)
- Type III Wald ANOVA tables (`*_anova.csv` / `*_anova_typeIII.csv`)
- Prediction grids used to build the figures
- Saved model objects (`.rds`)


### A note on numerical reproducibility

`sars::sar_power()` (used by `02_SAR_AAR.R` and `02b_figureS1_render.R`) selects its optimizer's starting values using unseeded random draws internally. Re-running these two scripts can therefore produce SAR/AAR coefficient estimates that differ from a previous run at the 3rd–4th significant digit — this is expected optimizer noise from the `sars` package itself, not a sign of a code or data change, and does not affect any conclusion at the precision reported in the manuscript or figures.

------------------------------------------------------------------------

## Remote sensing and data provenance

- **Imagery sources:** Sentinel-2 SR (20 m), Sentinel-1 GRD (VV/VH), and SRTM DEM
- **Processing environment:** Google Earth Engine
- **Classification:** Random Forest (overall accuracy > 95%, κ > 0.9)
- **Foundation-species (FS) cover:** Defined as emergent/floating vegetation class

All layers were resampled to 20 m resolution and reprojected to EPSG:4326 for consistency.

------------------------------------------------------------------------

## Citation

If using this repository, please cite the following software:

R Core Team (2024). R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/

Brooks ME, Kristensen K, van Benthem KJ, et al. (2017). glmmTMB balances speed and flexibility among packages for zero-inflated generalized linear mixed modeling. The R Journal, 9(2): 378–400.

Fox J, Weisberg S. (2019). An R Companion to Applied Regression. 3rd ed. Sage Publications.

Lüdecke D. (2018). ggeffects: Tidy data frames of marginal effects from regression models. Journal of Open Source Software, 3(26): 772.

Matthews TJ, Triantis KA, Whittaker RJ, Guilhaumon F. (2019). sars: An R package for fitting, evaluating and comparing species–area relationship models. Ecography, 42(8): 1446–1455.

------------------------------------------------------------------------

## License

- **Code:** MIT License
- **Text and figures:** CC BY 4.0

------------------------------------------------------------------------

## Contact

For questions, please contact authors.