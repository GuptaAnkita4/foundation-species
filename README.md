# Foundation species amplify biodiversity–area relationships and modify food-web structure in an aquatic ecosystem

**R version:** 4.3.1\
**Repository DOI:** (to be added upon acceptance)

------------------------------------------------------------------------

## Overview

This repository contains analysis and figure-generation code for the
manuscript:

> *Foundation species amplify biodiversity–area relationships and modify
> food-web structure in an aquatic ecosystem.*

The project investigates how foundation species (emergent vegetation)
and trophic-level structure jointly shape biodiversity–area and
abundance–area relationships (SARs and AARs) across small, seasonal
wetlands in an irrigated agricultural landscape.

------------------------------------------------------------------------

## Repository structure

foundation-species/ 
├── data/ \# Input datasets (not tracked if large)
├── scripts/ \# Analysis scripts (listed below) 
├── src/ \# Helper functions and utilities 
├── reports/ \# Manuscript text, figures, tables, supplements 
├── Figures/ \# Auto-generated figures 
├── results/ \# Model outputs and diagnostics 
├── foundation-species.Rproj 
├──README.md 
└── .gitignore, .gitattributes

> All file paths use `{here}` so scripts can be run directly from the
> repository root.

------------------------------------------------------------------------

## Software & dependencies

-   **R version:** 4.3.1\
-   **Key packages:**\
    `glmmTMB`, `car`, `emmeans`, `ggeffects`, `DHARMa`, `patchwork`,\
    `dplyr`, `ggplot2`, `tidyr`, `readr`, `fs`, `here`, `purrr`,\
    `GGally`, `cowplot`, `sars`, `multcompView`, `MASS`

### Recommended reproducible setup

``` r
install.packages("renv")
renv::init()
renv::restore()
```

| Script                                       | Description                                                                                       | Main Output                                 |
|:-----------------|:------------------------------------|:-----------------|
| `00_create_datasets.R`                       | Imports, cleans, and prepares data on species, trophic level, wetland area, and vegetation cover. | Processed `.csv` files in `data/processed/` |
| `01_eda_waterbirds.R`                        | Exploratory analyses of richness, abundance, and foundation-species cover distributions.          | Summary plots and EDA figures               |
| `02_SAR_AARR.R`                              | Fits baseline Species–Area and Abundance–Area Relationships (SARs & AARs).                        | Fig. 1                                      |
| `03_species-area-foundationspecies-models.R` | Tests categorical foundation-species (FS) cover effects on richness and abundance using GLMMs.    | Fig. 2                                      |
| `03b_analysis_wv_interactions_continuous.R`  | Evaluates continuous vegetation–area interactions and generates size-class-specific predictions.  |                                             |
| `04_area-trophic-interaction.R`              | Fits trophic-level SAR and AAR models (herbivores, omnivores, carnivores).                        | Fig. 4                                      |
| `05_fseffects_areaclass.R`                   | Tests size × FS interactions using discrete wetland size classes.                                 | Fig. 3                                      |
| `06_trophic-reponse-to-fs.R`                 | Models trophic-level presence and abundance responses to FS cover.                                | Supplementary Fig.                          |
| `07_biomass_fseffects_supp.R`                | Exploratory analysis of biomass–area–FS relationships.                                            |                                             |


------------------------------------------------------------------------

## How to reproduce results

### Prepare data

Input data files under data/processed.
(Raw data are available upon request.)

### Run analyses

Execute scripts sequentially from the repository root:

Rscript scripts/00_create_datasets.R

Rscript scripts/01_eda_waterbirds.R

Rscript scripts/02_SAR_AARR.R

Rscript scripts/03_species-area-foundationspecies-models.R

Rscript scripts/03b_analysis_wv_interactions_continuous.R

Rscript scripts/04_area-trophic-interaction.R

Rscript scripts/05_fseffects_areaclass.R

Rscript scripts/06_trophic-reponse-to-fs.R

Rscript scripts/07_biomass_fseffects_supp.R

Each script writes outputs to Figures/ and results/<YYYY-MM-DD>/.

------------------------------------------------------------------------

## Outputs

### Figures/

Fig1_SAR_AAR.png

Fig2_WV_interaction.png

Fig3_hotspot_size.png

Fig4_trophic_SAR.png

### results/

Model coefficient tables (*_coefficients.csv)

Type III Wald ANOVA tables (*_anova.csv)

Model objects (.rds)

Prediction grids for figures

### Model details

**Response variables**: Species richness and total abundance

**Predictors**: Log₂(wetland area), vegetation cover, trophic level

**Model type**: Negative binomial GLMMs (glmmTMB, family = nbinom2, link = log)

**Statistical tests**: Type III Wald χ² (car::Anova)

**Pairwise contrasts**: Tukey-adjusted estimated marginal means (emmeans)

**Diagnostics**: Residual uniformity and dispersion tests via DHARMa

### Interpretation:

Exponentiated coefficients represent multiplicative changes in richness or abundance per doubling of wetland area.

------------------------------------------------------------------------

## Remote sensing and data provenance

**Imagery sources**: Sentinel-2 SR (20 m), Sentinel-1 GRD (VV/VH), and SRTM DEM

**Processing environment**: Google Earth Engine

**Classification**: Random Forest (overall accuracy > 95%, κ > 0.9)

**Foundation-species (FS) cover**: Defined as emergent/floating vegetation class

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

**Code**: MIT License

**Text and figures**: CC BY 4.0

**Reproducibility checklist**

 Uses here::here() for file paths

 Includes renv.lock for package reproducibility

 All figures and tables generated automatically

 Random seeds fixed where applicable

 No manual intermediate editing required

------------------------------------------------------------------------

## Contact

For questions or reproducibility issues, please contact authors.
