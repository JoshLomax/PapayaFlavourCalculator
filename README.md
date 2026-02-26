# Papaya Flavour Predictor

A Shiny web application for predicting consumer liking scores in papaya fruit based on metabolite concentrations.

## About

This app uses XGBoost machine learning models trained on sensory panel data from 125 untrained consumers across nine distinct papaya genotypes. Users can input metabolite concentrations and receive predicted consumer liking scores.

Two models are available:
- **All Papaya** — trained on red- and yellow-fleshed genotypes
- **Red-Fleshed Only** — trained on red-fleshed genotypes

## Published Article

Lomax, J., Ford, R., Bar, I. (2026). Metabolomic modelling of sensory characteristics and consumer liking in papaya fruit. *Food Chemistry*. https://doi.org/10.1016/j.foodchem.2026.148323

## Repository Contents

```
scripts/Shiny_app_script.R   # Shiny app source code
app_data_files/              # Data files required to run the app
  xgb_combos.RDS             # Trained XGBoost models and performance metrics
  model_data.RDS             # Training data (42 samples)
  all_raw_data.RDS           # Full raw dataset (118 samples)
  papaya_flavour_wheel.png   # Flavour wheel image
  Funding_logo.png           # Hort Innovation logo
  Univiersity_logo.png       # Griffith University logo
```

## Running Locally

```r
shiny::runApp("scripts/Shiny_app_script.R")
```

Required R packages: `shiny`, `bslib`, `DT`, `gt`, `ggplot2`, `openxlsx`, `xgboost`, `dplyr`

## Authors

Joshua Lomax, Rebecca Ford, Ido Bar — Griffith University

Funded by Hort Innovation grant AS19003 (Genetics of fruit sensory preferences).

---

*Note: The Shiny application code was developed with the assistance of [Claude Sonnet 4.6](https://www.anthropic.com/claude) (Anthropic).*
