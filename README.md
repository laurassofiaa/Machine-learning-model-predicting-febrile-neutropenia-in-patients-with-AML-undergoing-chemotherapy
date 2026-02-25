# Predicting Febrile Neutropenia in AML Patients Undergoing Chemotherapy

A machine learning pipeline for binary classification of febrile neutropenia (FN) events during intensive chemotherapy for acute myeloid leukaemia (AML), using XGBoost with SHAP-based feature selection and interpretability analysis.

## Repository Structure

| File | Description |
|------|-------------|
| `AML_infection_model_data.ipynb` | Data preprocessing and preparation for modelling |
| `AML_XGB_classifier_all_variables.ipynb` | XGBoost classifier trained on all features |
| `SHAP_XGBClassifier_plot_function.R` | SHAP summary plots for the all-variable model |
| `SHAP_contributions.R` | Per-feature SHAP contribution plots (all variables) |
| `AML_XGB_classifier_SHAP_top_variables.ipynb` | Refined XGBoost classifier using top SHAP-selected features |
| `TOP_SHAP_XGBClassifier_plot_function.R` | SHAP summary plots for the top-variable model |
| `TOPXGB_shap_contibutions_.R` | Per-feature SHAP contribution plots (top variables) |
| `y_is_zero_topshap_contibutions.R` | SHAP contributions for observations where y = 0 |
| `TOPXGB_SHAP_wrongly_classified_test.R` | SHAP analysis of misclassified observations |
| `TOPXGB_SHAP_wrongly_classified_density_test.R` | Density-based SHAP analysis of misclassified observations |

## Running Order

The pipeline should be run in the following order. Steps grouped together belong to the same analysis stage.

### Stage 1 -- Data preparation

1. **`AML_infection_model_data.ipynb`** -- Preprocess raw data and export model-ready datasets.

### Stage 2 -- Full model (all variables)

2. **`AML_XGB_classifier_all_variables.ipynb`** -- Train and evaluate XGBoost with all features.
3. **`SHAP_XGBClassifier_plot_function.R`** -- Generate SHAP summary plots.
4. **`SHAP_contributions.R`** -- Generate per-feature SHAP contribution plots.

### Stage 3 -- Refined model (top SHAP-selected variables)

5. **`AML_XGB_classifier_SHAP_top_variables.ipynb`** -- Train and evaluate XGBoost with SHAP-selected features.
6. **`TOP_SHAP_XGBClassifier_plot_function.R`** -- Generate SHAP summary plots.
7. **`TOPXGB_shap_contibutions_.R`** -- Generate per-feature SHAP contribution plots.
8. **`y_is_zero_topshap_contibutions.R`** -- Analyse SHAP contributions for y = 0 observations.
9. **`TOPXGB_SHAP_wrongly_classified_test.R`** -- Analyse SHAP values of misclassified observations.
10. **`TOPXGB_SHAP_wrongly_classified_density_test.R`** -- Density-based analysis of misclassified observations.

> **Note:** If SHAP values for y = 0 are not stored, step 5 may fail. In that case, run step 8 (`y_is_zero_topshap_contibutions.R`) first, then re-run steps 5--8.

## Requirements

**Python** (notebooks): `pandas`, `numpy`, `matplotlib`, `seaborn`, `xgboost`, `scikit-learn`, `imblearn`, `scipy`, `shap`

**R** (scripts): `xgboost`, `data.table`, `ggplot2`, `shapviz`, `jsonlite`, `cli`

## Data

All file paths in the code use placeholder paths (`/path/to/...`). Update these to point to your local data directory before running.
