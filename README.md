# Predictors of advance care planning in Parkinson's disease

Analysis code for a cross-sectional study of people with Parkinson's disease (N = 94). It examines which sociodemographic, clinical and attitudinal factors are associated with having

1. an **advance directive** (AD),
2. a **power of attorney / healthcare proxy** (PoA),
3. **both** documents.

Authors: Anna Pedrosa, David Pedrosa

## Repository contents

| File | Purpose |
|---|---|
| `AnalyseVorsorge_v3_4.R` | Complete analysis: data preparation, descriptive tables, models, figures, export |
| `README.md` | This file |
| `CHANGELOG.md` | Version history |
| `LLM_STATEMENT.md` | Statement on the use of a large language model |

The data file (`Datenmatrix_end-of-life_final_version11.10.23.xlsx`, first sheet `raw_data`) contains patient data and is **not** included. Data may be available from the authors on reasonable request and subject to data protection regulations.

## Requirements

R ≥ 4.1 and the packages `dplyr`, `tidyr`, `ggplot2`, `readxl`, `tableone`, `logistf`, `glmnet`, `pROC`. With `PPR_INSTALL_PACKAGES <- TRUE` (section 0) missing packages are installed automatically.

## How to run

1. Adjust the settings in section 0 if needed, and the working directory in section 2 (the branch for unknown users points to the data folder on the first author's machine).
2. Restart R and source the whole file. Sections depend on objects created earlier, so running them out of order fails.
3. Output is written to `results_final/` inside the working directory.

Runtime is dominated by the resampling steps. Note that both the optimism correction and the stability screen use `min(B_BOOT, 500)` resamples.

## Analysis overview

**Data preparation.** PDQ-39 domain scores and summary index (mean of the eight 0–100 domain scores). Predictors are dichotomised as documented in section 5; "not applicable" answers are kept as their own category (e.g. cohabitation for people in assisted living or nursing homes) rather than treated as missing. Continuous predictors are rescaled to reporting units (age per 5 years, LEDD per 100 mg, UPDRS, PDQ-39 and BDI per 5 points). The Charlson index (without age) is dichotomised (none vs. any), as 69 of 94 participants score 0.

**Primary analysis.** Firth-penalised logistic regression on a fixed set of seven predictors (age, sex, marital status, higher education, stated preferred place of death, LEDD, comorbidity) using all available participants. Odds ratios are reported with profile-likelihood 95% confidence intervals and penalised likelihood-ratio p-values. Discrimination is reported as apparent and bootstrap optimism-corrected AUC (Efron/Harrell).

**Secondary analyses.**
- Unadjusted (univariable) Firth models for all 21 candidate predictors.
- Firth model with all 21 candidate predictors (exploratory; low events per parameter).
- Selection stability of the candidate predictors using lasso-penalised logistic regression on repeated half-sample subsamples, with λ re-tuned by cross-validation in each subsample. Four pure-noise variables are carried through the screen; the highest selection frequency reached by noise ("noise floor") serves as a data-driven reference for interpreting the other predictors.
- Correlation heatmap of the candidate predictors, needed to interpret the stability screen (correlated predictors share selection frequency).
- Sensitivity analysis adding the respective other document as a covariate to the AD and PoA models. It is excluded from the primary models because only 10 of 94 participants hold one document without the other.

## Outputs (`results_final/`)

| File | Content |
|---|---|
| `Table1_raw.csv`, `Table1_by_<outcome>.csv` | Descriptive statistics, overall and stratified by outcome |
| `primaryFirth_<outcome>.csv` | Primary adjusted model |
| `univariable_<outcome>.csv` | Unadjusted associations |
| `fullFirth_<outcome>.csv` | Exploratory 21-predictor model |
| `enetStability_<outcome>.csv` | Selection frequencies, median coefficient when selected, direction |
| `sensitivity_*.csv` | Models including the other document |
| `performance_summary.csv` | N, events, EPP, apparent and optimism-corrected AUC |
| `correlation_heatmap.png`, `correlation_matrix.csv` | Predictor correlations |
| `forest_<outcome>.png`, `stability_<outcome>.png` | Figures |

## Limitations

With N = 94 the number of events per parameter is about 5 in the primary model and well below 2 in the 21-predictor model. The exploratory model and the stability screen are hypothesis-generating; the selection frequencies indicate robustness of selection, not effect size.

## License

The code is released under the MIT License (see `LICENSE`). This applies to the code only, not to the data.
