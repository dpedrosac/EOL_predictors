# Changelog

## 3.4 – 2026-10-24

### Fixed
- `PPR_INSTALL_PACKAGES` had no effect: the script read an environment variable of the same name instead of the R setting in section 0.
- Stability screen: the median coefficient was taken over all resamples, including those in which the term was not selected (coefficient 0). For any term selected in fewer than half of the resamples the median was therefore 0, and the direction rule labelled it as negative. The median is now taken over the resamples in which the term was selected; terms never selected get direction `0`. This changes `median_coef`, `direction` and the colours of the stability plots, not the selection percentages.

- Stability screen: glmnet convergence warnings at the unpenalised end of the lambda path. The path now stops at 1% of the largest lambda (`STABILITY_LMR`), and the console reports how often the chosen lambda lies at the end of the path.
- Stratified descriptive tables: p-values for continuous variables were all NA because the data were passed to `tableone` as a tibble; now passed as a data frame.
- Forest plot: deprecated `geom_errorbarh()` replaced by `geom_errorbar(orientation = "y")`.

### Documentation
- Sample-specific numbers and interpretive text removed from the script (comments, console verdicts, interpretation guide); the script now only computes and reports.
- Header, rationale and interpretation guide brought in line with the code: script name and version, reporting units (5 rather than 10 points for UPDRS and PDQ), alpha fixed at 1 rather than re-tuned, seven (not five) primary predictors, 21 (not 20) exploratory predictors, rationale for `death_preference_dich` in the primary set.

## 3.0 – 3.3

- Primary analysis changed to Firth-penalised logistic regression on a fixed set of predictors using all participants.
- 80/20 train/test split replaced by bootstrap optimism-corrected AUC.
- ElasticNet no longer reported as coefficients; used as a selection-stability screen only, with λ re-tuned within each resample.
- Noise calibration variables added to the stability screen (3.1).
- Reporting units for continuous predictors; disease severity represented by LEDD instead of UPDRS-III; Charlson index dichotomised (3.1).
- Correlation heatmap of candidate predictors (3.1).
- Random forest removed.
- The other directive excluded from the primary AD and PoA models; available as sensitivity analysis.
- Third outcome is "both documents" (not "either document"); naming made consistent throughout.
- All fits use one design-matrix builder, so factor coding is identical across fits, resamples and predictions.

## 2.x

- Script rebuilt from `AnalyseVorsorge11_2025.R`: three outcomes, each analysed with a full GLM, bidirectional stepwise GLM (AIC) and ElasticNet via `caret`, 80/20 split with repeated 5-fold cross-validation, bootstrap coefficient summaries.
- Coding errors found and corrected against the raw data:
  - sex was coded with the wrong values (1/2 instead of 0/1), labelling all men as women and setting all women to missing;
  - religion code 8 ("no religion") was not mapped and became missing for 35 participants;
  - cohabitation contained no code 0, so the variable was constant; "not applicable" answers of people in institutional care are now their own category;
  - residential location code 4 (n = 11) is treated as rural;
  - the PDQ-39 summary index evaluated to 0 for all participants; domain scores are now computed from explicit item lists.
- `housing_dich` removed from the models: it is identical to the "institutional" cohabitation category.
- Charlson index without age used, to avoid counting age twice.
- 2.2: Firth regression introduced (`AnalyseVorsorge_09_07_26_Firth.R`).

## 1.0

- `AnalyseVorsorge11_2025.R`: descriptive statistics and lasso regression for advance directive and power of attorney.
