# ===========================================================
# Project: Parkinson's EOLC Planning Analysis
# Script:  AnalyseVorsorge_ThreeOutcomeModels.R
# Authors: Anna and David Pedrosa
# Date:    2026-06-26
# Version: 2.0 (cleaned; extended to 3 outcomes)
#
# What changed vs. AnalyseVorsorge11_2025.R:
#   - Fixed a syntax error in `codebook` (unescaped quotes around
#     "Rarely" and the religion/worldview entry) that would have
#     stopped the whole script from sourcing.
#   - The labelling loop pointed at an undefined object called `data`,
#     so it never actually ran. I did NOT just redirect it at `eol`,
#     because the dichotomisation step below depends on the raw
#     numeric codes (e.g. gender == 1). Instead, labels are applied
#     to a separate `eol_display` copy, used only for the descriptive
#     table. `eol` itself stays numeric throughout.
#   - Removed a duplicated, redundant single-GLM scratch block.
#   - Removed the "modelSAFETY" results-export block and the broken
#     imputation/FIMA fragment - confirmed both belong to a different,
#     unrelated study/dataset.
#   - Removed the old single-shot LASSO code, per your request.
#   - functionsUsed.R is no longer sourced. results_modelSAFETYmod()
#     and print_AUC() are rebuilt below as self-contained functions
#     (the originals were tied to the unrelated SAFETY study, so
#     reusing them as-is wasn't safe).
#   - The modeling section now runs for THREE outcomes:
#       1) advance_directive_dich   ("AD only")
#       2) power_of_attorney_dich   ("HP only")
#       3) both_directives_dich       ("both existing")
#     each with Full GLM / Stepwise GLM (both directions) / Elastic Net,
#     compared by test-set AUC.
# ===========================================================

# -------------------------------
# 1. Load packages and set directories
# -------------------------------
# (packages.R removed - all needed libraries are loaded directly below instead)

required_pkgs <- c(
  "broom",
  "caret",
  "dplyr",
  "glmnet",
  "MASS",
  "pROC",
  "readxl",
  "tableone",
  "tidyr",
  "tidyverse"
)

install_missing <- isTRUE(as.logical(Sys.getenv("PPR_INSTALL_PACKAGES", "FALSE")))

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  msg <- paste0(
    "Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\n",
    "Install them first, e.g.:\n",
    "  install.packages(c(", paste(sprintf('"%s"', missing_pkgs), collapse = ", "), "))\n",
    "Or set environment variable PPR_INSTALL_PACKAGES=TRUE to auto-install."
  )
  if (!install_missing) {
    stop(msg, call. = FALSE)
  }
  install.packages(missing_pkgs, dependencies = TRUE)
}

suppressPackageStartupMessages({
  for (p in required_pkgs) {
    library(p, character.only = TRUE)
  }
})


#############################################################
## Preamble: setup working directory, load packages and helpers
#############################################################

username <- Sys.info()[["login"]]

if (username == "dpedr") {
  wdir <- "D:/EOL_predictors/"
  data_dir <- file.path(wdir, "data")
} else if (username == "david") {
  wdir <- "/media/storage/EOL_predictors/"
  data_dir <- file.path(wdir, "data")
} else {
  wdir     <- "/Users/annapedrosa/Documents/dokumente/Arbeit_Anna/Aktuelle Projekte/StudiePoD_PD"
  data_dir <- wdir
}

setwd(wdir)
dir.create(file.path(wdir, "results"), showWarnings = FALSE)

# -------------------------------
# 2. Load data and check variables
# -------------------------------
eol <- readxl::read_excel(
  file.path(data_dir, "Datenmatrix_end-of-life_final_version11.10.23.xlsx")
)

# Trims stray whitespace from column names (the real file has " Charlson_withage" with a
# leading space) and fixes a typo in the source spreadsheet (Thoughts_dicussed_with -> ...discussed...)
# so the rest of the script can refer to columns by their intended names.
names(eol) <- trimws(names(eol))
eol <- eol %>% dplyr::rename(Thoughts_discussed_with = Thoughts_dicussed_with)

expected_vars <- c(
  "gender", "age", "age_at_diagnosis", "duration",
  "nationality", "marital_status", "religion_belief_worldview",
  "Housing_situation", "Cohabitation", "Nursing_support",
  "Residential_location", "Professional_education",
  "Advance_directive", "Power_of_attorney",
  "thoughts_about_end_of_life_wishes", "Sharing_of_thoughts",
  "Thoughts_discussed_with",
  "PDQ_1", "PDQ_10", "PDQ_39", # range check markers - NOT the full item list, see note below
  "preferred_place_of_care", "preferred_place_of_death",
  "UPDRS_sum", "bdi_score", "MOCA_score",
  "Charlson_withoutage", "Charlson_withage", "dbs"
)

missing_vars <- setdiff(expected_vars, names(eol))
extra_vars   <- setdiff(names(eol), expected_vars)

if (length(missing_vars) > 0) {
  warning("The following expected variables are missing from the dataset:\n",
          paste(missing_vars, collapse = ", "))
} else {
  message("All expected variables are present.")
}

if (length(extra_vars) > 0) {
  message("The dataset contains extra variables not used in this script:\n",
          paste(extra_vars, collapse = ", "))
}

# -------------------------------
# 3. PDQ subscores
# -------------------------------
# Rewritten to use explicit column name vectors rather than dplyr range-select syntax
# (PDQ_1:PDQ_10 etc.) or matches(), both of which can silently produce zeros inside
# mutate() depending on package load order and R version. Explicit vectors are unambiguous
# and work reliably regardless of column order in the source spreadsheet.

eol <- eol %>%
  mutate(
    PDQ_mobility      = 100/40 * rowSums(across(all_of(paste0("PDQ_", 1:10)))),
    PDQ_ADL           = 100/24 * rowSums(across(all_of(paste0("PDQ_", 11:16)))),
    PDQ_emotional     = 100/24 * rowSums(across(all_of(paste0("PDQ_", 17:22)))),
    PDQ_stigma        = 100/16 * rowSums(across(all_of(paste0("PDQ_", 23:26)))),
    PDQ_support       = 100/12 * rowSums(across(all_of(paste0("PDQ_", 27:29)))),
    PDQ_cognition     = 100/16 * rowSums(across(all_of(paste0("PDQ_", 30:33)))),
    PDQ_communication = 100/12 * rowSums(across(all_of(paste0("PDQ_", 34:36)))),
    PDQ_physical      = 100/12 * rowSums(across(all_of(paste0("PDQ_", 37:39)))),
    PDQ_score         = (PDQ_mobility + PDQ_ADL + PDQ_emotional + PDQ_stigma +
                         PDQ_support + PDQ_cognition + PDQ_communication + PDQ_physical) / 8
  )

# Quick sanity check - should print mean ~25.5, SD ~13.7
message(sprintf("PDQ_score check: mean = %.2f, SD = %.2f (expected ~25.5 / ~13.7)",
                mean(eol$PDQ_score, na.rm = TRUE),
                sd(eol$PDQ_score, na.rm = TRUE)))

# -------------------------------
# 4. Codebook (fixed quoting) - used ONLY for a labelled display copy
# -------------------------------
codebook <- list(
  gender = c("0" = "Female", "1" = "Male", "2" = "Diverse"),
  marital_status = c("0" = "Single", "1" = "Married", "2" = "Registered civil partnership",
                      "3" = "Divorced", "4" = "Dissolved registered partnership", "5" = "Widowed"),
  religion_belief_worldview = c("0" = "Catholic", "1" = "Protestant", "2" = "Other Christian community",
                                 "3" = "Judaism", "4" = "Islam", "5" = "Buddhism", "6" = "Hinduism",
                                 "7" = "Other religion, belief, worldview", "8" = "No religion, belief, worldview"),
  Housing_situation = c("0" = "In own home", "1" = "With children", "2" = "In assisted living",
                         "3" = "Nursing home/senior residence", "4" = "Other"),
  Nursing_support = c("0" = "None", "1" = "Family carer", "2" = "Home care service",
                       "3" = "Personal carer", "4" = "Other"),
  Residential_location = c("0" = "City", "1" = "Middle town", "2" = "Small town",
                            "3" = "Rural town/municipality"),
  School_graduation = c("0" = "None", "1" = "Hauptschule", "2" = "Realschule",
                         "3" = "Fachhochschulreife", "4" = "Abitur", "5" = "Other"),
  Professional_education = c("0" = "None", "1" = "Recognised vocational training",
                              "2" = "Fachhochschule", "3" = "University", "4" = "Other"),
  Advance_directive = c("0" = "No", "1" = "Yes", "2" = "Maybe"),
  Power_of_attorney = c("0" = "No", "1" = "Yes", "2" = "Maybe"),
  palliative_care_knowledge = c("0" = "No", "1" = "Yes"),
  hospice_knowledge = c("0" = "No", "1" = "Yes"),
  thoughts_about_end_of_life_wishes = c("0" = "Never", "1" = "Rarely", "2" = "Now and then", "3" = "Frequently"),
  Sharing_of_thoughts = c("0" = "Never", "1" = "Rarely", "2" = "Now and then", "3" = "Frequently", "4" = "I don't know"),
  Thoughts_discussed_with = c("0" = "Partner", "1" = "Children", "2" = "GP", "3" = "Neurologist",
                               "4" = "Friend", "5" = "Other family member", "6" = "Other"),
  asked_about_end_of_life_wishes = c("0" = "Never", "1" = "Rarely", "2" = "Now and then", "3" = "Frequently", "4" = "I don't know"),
  asked_by_whom = c("0" = "Partner", "1" = "Children", "2" = "GP", "3" = "Neurologist",
                     "4" = "Friend", "5" = "Other family member", "6" = "Other"),
  preferred_place_of_care = c("0" = "Own household", "1" = "Nursing home", "2" = "Hospital", "3" = "Hospice",
                              "4" = "Relative's household", "5" = "Other", "6" = "Not important", "7" = "Palliative care unit"),
  preferred_place_of_death = c("0" = "Own household", "1" = "Nursing home", "2" = "Hospital", "3" = "Hospice",
                               "4" = "Relative's household", "5" = "Other", "6" = "Not important", "7" = "Palliative care unit")
)

# -------------------------------
# 4b. Data cleaning: Residential_location == 4 is not in the documented codebook
# (0=City, 1=Middle town, 2=Small town, 3=Rural town/municipality) and occurred in 11/94
# people. Recoded to 3 (rural) per confirmation - revisit if the raw data gets rechecked.
# Done BEFORE eol_display is created (below) so both the raw descriptive table and the
# modeling data reflect this fix consistently - doing it after left eol_display's 11 people
# as unlabeled NA, since the codebook below only defines labels for 0-3.
# -------------------------------
eol$Residential_location[eol$Residential_location == 4] <- 3

eol_display <- eol  # labelled copy, for the raw/original descriptive table only

for (var in names(codebook)) {
  if (var %in% names(eol_display)) {
    eol_display[[var]] <- factor(eol_display[[var]],
                                  levels = names(codebook[[var]]),
                                  labels = codebook[[var]])
  }
}

# -------------------------------
# 5. Dichotomisation / recoding (unchanged - operates on numeric `eol`)
# -------------------------------
eol <- eol %>%
  mutate(
    # Fixed CRITICAL bug: the old logic (gender==1~Female, gender==2~Male) didn't match
    # either the actual data OR this script's own codebook above (0=Female, 1=Male, 2=Diverse).
    # The real data only contains codes 0 (n=28) and 1 (n=66) - no code 2 occurs. Under the old
    # logic, every gender==1 person was mislabeled "Female", and every gender==0 person matched
    # no condition and became NA, dropping all 28 of them from every single model via
    # complete-case deletion. The variable that survived had only ONE level ("Female") across
    # the whole dataset - this is the direct cause of the "contrasts can be applied only to
    # factors with 2 or more levels" error, and would have failed even without any train/test
    # split. Corrected to match the codebook.
    gender_dich = case_when(
      gender == 0 ~ 0,
      gender == 1 ~ 1,
      gender == 2 ~ NA_real_,  # "Diverse" - not present in this data, but handled explicitly
      TRUE ~ NA_real_
    ),
    gender_dich = factor(gender_dich, levels = c(0, 1), labels = c("Female", "Male")),

    marital_status_dich = case_when(
      marital_status %in% c(1, 2) ~ 1,
      marital_status %in% c(0, 3, 4, 5) ~ 0,
      TRUE ~ NA_real_
    ),
    marital_status_dich = factor(marital_status_dich, levels = c(0, 1),
                                  labels = c("Not married", "Married")),

    nationality_dich = case_when(
      nationality == 0 ~ 0,
      nationality %in% 1:11 ~ 1,
      TRUE ~ NA_real_
    ),
    nationality_dich = factor(nationality_dich, levels = c(0, 1),
                               labels = c("German", "Non-German")),

    # Fixed: the original logic (0:6 ~ "Yes", 7 ~ "No") left code 8 ("No religion, belief,
    # worldview") completely unhandled - it silently became NA for 35/94 people (37%) here.
    # Corrected per discussion: any stated affiliation (0-7, including "Other") vs. no
    # religion at all (8).
    spiritual_orientation_dich = case_when(
      religion_belief_worldview %in% 0:7 ~ 0,
      religion_belief_worldview == 8 ~ 1,
      TRUE ~ NA_real_
    ),
    spiritual_orientation_dich = factor(spiritual_orientation_dich,
                                         levels = c(0, 1),
                                         labels = c("Yes", "No")),

    housing_dich = case_when(
      Housing_situation == 0 ~ 0,
      Housing_situation %in% 1:4 ~ 1,
      TRUE ~ NA_real_
    ),
    housing_dich = factor(housing_dich, levels = c(0, 1), labels = c("Yes", "No")),

    nursing_support_dich = case_when(
      Nursing_support == 0 ~ 0,
      Nursing_support %in% 1:4 ~ 1,
      TRUE ~ NA_real_
    ),
    nursing_support_dich = factor(nursing_support_dich, levels = c(0, 1),
                                   labels = c("No", "Yes")),

    # Fixed: the real data has no Cohabitation == 0 at all (only 1-6 occur), so the original
    # "0 = alone, >0 = with someone" split produced a constant column. Corrected and confirmed:
    # code 1 = alone ("No"), codes 2-6 = living with someone ("Yes").
    # The 4 remaining NAs are NOT random missingness - they're all people in assisted living
    # or a nursing home (Housing_situation 2/3), where "who do you live with" doesn't apply.
    # Coded as an explicit "Institutional" category instead of NA, so complete-case deletion
    # doesn't silently drop these 4 people from every model.
    cohabitation_dich = dplyr::case_when(
      Housing_situation %in% c(2, 3) ~ "Institutional",
      Cohabitation == 1 ~ "No",
      Cohabitation %in% 2:6 ~ "Yes",
      TRUE ~ NA_character_
    ),
    cohabitation_dich = factor(cohabitation_dich, levels = c("No", "Yes", "Institutional")),

    # Residential_location == 4 is recoded to 3 above (step 4b), so the simple split below
    # is sufficient.
    rural_dich = case_when(
      Residential_location == 3 ~ 1,
      Residential_location %in% 0:2 ~ 0,
      TRUE ~ NA_real_
    ),
    rural_dich = factor(rural_dich, levels = c(0, 1), labels = c("No", "Yes")),

    higher_education_dich = case_when(
      Professional_education %in% c(2, 3) ~ 1,
      Professional_education %in% c(0, 1, 4) ~ 0,
      TRUE ~ NA_real_
    ),
    higher_education_dich = factor(higher_education_dich,
                                    levels = c(0, 1),
                                    labels = c("No", "Yes")),

    advance_directive_dich = case_when(
      Advance_directive == 1 ~ 1,
      Advance_directive %in% c(0, 2) ~ 0,
      TRUE ~ NA_real_
    ),
    advance_directive_dich = factor(advance_directive_dich,
                                     levels = c(0, 1),
                                     labels = c("No", "Yes")),

    power_of_attorney_dich = case_when(
      Power_of_attorney == 1 ~ 1,
      Power_of_attorney %in% c(0, 2) ~ 0,
      TRUE ~ NA_real_
    ),
    power_of_attorney_dich = factor(power_of_attorney_dich,
                                     levels = c(0, 1),
                                     labels = c("No", "Yes")),

    end_of_life_thoughts_dich = case_when(
      thoughts_about_end_of_life_wishes == 0 ~ 0,
      thoughts_about_end_of_life_wishes %in% 1:3 ~ 1,
      TRUE ~ NA_real_
    ),
    end_of_life_thoughts_dich = factor(end_of_life_thoughts_dich,
                                        levels = c(0, 1),
                                        labels = c("No", "Yes")),

    sharing_thoughts_dich = case_when(
      Sharing_of_thoughts %in% c(1, 2, 3) ~ 1,
      Sharing_of_thoughts %in% c(0, 4) ~ 0,
      TRUE ~ NA_real_
    ),
    sharing_thoughts_dich = factor(sharing_thoughts_dich,
                                    levels = c(0, 1),
                                    labels = c("No", "Yes")),

    thoughts_discussed_with_cat = case_when(
      Thoughts_discussed_with %in% c(0, 1, 4, 5) ~ "Family/Friends",
      Thoughts_discussed_with %in% c(2, 3)       ~ "Physician",
      Thoughts_discussed_with == 6               ~ "Other",
      TRUE ~ NA_character_
    ),
    thoughts_discussed_with_cat = factor(thoughts_discussed_with_cat,
                                          levels = c("Family/Friends", "Physician", "Other"))
  )

eol <- eol %>%
  mutate(
    dbs_dich = case_when(
      dbs == 0 ~ 0,
      dbs == 1 ~ 1,
      TRUE ~ NA_real_
    ),
    dbs_dich = factor(dbs_dich, levels = c(0, 1), labels = c("No", "Yes"))
  )

# Both directives defined here (after advance_directive_dich and power_of_attorney_dich
# exist) so it's available for the descriptive table in section 7 below.
eol <- eol %>%
  mutate(
    both_directives_dich = case_when(
      advance_directive_dich == "Yes" & power_of_attorney_dich == "Yes" ~ 1,
      advance_directive_dich == "No"  | power_of_attorney_dich == "No"  ~ 0,
      TRUE ~ NA_real_
    ),
    both_directives_dich = factor(both_directives_dich,
                                   levels = c(0, 1),
                                   labels = c("No", "Yes"))
  )

# --- Existence of a stated preference for place of care / place of death ---
# Coded as: Yes = any specific location chosen (codes 0-5, 7); No = "Not important" (code 6).
# True missing responses stay NA. This captures WHETHER a preference exists, not WHICH
# location was chosen.
eol <- eol %>%
  mutate(
    care_preference_dich = case_when(
      preferred_place_of_care %in% c(0, 1, 2, 3, 4, 5, 7) ~ 1,
      preferred_place_of_care == 6 ~ 0,
      TRUE ~ NA_real_
    ),
    care_preference_dich = factor(care_preference_dich, levels = c(0, 1), labels = c("No", "Yes")),

    death_preference_dich = case_when(
      preferred_place_of_death %in% c(0, 1, 2, 3, 4, 5, 7) ~ 1,
      preferred_place_of_death == 6 ~ 0,
      TRUE ~ NA_real_
    ),
    death_preference_dich = factor(death_preference_dich, levels = c(0, 1), labels = c("No", "Yes"))
  )

# -------------------------------
# 6. Descriptive table (raw, labelled) - uses eol_display
# -------------------------------
allVars_original <- c(
  "gender", "age", "age_at_diagnosis", "duration",
  "nationality", "marital_status", "religion_belief_worldview",
  "Housing_situation", "Cohabitation", "Nursing_support",
  "Residential_location", "Professional_education",
  "Advance_directive", "Power_of_attorney",
  "thoughts_about_end_of_life_wishes", "Sharing_of_thoughts",
  "Thoughts_discussed_with", "PDQ_score", "UPDRS_sum", "bdi_score",
  "MOCA_score", "Charlson_withoutage", "Charlson_withage", "dbs",
  "preferred_place_of_care", "preferred_place_of_death"
)

catVars_original <- c(
  "gender", "nationality", "marital_status",
  "religion_belief_worldview", "Housing_situation",
  "Cohabitation", "Nursing_support", "Residential_location",
  "Professional_education", "Advance_directive", "Power_of_attorney",
  "thoughts_about_end_of_life_wishes", "Sharing_of_thoughts",
  "Thoughts_discussed_with", "dbs", "preferred_place_of_care", "preferred_place_of_death"
)

table_descriptive <- CreateTableOne(
  vars       = allVars_original,
  data       = eol_display,
  factorVars = catVars_original
)
print(table_descriptive)

# -------------------------------
# 7. Analytic table (dichotomised) - uses eol
# -------------------------------
allVars_dich <- c(
  "gender_dich", "age", "age_at_diagnosis", "duration",
  "nationality_dich", "marital_status_dich", "spiritual_orientation_dich",
  "housing_dich", "cohabitation_dich", "nursing_support_dich", "rural_dich",
  "higher_education_dich", "advance_directive_dich", "power_of_attorney_dich",
  "both_directives_dich",
  "end_of_life_thoughts_dich", "sharing_thoughts_dich",
  "thoughts_discussed_with_cat", "PDQ_score", "UPDRS_sum", "bdi_score",
  "MOCA_score", "Charlson_withoutage", "Charlson_withage", "dbs",
  "care_preference_dich", "death_preference_dich"
)

catVars_dich <- c(
  "gender_dich", "nationality_dich", "marital_status_dich",
  "spiritual_orientation_dich", "housing_dich", "cohabitation_dich",
  "nursing_support_dich", "rural_dich", "higher_education_dich",
  "advance_directive_dich", "power_of_attorney_dich", "both_directives_dich",
  "end_of_life_thoughts_dich", "sharing_thoughts_dich",
  "thoughts_discussed_with_cat", "dbs", "care_preference_dich", "death_preference_dich"
)

table_dich <- CreateTableOne(vars = allVars_dich, data = eol, factorVars = catVars_dich)
print(table_dich)

# ================================================================================================
# 9. Modeling helpers: Full GLM / Stepwise GLM (both directions) / Elastic Net
# ================================================================================================
# These replace functionsUsed.R's results_modelSAFETYmod()/print_AUC(), which were built for an
# unrelated study and dataset. Everything the pipeline needs now lives in this script.

fit_three_models <- function(outcome, predictors, data, train_control,
                              lambda_grid = seq(0.0001, 1, length.out = 100),
                              alpha_grid  = seq(0, 1, length.out = 11),
                              p_split = 0.8, seed = 123) {

  set.seed(seed)

  model_data <- data %>%
    dplyr::select(dplyr::all_of(c(outcome, predictors))) %>%
    tidyr::drop_na(dplyr::all_of(outcome))   # predictor-level NAs are handled per-model via na.action

  form <- as.formula(paste(outcome, "~", paste(predictors, collapse = " + ")))

  index      <- caret::createDataPartition(model_data[[outcome]], p = p_split, list = FALSE)
  train_data <- model_data[index, ]
  test_data  <- model_data[-index, ]

  fit_full <- caret::train(
    form, data = train_data, method = "glm", family = binomial,
    trControl = train_control, metric = "logLoss", na.action = na.omit
  )

  fit_step <- caret::train(
    form, data = train_data, method = "glmStepAIC", family = binomial,
    trControl = train_control, metric = "logLoss", na.action = na.omit,
    direction = "both", trace = FALSE
  )

  grid_en <- expand.grid(alpha = alpha_grid, lambda = lambda_grid)
  fit_pen <- caret::train(
    form, data = train_data, method = "glmnet", family = "binomial",
    trControl = train_control, tuneGrid = grid_en, metric = "logLoss", na.action = na.omit
  )

  list(full = fit_full, step = fit_step, pen = fit_pen,
       model_data = model_data, train_data = train_data, test_data = test_data, outcome = outcome)
}

# Accuracy / AUC / logLoss for one fitted model, evaluated on its held-out test set.
# This matches the published methods, which compared models on all three metrics, not AUC alone.
# Accuracy / AUC / logLoss for one fitted model, evaluated on its held-out test set.
# This matches the published methods, which compared models on all three metrics, not AUC alone.
evaluate_test_metrics <- function(fit, test_data, outcome) {
  # predict() silently DROPS rows with NA in any predictor used by the fitted model (because
  # na.action = na.omit was used when fitting), but `truth` was being read from the full,
  # un-filtered test_data - causing probs/preds to come back shorter than truth and breaking
  # pROC::roc() with a length-mismatch error. Restricting test_data to complete cases on the
  # model's predictors up front keeps probs/preds/truth aligned to the same rows.
  model_vars     <- all.vars(formula(fit))[-1]  # drops the outcome (first term) from the formula
  test_data_cc   <- test_data[stats::complete.cases(test_data[, model_vars]), ]

  probs <- predict(fit, newdata = test_data_cc, type = "prob")[, "Yes"]
  preds <- predict(fit, newdata = test_data_cc, type = "raw")
  truth <- test_data_cc[[outcome]]

  roc_obj <- pROC::roc(response = truth, predictor = probs,
                        levels = c("No", "Yes"), direction = "<", quiet = TRUE)
  auc_val <- as.numeric(pROC::auc(roc_obj))
  acc_val <- mean(preds == truth, na.rm = TRUE)

  eps         <- 1e-15
  p_clamped   <- pmin(pmax(probs, eps), 1 - eps)
  y_binary    <- as.numeric(truth == "Yes")
  logloss_val <- -mean(y_binary * log(p_clamped) + (1 - y_binary) * log(1 - p_clamped))

  list(roc = roc_obj, AUC = auc_val, Accuracy = acc_val, logLoss = logloss_val)
}

# Combined ROC plot for Full / Stepwise / Elastic Net, with AUC + 95% CI in the legend
plot_auc_comparison <- function(roc_list, title) {
  legend_labels <- vapply(names(roc_list), function(nm) {
    auc_val <- as.numeric(pROC::auc(roc_list[[nm]]))
    # ci.auc() can fail ("need at least two non-NA values to interpolate") when a model hits
    # perfect/near-perfect separation - all predicted probabilities collapse toward 0/1, giving
    # a degenerate ROC curve with no interpolation room. With ~20 predictors and only ~75
    # training rows, this is a plausible outcome for the Full (and sometimes Stepwise) model,
    # not a bug - it's itself a sign that model is overfit/unstable. Caught here so one
    # degenerate model doesn't crash the whole pipeline; the AUC point estimate is still shown.
    ci <- tryCatch(as.numeric(pROC::ci.auc(roc_list[[nm]])), error = function(e) {
      warning(sprintf(
        "ci.auc() failed for model '%s' (%s) - likely perfect/near-perfect separation given ~20 predictors vs N=94. AUC point estimate shown without a CI.",
        nm, title
      ))
      c(NA_real_, auc_val, NA_real_)
    })
    if (anyNA(ci)) {
      sprintf("%s (AUC = %.2f [CI unavailable])", nm, auc_val)
    } else {
      sprintf("%s (AUC = %.2f [%.2f-%.2f])", nm, ci[2], ci[1], ci[3])
    }
  }, character(1))
  names(roc_list) <- legend_labels

  pROC::ggroc(roc_list, legacy.axes = TRUE) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey70") +
    ggplot2::labs(title = title, x = "1 - Specificity", y = "Sensitivity", colour = "Model") +
    ggplot2::theme_minimal()
}

# All terms, p-values, and parametric (profile-likelihood) 95% CIs for a glm-type fit (Full / Stepwise)
extract_glm_coefs <- function(fit, model_label, alpha = 0.05) {
  # Profile-likelihood CIs (the default for glm via MASS::confint) can also fail to converge
  # under separation - same root cause as the ci.auc() issue above. Falls back to coefficients
  # and p-values without CIs rather than crashing the whole pipeline.
  tidy_out <- tryCatch(
    broom::tidy(fit$finalModel, conf.int = TRUE, conf.level = 1 - alpha),
    error = function(e) {
      warning(sprintf(
        "Profile-likelihood CIs failed for model '%s' - likely separation given ~20 predictors vs N=94. Coefficients/p-values shown without CIs.",
        model_label
      ))
      broom::tidy(fit$finalModel, conf.int = FALSE)
    }
  )
  tidy_out %>%
    dplyr::mutate(model = model_label, significant = p.value < alpha)
}

# Non-zero terms for the Elastic Net fit at its cross-validated best lambda
extract_glmnet_coefs <- function(fit, model_label) {
  coefs <- as.matrix(coef(fit$finalModel, s = fit$bestTune$lambda))
  data.frame(term = rownames(coefs), estimate = coefs[, 1], model = model_label) %>%
    dplyr::filter(estimate != 0)
}

# --------------------------------------------------------------------------------------------------
# Bootstrap coefficient stability / 95% CIs (default 1000 resamples), applied to the Stepwise and
# Elastic Net models only - matches the published methods (the Full model's CIs come from the
# parametric profile-likelihood interval in extract_glm_coefs() above instead). Resampling is done
# WITH replacement from the FULL dataset (not just the 80% training partition) - confirmed by the
# results_bootstrap() call in the SAFETY study, which passes the full pre-split data. The 80/20
# split is reserved for the Accuracy/AUC/logLoss performance comparison; this bootstrap is a
# separate analysis aimed at coefficient stability/uncertainty using all available data. Each
# resample is refit with the identical procedure: AIC-based bidirectional stepwise selection, or
# glmnet at the alpha/lambda already chosen via cross-validation on the real (non-bootstrapped) data.
# --------------------------------------------------------------------------------------------------

bootstrap_stepwise_coefs <- function(source_data, outcome, predictors, B = 1000, seed = 123) {
  set.seed(seed)
  form <- as.formula(paste(outcome, "~", paste(predictors, collapse = " + ")))
  n <- nrow(source_data)

  boot_list   <- vector("list", B)
  error_msgs  <- character(0)
  for (b in seq_len(B)) {
    boot_data <- source_data[sample.int(n, n, replace = TRUE), ]
    fit_b <- tryCatch({
      full_b <- glm(form, data = boot_data, family = binomial)
      MASS::stepAIC(full_b, direction = "both", trace = FALSE)
    }, error = function(e) {
      error_msgs[[length(error_msgs) + 1]] <<- conditionMessage(e)
      NULL
    })
    if (!is.null(fit_b)) {
      cf <- coef(fit_b)
      boot_list[[b]] <- data.frame(term = names(cf), estimate = as.numeric(cf))
    }
  }

  n_success <- sum(!vapply(boot_list, is.null, logical(1)))
  message(sprintf("bootstrap_stepwise_coefs: %d/%d resamples succeeded.", n_success, B))
  if (n_success == 0) {
    warning("bootstrap_stepwise_coefs: ALL resamples failed. First distinct error message(s):\n  - ",
            paste(unique(error_msgs)[seq_len(min(3, length(unique(error_msgs))))], collapse = "\n  - "))
    return(data.frame(term = character(0), mean = numeric(0), median = numeric(0),
                       sd = numeric(0), ci_lower = numeric(0), ci_upper = numeric(0),
                       selected_pct = numeric(0)))
  }

  dplyr::bind_rows(boot_list) %>%
    dplyr::group_by(term) %>%
    dplyr::summarise(
      mean         = mean(estimate),
      median       = stats::median(estimate),
      sd           = stats::sd(estimate),
      ci_lower     = stats::quantile(estimate, .025),
      ci_upper     = stats::quantile(estimate, .975),
      selected_pct = dplyr::n() / B * 100,   # % of all B resamples in which this term survived stepAIC
      .groups = "drop"
    )
}

bootstrap_glmnet_coefs <- function(source_data, outcome, predictors, alpha, lambda, B = 1000, seed = 123) {
  set.seed(seed)
  form <- as.formula(paste(outcome, "~", paste(predictors, collapse = " + ")))
  n <- nrow(source_data)

  boot_list  <- vector("list", B)
  error_msgs <- character(0)
  for (b in seq_len(B)) {
    boot_data <- source_data[sample.int(n, n, replace = TRUE), ]
    fit_b <- tryCatch({
      X_b <- model.matrix(form, data = boot_data)[, -1, drop = FALSE]  # glmnet adds its own intercept
      y_b <- boot_data[[outcome]]
      glmnet::glmnet(X_b, y_b, family = "binomial", alpha = alpha, lambda = lambda)
    }, error = function(e) {
      error_msgs[[length(error_msgs) + 1]] <<- conditionMessage(e)
      NULL
    })
    if (!is.null(fit_b)) {
      cf <- as.matrix(coef(fit_b, s = lambda))
      boot_list[[b]] <- data.frame(term = rownames(cf), estimate = cf[, 1])
    }
  }

  n_success <- sum(!vapply(boot_list, is.null, logical(1)))
  message(sprintf("bootstrap_glmnet_coefs: %d/%d resamples succeeded.", n_success, B))
  if (n_success == 0) {
    warning("bootstrap_glmnet_coefs: ALL resamples failed. First distinct error message(s):\n  - ",
            paste(unique(error_msgs)[seq_len(min(3, length(unique(error_msgs))))], collapse = "\n  - "))
    return(data.frame(term = character(0), mean = numeric(0), median = numeric(0),
                       sd = numeric(0), ci_lower = numeric(0), ci_upper = numeric(0),
                       nonzero_pct = numeric(0)))
  }

  dplyr::bind_rows(boot_list) %>%
    dplyr::group_by(term) %>%
    dplyr::summarise(
      mean        = mean(estimate),
      median      = stats::median(estimate),
      sd          = stats::sd(estimate),
      ci_lower    = stats::quantile(estimate, .025),
      ci_upper    = stats::quantile(estimate, .975),
      nonzero_pct = mean(estimate != 0) * 100,   # % of all B resamples in which this term stayed non-zero
      .groups = "drop"
    )
}

# Runs all three models for one outcome: fits, evaluates (Accuracy/AUC/logLoss), plots ROC,
# extracts coefficients (with parametric CIs for Full/Stepwise), and bootstraps CIs for
# Stepwise/ElasticNet coefficients.
run_pipeline <- function(outcome, predictors, data, train_control, label, B = 1000) {

  fits <- fit_three_models(outcome, predictors, data, train_control)

  m_full <- evaluate_test_metrics(fits$full, fits$test_data, outcome)
  m_step <- evaluate_test_metrics(fits$step, fits$test_data, outcome)
  m_pen  <- evaluate_test_metrics(fits$pen,  fits$test_data, outcome)

  metrics <- data.frame(
    model    = c("Full", "Stepwise", "ElasticNet"),
    Accuracy = c(m_full$Accuracy, m_step$Accuracy, m_pen$Accuracy),
    AUC      = c(m_full$AUC,      m_step$AUC,      m_pen$AUC),
    logLoss  = c(m_full$logLoss,  m_step$logLoss,  m_pen$logLoss)
  )

  roc_list <- list(Full = m_full$roc, Stepwise = m_step$roc, ElasticNet = m_pen$roc)
  fig <- plot_auc_comparison(roc_list, title = paste("ROC comparison -", label))

  coefs <- dplyr::bind_rows(
    extract_glm_coefs(fits$full, "Full"),
    extract_glm_coefs(fits$step, "Stepwise"),
    extract_glmnet_coefs(fits$pen, "ElasticNet")
  )

  # Bootstrap from the FULL dataset (model_data), not just the training partition - see note
  # above the bootstrap functions. Cached to disk since 1000 resamples x 2 methods is slow;
  # delete the .Rdata file in results/ to force a recompute (e.g. after changing predictors).
  cache_file <- file.path(wdir, "results", paste0("CIdataBootstrap_", outcome, ".Rdata"))
  if (file.exists(cache_file)) {
    load(cache_file)   # restores boot_step, boot_pen
  } else {
    boot_step <- bootstrap_stepwise_coefs(fits$model_data, outcome, predictors, B = B)
    boot_pen  <- bootstrap_glmnet_coefs(
      fits$model_data, outcome, predictors,
      alpha  = fits$pen$bestTune$alpha,
      lambda = fits$pen$bestTune$lambda,
      B = B
    )
    save(boot_step, boot_pen, file = cache_file)
  }
  boot_ci <- dplyr::bind_rows(
    dplyr::mutate(boot_step, model = "Stepwise"),
    dplyr::mutate(boot_pen,  model = "ElasticNet")
  )

  list(fits = fits, metrics = metrics, roc = roc_list, plot = fig, coefs = coefs, boot_ci = boot_ci)
}

# ================================================================================================
# 10. Predictor sets per outcome
# ================================================================================================
# both_directives_dich is DEFINED from advance_directive_dich & power_of_attorney_dich (step 8),
# so neither can be used as a predictor for that outcome - that would leak the outcome into the
# predictors and produce a meaningless, near-perfect "model". For the AD-only and HP-only models,
# the *other* directive variable is kept as a covariate, since the two are correlated but not
# identical (mirrors the approach in your earlier study).

# Using Charlson_withoutage (not Charlson_withage) since the latter already bakes age
# into the comorbidity score - using both alongside a standalone `age` term would have
# double-counted age's contribution.
predictors_common <- c(
  "age", "gender_dich", "nationality_dich", "marital_status_dich",
  "spiritual_orientation_dich", "housing_dich", "cohabitation_dich",
  "nursing_support_dich", "rural_dich", "higher_education_dich",
  "end_of_life_thoughts_dich", "sharing_thoughts_dich",
  "care_preference_dich", "death_preference_dich",
  "duration", "dbs_dich", "UPDRS_sum", "PDQ_score", "MOCA_score",
  "bdi_score", "Charlson_withoutage"
)

predictors_AD  <- c(predictors_common, "power_of_attorney_dich")
predictors_HP  <- c(predictors_common, "advance_directive_dich")
predictors_ANY <- predictors_common

# ================================================================================================
# 11. trainControl (shared across all 3 outcomes x 3 models)
# ================================================================================================
# summaryFunction = mnLogLoss (not twoClassSummary/ROC) - matches the trainControl used in your
# other studies on this dataset family, where log-loss drove the CV-based hyperparameter selection.
train_control <- caret::trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 10,
  classProbs = TRUE,
  summaryFunction = caret::mnLogLoss,
  savePredictions = "final",
  verboseIter = TRUE
)

# ================================================================================================
# 12. Run the three outcomes
# ================================================================================================
# SAMPLE SIZE CAVEAT: With N=94 (~75 in the training partition after the 80/20 split) and
# ~20-21 predictors per model, events-per-predictor is well below standard guidance (Peduzzi
# et al. recommend >=10, ideally >=20). Deliberate decision: keep the full predictor set rather
# than trim it based on a priori judgment, and treat ElasticNet (built for exactly this
# high-predictor/limited-N regime) as the PRIMARY result. Full GLM and Stepwise GLM are reported
# as EXPLORATORY only - expect convergence warnings, separation (perfect/near-perfect AUC with
# no usable CI - see the tryCatch warnings above), and unstable coefficient selection. The
# bootstrap `selected_pct` (Stepwise) and `nonzero_pct` (ElasticNet) columns in boot_ci, below,
# quantify exactly how unstable that selection is across 1000 resamples - useful evidence for
# this caveat when reporting results.
message(paste(strwrap(
  "Sample size caveat: N=94, ~20-21 predictors per model. Full/Stepwise GLM results are
   exploratory due to low events-per-predictor; ElasticNet is the primary result.",
  width = 100), collapse = "\n"))

# B_BOOT controls the number of bootstrap resamples used for the Stepwise/ElasticNet coefficient
# CIs. The published methods used 1000. Each run_pipeline() call fits roughly 2 * B_BOOT extra
# models on top of the CV-tuned Full/Stepwise/ElasticNet fits, so this is by far the slowest part
# of the script. Lower it (e.g. to 100) for a quick test run, then set back to 1000 for the
# final analysis.
B_BOOT <- 1000

res_AD  <- run_pipeline("advance_directive_dich", predictors_AD,  eol, train_control, "Advance Directive", B = B_BOOT)
res_HP  <- run_pipeline("power_of_attorney_dich", predictors_HP,  eol, train_control, "Healthcare Proxy", B = B_BOOT)
res_ANY <- run_pipeline("both_directives_dich",      predictors_ANY, eol, train_control, "Both Directives (AD and HP)", B = B_BOOT)

res_AD$plot
res_HP$plot
res_ANY$plot

# ================================================================================================
# 13. Console + file summary of Accuracy / AUC / logLoss across all 9 fitted models
# ================================================================================================
metric_summary <- dplyr::bind_rows(
  dplyr::mutate(res_AD$metrics,  outcome = "Advance Directive"),
  dplyr::mutate(res_HP$metrics,  outcome = "Healthcare Proxy"),
  dplyr::mutate(res_ANY$metrics, outcome = "Both Directives")
)
print(metric_summary)
write.csv(metric_summary, file.path(wdir, "results", "TestSetMetrics_allModels.csv"), row.names = FALSE)

# ================================================================================================
# 14. Export coefficient tables (with parametric CIs for Full/Stepwise) + ROC plots
# ================================================================================================
write.csv(res_AD$coefs,  file.path(wdir, "results", "coefficients_AdvanceDirective.csv"), row.names = FALSE)
write.csv(res_HP$coefs,  file.path(wdir, "results", "coefficients_HealthcareProxy.csv"),  row.names = FALSE)
write.csv(res_ANY$coefs, file.path(wdir, "results", "coefficients_AnyDirective.csv"),     row.names = FALSE)

ggplot2::ggsave(file.path(wdir, "results", "ROC_AdvanceDirective.png"), res_AD$plot,  width = 6, height = 5)
ggplot2::ggsave(file.path(wdir, "results", "ROC_HealthcareProxy.png"), res_HP$plot,  width = 6, height = 5)
ggplot2::ggsave(file.path(wdir, "results", "ROC_AnyDirective.png"),    res_ANY$plot, width = 6, height = 5)

# ================================================================================================
# 15. Export bootstrap coefficient stability / 95% CIs (Stepwise & ElasticNet only)
# ================================================================================================
write.csv(res_AD$boot_ci,  file.path(wdir, "results", "bootstrapCI_AdvanceDirective.csv"), row.names = FALSE)
write.csv(res_HP$boot_ci,  file.path(wdir, "results", "bootstrapCI_HealthcareProxy.csv"),  row.names = FALSE)
write.csv(res_ANY$boot_ci, file.path(wdir, "results", "bootstrapCI_AnyDirective.csv"),     row.names = FALSE)

# ================================================================================================
# 16. Print readable results summary to console for each outcome
# ================================================================================================
# Shows: test-set performance, significant Full/Stepwise predictors, ElasticNet non-zero
# coefficients, and bootstrap stability - all in one place per outcome.
# Remember: Full/Stepwise results are EXPLORATORY given N=94 vs ~21 predictors (see section 12
# caveat). ElasticNet + bootstrap stability are the primary, trustworthy results.

summarize_outcome <- function(res, label) {
  cat("\n========================================\n")
  cat(label, "\n")
  cat("========================================\n\n")

  cat("--- Test-set performance (Accuracy / AUC / logLoss) ---\n")
  print(as.data.frame(res$metrics))

  cat("\n--- Full GLM: significant predictors (p < .05) ---\n")
  full_sig <- res$coefs %>%
    dplyr::filter(model == "Full", significant == TRUE, term != "(Intercept)") %>%
    dplyr::select(term, estimate, p.value, conf.low, conf.high)
  if (nrow(full_sig) == 0) cat("  (none)\n") else print(as.data.frame(full_sig))

  cat("\n--- Stepwise GLM: retained predictors (any p-value) ---\n")
  step_coefs <- res$coefs %>%
    dplyr::filter(model == "Stepwise", term != "(Intercept)") %>%
    dplyr::select(term, estimate, p.value, significant)
  if (nrow(step_coefs) == 0) cat("  (none retained)\n") else print(as.data.frame(step_coefs))

  cat("\n--- Stepwise: bootstrap selection stability (top 10 by selected_pct) ---\n")
  step_boot <- res$boot_ci %>%
    dplyr::filter(model == "Stepwise", term != "(Intercept)") %>%
    dplyr::arrange(dplyr::desc(selected_pct)) %>%
    dplyr::select(term, mean, ci_lower, ci_upper, selected_pct) %>%
    head(10)
  if (nrow(step_boot) == 0) cat("  (no bootstrap data)\n") else print(as.data.frame(step_boot))

  cat("\n--- ElasticNet: non-zero coefficients (PRIMARY RESULT) ---\n")
  pen_coef <- res$coefs %>%
    dplyr::filter(model == "ElasticNet", term != "(Intercept)") %>%
    dplyr::select(term, estimate) %>%
    dplyr::arrange(dplyr::desc(abs(estimate)))
  if (nrow(pen_coef) == 0) cat("  (all shrunk to zero)\n") else print(as.data.frame(pen_coef))

  cat("\n--- ElasticNet: bootstrap stability (top 10 by nonzero_pct) ---\n")
  pen_boot <- res$boot_ci %>%
    dplyr::filter(model == "ElasticNet", term != "(Intercept)") %>%
    dplyr::arrange(dplyr::desc(nonzero_pct)) %>%
    dplyr::select(term, mean, ci_lower, ci_upper, nonzero_pct) %>%
    head(10)
  if (nrow(pen_boot) == 0) cat("  (no bootstrap data)\n") else print(as.data.frame(pen_boot))

  cat("\n")
}

summarize_outcome(res_AD,  "ADVANCE DIRECTIVE (AD)")
summarize_outcome(res_HP,  "HEALTHCARE PROXY (HP)")
summarize_outcome(res_ANY, "BOTH DIRECTIVES (AD and HP)")

