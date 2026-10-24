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
#
# Version 2.1 (this edit):
#   - The original Full GLM / Stepwise GLM (MASS::stepAIC) coefficient tables showed classic
#     (quasi-)complete separation: coefficients in the tens-to-hundreds, SEs 1000x+ larger than
#     the estimate, and p-values pinned near 0.999 across the board. That is the ordinary-MLE
#     logistic regression failing to converge to a finite optimum, not evidence of weak effects.
#   - "Full GLM" and "Stepwise GLM" are now fit with Firth's penalized likelihood
#     (package `logistf`) instead of ordinary glm()/MASS::stepAIC(). Firth's penalty keeps
#     coefficients and SEs finite under separation and gives valid profile-likelihood-based CIs
#     and p-values in exactly this situation, without ElasticNet's shrink-everything-toward-zero
#     behaviour. "Stepwise" now means logistf::backward() (penalized-likelihood-ratio-based
#     backward elimination) rather than AIC-based stepAIC.
#   - Model labels "Full" and "Stepwise" are kept as-is throughout the script (metrics, coefs,
#     boot_ci, plots, summarize_outcome) so downstream filtering/reporting code did not need to
#     change - just be aware both now refer to Firth-penalized fits, not ordinary glm().
#   - ElasticNet is unchanged.
#
# Version 2.2 (this edit):
#   - Dropped "Stepwise" (logistf::backward()) entirely. It repeatedly hit new low-level numerical
#     failures as different candidate submodels were tried during elimination - PL-CI
#     non-convergence, "selecting among factor variables not supported", and (twice, for different
#     predictor subsets) LAPACK rank-deficiency errors inside logistf.fit(). Each was individually
#     fixable, but the pattern itself is the real finding: backward selection over ~22 candidate
#     submodels is inherently unstable at N=94 with several near-collinear/rare-category
#     predictors, Firth-penalized or not. Patching each crash doesn't fix that instability.
#   - The pipeline now compares TWO models per outcome: Full GLM (Firth, all predictors, no
#     selection) and ElasticNet (unchanged) - not three. All Stepwise-specific code
#     (fit_step/backward() calls, bootstrap_firth_coefs(), the "Stepwise" bootstrap CI block, and
#     the corresponding summarize_outcome() section) has been removed rather than left dead/unused.
#   - Full GLM's CIs still come from Firth's own profile-likelihood interval (extract_firth_coefs())
#     - no bootstrap needed for Full, matching the original design's treatment of Full vs Stepwise.
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
  "logistf",
  "randomForest",
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
# 9. Modeling helpers: Full GLM (Firth) / Elastic Net / Random Forest (accuracy benchmark only)
# ================================================================================================
# These replace functionsUsed.R's results_modelSAFETYmod()/print_AUC(), which were built for an
# unrelated study and dataset. Everything the pipeline needs now lives in this script.

# logistf::backward()/add1.logistf() cannot select among multi-level factor predictors directly
# ("Selecting among factor variables is not supported") - they need already-numeric (0/1
# dummy-coded) columns to add/drop one at a time, the same representation glmnet/caret already
# use internally for the ElasticNet model. This helper expands factor predictors to their 0/1
# dummy columns via model.matrix() (same encoding findLinearCombos()/the correlation check used
# earlier), and rebuilds the formula against those dummy column names instead of the original
# factor names. Column names from model.matrix() can contain characters that aren't valid bare R
# names (e.g. "nationality_dichNon-German" has a hyphen) - check.names = FALSE preserves the exact
# names so they still match between the data.frame and the backtick-quoted formula terms.
make_dummy_design <- function(data, outcome, predictors) {
  # Filter to complete cases FIRST, explicitly - model.matrix() silently drops any row with an NA
  # in a predictor via its default na.action, which previously caused data_dummy[[outcome]] <-
  # data[[outcome]] to fail with a row-count mismatch (X had fewer rows than the original `data`
  # whenever a predictor had a missing value, e.g. one NA row in a test set: 18 rows in, 17 out).
  # Filtering here up front means model.matrix() has nothing left to silently drop, so row counts
  # always match by construction.
  data_cc <- data[stats::complete.cases(data[, c(outcome, predictors)]), , drop = FALSE]

  rhs_form <- as.formula(paste("~", paste(predictors, collapse = " + ")))
  X <- model.matrix(rhs_form, data = data_cc)[, -1, drop = FALSE]  # drop intercept column

  # Sanitize to syntactically valid R names (e.g. "nationality_dichNon-German" ->
  # "nationality_dichNon.German") instead of relying on backtick-quoting in the formula text.
  # Backtick-quoting worked for simple cases but did not survive logistf's internal
  # model.frame()/predict() reconstruction cleanly ("undefined columns selected" errors) - R
  # model internals are inconsistent about round-tripping non-syntactic names even when
  # backtick-quoted. Sanitizing up front removes the whole problem instead of chasing it further.
  dummy_names_original <- colnames(X)
  dummy_names_safe      <- make.names(dummy_names_original, unique = TRUE)
  colnames(X) <- dummy_names_safe

  data_dummy <- as.data.frame(X, check.names = FALSE)
  data_dummy[[outcome]] <- data_cc[[outcome]]

  form_dummy <- as.formula(paste(outcome, "~", paste(dummy_names_safe, collapse = " + ")))

  # name_map translates a sanitized column name (what will show up in coef()/formula(fit)) back
  # to the original, readable dummy label (e.g. "nationality_dichNon.German" ->
  # "nationality_dichNon-German") - used by extract_firth_coefs() so coefficient tables show
  # readable term names despite the internal sanitization.
  name_map <- stats::setNames(dummy_names_original, dummy_names_safe)

  list(data = data_dummy, predictors = dummy_names_safe, formula = form_dummy, name_map = name_map)
}

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

  # Full GLM is fit with Firth's penalized likelihood (logistf) instead of ordinary glm(), which
  # broke down under (quasi-)complete separation given ~20 predictors vs N=94. See the version 2.1
  # note at the top of this script.
  # complete.cases() up front: logistf doesn't take an na.action argument the way caret::train
  # does, so NA handling has to be done explicitly here to match the complete-case behaviour
  # na.omit gave the other models.
  train_data_cc <- train_data[stats::complete.cases(train_data[, c(outcome, predictors)]), ]

  # Dummy-code factor predictors - keeps term naming/encoding identical to what ElasticNet already
  # uses internally, and is retained from the (now-removed) Stepwise/backward() code path since it
  # does no harm here and keeps coefficient naming consistent between Full and ElasticNet.
  dd_train <- make_dummy_design(train_data_cc, outcome, predictors)

  fit_full <- tryCatch(
    logistf::logistf(
      dd_train$formula, data = dd_train$data, pl = TRUE, na.action = na.omit,
      # Default control maxit = 25 for the main Newton-Raphson fit was too low too - "Maximum
      # number of iterations for full model exceeded" was seen in the last run alongside the
      # PL-CI non-convergence below. Both are flat-likelihood symptoms, not bugs: with ~22
      # predictors vs N=94, both the point-estimate fit and the CI search legitimately need more
      # room. If warnings persist for some terms even at these higher maxit values, that is itself
      # informative (those effects are close to unidentifiable at this N) rather than something to
      # keep chasing with ever-higher maxit.
      control   = logistf::logistf.control(maxit = 1000),
      plcontrol = logistf::logistpl.control(maxit = 1000)
    ),
    error = function(e) {
      stop(sprintf(paste(
        "logistf() failed for outcome '%s' with:\n  %s\n",
        "If this mentions 'not positive definite' / LAPACK dpotrf, the predictor matrix is likely",
        "rank-deficient (exact or near-exact linear dependency between predictors), not the usual",
        "separation problem Firth otherwise handles. Diagnose with:\n",
        "  X <- model.matrix(%s, data = train_data_cc)[, -1, drop = FALSE]\n",
        "  caret::findLinearCombos(X)   # identifies the dependent column(s) to drop/merge\n",
        "(This happened once already for advance_directive_dich, caused by housing_dich being an",
        "exact duplicate of cohabitation_dichInstitutional in this dataset - see predictors_common",
        "comment in section 10 - and was fixed by dropping housing_dich. A different outcome or",
        "predictor set can still hit a different rank-deficiency and will need the same diagnosis.)",
        sep = "\n"
      ), outcome, conditionMessage(e), deparse(dd_train$formula)), call. = FALSE)
    }
  )

  grid_en <- expand.grid(alpha = alpha_grid, lambda = lambda_grid)
  fit_pen <- caret::train(
    form, data = train_data, method = "glmnet", family = "binomial",
    trControl = train_control, tuneGrid = grid_en, metric = "logLoss", na.action = na.omit
  )

  # Random Forest - added purely as an accuracy-comparison benchmark (not for coefficient
  # interpretation - RF has no coefficients/p-values, so it is deliberately NOT added to `coefs`/
  # `boot_ci`/plot_significant_effects(); only to `metrics`/roc_list in run_pipeline()). Motivation:
  # unlike Full/ElasticNet, RF captures interactions and non-linearities implicitly without them
  # being specified - useful as a check on whether a flexible, interaction-capturing model beats
  # the additive models by much on held-out AUC. At N=94 it's also easy for RF to overfit, so a
  # "wins" comparison should be read as a diagnostic hint (worth chasing down which interaction
  # might matter, e.g. via glinternet/hierNet), not as a signal to replace the additive models.
  # method = "rf" handles factors directly like glmnet's formula interface does - no dummy-coding
  # needed. No tuneGrid given - caret's default mtry grid (a handful of values scaled to predictor
  # count) is adequate here; metric = "logLoss" matches train_control's mnLogLoss summaryFunction,
  # same convention as fit_pen above - keep both in sync if you change train_control's
  # summaryFunction (e.g. to defaultSummary/twoClassSummary, as discussed for the LOOCV issue).
  fit_rf <- caret::train(
    form, data = train_data, method = "rf",
    trControl = train_control, metric = "logLoss", na.action = na.omit, importance = TRUE
  )

  list(full = fit_full, pen = fit_pen, rf = fit_rf, name_map = dd_train$name_map,
       model_data = model_data, train_data = train_data_cc, test_data = test_data, outcome = outcome)
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

# Same metrics as evaluate_test_metrics(), but for logistf objects (Firth Full / Firth Backward).
# logistf fits aren't caret::train objects - predict.logistf() returns probabilities directly
# rather than via type = "prob"/"raw", so this needs its own small function instead of reusing
# evaluate_test_metrics() as-is.
# Same metrics as evaluate_test_metrics(), but for logistf objects (Firth Full / Firth Backward).
# logistf fits aren't caret::train objects - predict.logistf() returns probabilities directly
# rather than via type = "prob"/"raw", so this needs its own small function instead of reusing
# evaluate_test_metrics() as-is. `predictors` (original, non-dummy) is needed here so test_data can
# be dummy-coded the same way training data was in fit_three_models() - the fitted model's formula
# is expressed in terms of dummy columns (see make_dummy_design()), so newdata must match.
evaluate_firth_metrics <- function(fit, test_data, outcome, predictors) {
  dd_test      <- make_dummy_design(test_data, outcome, predictors)
  model_vars   <- all.vars(formula(fit))[-1]
  test_data_cc <- dd_test$data[stats::complete.cases(dd_test$data[, model_vars]), ]

  probs <- predict(fit, newdata = test_data_cc, type = "response")
  preds <- factor(ifelse(probs > 0.5, "Yes", "No"), levels = c("No", "Yes"))
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

# Combined ROC plot - works for any named list of roc objects (Full / ElasticNet / RandomForest),
# with AUC + 95% CI in the legend
plot_auc_comparison <- function(roc_list, title) {
  legend_labels <- vapply(names(roc_list), function(nm) {
    auc_val <- as.numeric(pROC::auc(roc_list[[nm]]))
    # ci.auc() can fail ("need at least two non-NA values to interpolate") when a model hits
    # perfect/near-perfect separation - all predicted probabilities collapse toward 0/1, giving
    # a degenerate ROC curve with no interpolation room. With ~20 predictors and only ~75
    # training rows, this is a plausible outcome for the Full model, not a bug - it's itself a
    # sign that model is overfit/unstable. Caught here so one degenerate model doesn't crash the
    # whole pipeline; the AUC point estimate is still shown.
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

# All terms, p-values, and penalized-likelihood 95% CIs for the Firth Full-GLM fit.
# Replaces the old extract_glm_coefs(), which pulled these from ordinary glm() via broom::tidy()
# and needed a tryCatch fallback because profile-likelihood CIs failed under separation. logistf
# objects carry valid CIs/p-values directly (fit$ci.lower/ci.upper/prob) precisely because Firth's
# penalty is designed not to break down in that situation - no fallback branch needed here.
# name_map (from make_dummy_design()) is optional and, if supplied, translates sanitized dummy
# term names (e.g. "nationality_dichNon.German") back to their original, readable form (e.g.
# "nationality_dichNon-German") for display - the model itself is unaffected either way.
extract_firth_coefs <- function(fit, model_label, alpha = 0.05, name_map = NULL) {
  term_names <- names(coef(fit))
  if (!is.null(name_map)) {
    term_names <- dplyr::if_else(term_names %in% names(name_map),
                                  name_map[term_names], term_names)
  }
  data.frame(
    term        = term_names,
    estimate    = as.numeric(coef(fit)),
    std.error   = sqrt(diag(fit$var)),
    p.value     = as.numeric(fit$prob),
    conf.low    = as.numeric(fit$ci.lower),
    conf.high   = as.numeric(fit$ci.upper),
    model       = model_label,
    significant = as.numeric(fit$prob) < alpha,
    row.names   = NULL
  )
}

# Non-zero terms for the Elastic Net fit at its cross-validated best lambda
extract_glmnet_coefs <- function(fit, model_label) {
  coefs <- as.matrix(coef(fit$finalModel, s = fit$bestTune$lambda))
  data.frame(term = rownames(coefs), estimate = coefs[, 1], model = model_label) %>%
    dplyr::filter(estimate != 0)
}

# Variable importance for the Random Forest fit. Unlike extract_firth_coefs()/extract_glmnet_coefs(),
# this has no p-value or direction of effect (positive/negative) - importance measures how much a
# predictor helps prediction accuracy overall, not which way it pushes the outcome. Kept as its own
# table (res$rf_importance) rather than merged into `coefs`, which downstream code
# (plot_significant_effects()/summarize_outcome()) assumes has a signed estimate + p-value.
extract_rf_importance <- function(fit, model_label = "RandomForest") {
  vi <- caret::varImp(fit)$importance
  # caret::varImp() for classification RF sometimes returns one importance column per class level
  # (e.g. "No", "Yes") instead of a single "Overall" column, depending on the importance type
  # available - average across columns in that case, giving one overall score per predictor
  # regardless of which format was returned, rather than silently reading just the first class.
  imp_val <- if ("Overall" %in% colnames(vi)) vi[["Overall"]] else rowMeans(vi)
  data.frame(term = rownames(vi), importance = imp_val, model = model_label, row.names = NULL) %>%
    dplyr::arrange(dplyr::desc(importance))
}

# --------------------------------------------------------------------------------------------------
# Bootstrap coefficient stability / 95% CIs (default 1000 resamples), applied to ElasticNet only -
# matches the published methods (Full GLM's CI comes from Firth's own parametric
# penalized-likelihood interval in extract_firth_coefs() above instead, so no bootstrap is needed
# for it). Resampling is done WITH replacement from the FULL dataset (not just the 80% training
# partition) - confirmed by the results_bootstrap() call in the SAFETY study, which passes the
# full pre-split data. The 80/20 split is reserved for the Accuracy/AUC/logLoss performance
# comparison; this bootstrap is a separate analysis aimed at coefficient stability/uncertainty
# using all available data. Each resample is refit at the alpha/lambda already chosen via
# cross-validation on the real (non-bootstrapped) data.
#
# (Version 2.2 removed the "Stepwise"/logistf::backward()-based bootstrap that used to live here -
# see the version 2.2 note at the top of this script for why Stepwise was dropped entirely.)
# --------------------------------------------------------------------------------------------------

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

# Runs both models for one outcome: fits, evaluates (Accuracy/AUC/logLoss), plots ROC, extracts
# coefficients (with parametric CIs for Full), and bootstraps CIs for ElasticNet coefficients.
# NOTE: "Full" below is a Firth-penalized fit (logistf), not ordinary glm() - see the version 2.1
# note at the top of this script. The model LABEL is unchanged so metrics/coefs/boot_ci/plots/
# summarize_outcome downstream don't need to know or care which fitting method produced it.
# "Stepwise" was removed entirely in version 2.2 - see that note for why.
run_pipeline <- function(outcome, predictors, data, train_control, label, B = 1000) {

  fits <- fit_three_models(outcome, predictors, data, train_control)

  m_full <- evaluate_firth_metrics(fits$full, fits$test_data, outcome, predictors)
  m_pen  <- evaluate_test_metrics(fits$pen,  fits$test_data, outcome)
  m_rf   <- evaluate_test_metrics(fits$rf,   fits$test_data, outcome)

  metrics <- data.frame(
    model    = c("Full", "ElasticNet", "RandomForest"),
    Accuracy = c(m_full$Accuracy, m_pen$Accuracy, m_rf$Accuracy),
    AUC      = c(m_full$AUC,      m_pen$AUC,      m_rf$AUC),
    logLoss  = c(m_full$logLoss,  m_pen$logLoss,  m_rf$logLoss)
  )

  roc_list <- list(Full = m_full$roc, ElasticNet = m_pen$roc, RandomForest = m_rf$roc)
  fig <- plot_auc_comparison(roc_list, title = paste("ROC comparison -", label))

  coefs <- dplyr::bind_rows(
    extract_firth_coefs(fits$full, "Full", name_map = fits$name_map),
    extract_glmnet_coefs(fits$pen, "ElasticNet")
  )
  # RandomForest is NOT added to `coefs` - it has no coefficients/p-values/direction of effect,
  # only variable importance (how much a predictor helps prediction accuracy, not whether it
  # increases or decreases the outcome). Keeping it separate avoids it being silently picked up by
  # plot_significant_effects()/summarize_outcome(), which assume a signed coefficient + p-value.
  rf_importance <- extract_rf_importance(fits$rf)

  # Bootstrap from the FULL dataset (model_data), not just the training partition - see note
  # above bootstrap_glmnet_coefs(). Cached to disk since 1000 resamples is slow; delete the
  # .Rdata file in results/ to force a recompute (e.g. after changing predictors).
  # Cache filename carries a "_firth" tag so any .Rdata files left over from the old
  # stepAIC-based bootstrap (version <= 2.0) are NOT silently reloaded here. Since version 2.2
  # also removed the Stepwise bootstrap, any "_firth" cache from version 2.1 is ALSO now stale
  # (it contains a boot_step object this version no longer produces or expects) - delete existing
  # CIdataBootstrap_<outcome>_firth.Rdata files once after upgrading to this version.
  cache_file <- file.path(wdir, "results", paste0("CIdataBootstrap_", outcome, "_firth.Rdata"))
  if (file.exists(cache_file)) {
    load(cache_file)   # restores boot_pen
  } else {
    boot_pen <- bootstrap_glmnet_coefs(
      fits$model_data, outcome, predictors,
      alpha  = fits$pen$bestTune$alpha,
      lambda = fits$pen$bestTune$lambda,
      B = B
    )
    save(boot_pen, file = cache_file)
  }
  boot_ci <- dplyr::mutate(boot_pen, model = "ElasticNet")

  list(fits = fits, metrics = metrics, roc = roc_list, plot = fig, coefs = coefs,
       boot_ci = boot_ci, rf_importance = rf_importance)
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
# housing_dich REMOVED (see note below) - was found to be exactly collinear with
# cohabitation_dichInstitutional (r = 1.00, caret::findLinearCombos() column pair 6/8), which
# caused logistf.fit() to fail with "LAPACK dpotrf failed: matrix is not positive definite".
# Root cause: both are derived from Housing_situation - housing_dich = "No" means
# Housing_situation %in% 1:4 (anything but "own home"), cohabitation_dich = "Institutional" means
# the narrower Housing_situation %in% c(2,3) (assisted living/nursing home). They should NOT be
# identical in general (Institutional is meant to be a strict subset of housing_dich == "No"),
# but with only 4 people total in that group in this dataset, none of them happen to be coded
# "With children" (1) or "Other" (4) - so the two columns collapse to be numerically identical
# here. cohabitation_dich carries strictly more information (alone / with someone / institutional)
# than housing_dich (own home / not), so nothing is lost by dropping housing_dich instead of
# cohabitation_dich. Re-run caret::findLinearCombos() on the full predictor matrix after this
# change to confirm no other rank-deficiencies remain before refitting.
predictors_common <- c(
  "age", "gender_dich", "nationality_dich", "marital_status_dich",
  "spiritual_orientation_dich", "cohabitation_dich",
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
# high-predictor/limited-N regime) as the PRIMARY result. Full GLM is still reported as
# EXPLORATORY - low EPP means wide CIs regardless of fitting method. As of version 2.1 it is fit
# with Firth's penalized likelihood (see the note at the top of this script) specifically because
# the ORIGINAL ordinary glm() version hit (quasi-)complete separation here - coefficients in the
# tens-to-hundreds, SEs 1000x+ the estimate, p-values pinned near .999. Firth's correction fixes
# the separation/non-convergence problem itself (coefficients and CIs are now finite and usable),
# but it does NOT fix the underlying low-EPP instability. ("Stepwise" was dropped entirely in
# version 2.2 after repeatedly hitting new numerical failures during backward elimination - see
# that note - so there is no bootstrap selection-instability column for it anymore; the
# `nonzero_pct` column in boot_ci, below, still quantifies ElasticNet's selection instability
# across 1000 resamples.)
message(paste(strwrap(
  "Sample size caveat: N=94, ~20-21 predictors per model. Full GLM (Firth-penalized as of v2.1)
   results are exploratory due to low events-per-predictor; ElasticNet is the primary result.
   Stepwise selection was dropped (v2.2) due to repeated numerical instability during backward
   elimination at this N.",
  width = 100), collapse = "\n"))

# B_BOOT controls the number of bootstrap resamples used for the ElasticNet coefficient CIs. The
# published methods used 1000. Each run_pipeline() call fits roughly B_BOOT extra models on top of
# the CV-tuned Full/ElasticNet fits, so this is by far the slowest part of the script. Lower it
# (e.g. to 100) for a quick test run, then set back to 1000 for the final analysis.
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
# 14. Export coefficient tables (with parametric CI for Full) + ROC plots
# ================================================================================================
write.csv(res_AD$coefs,  file.path(wdir, "results", "coefficients_AdvanceDirective.csv"), row.names = FALSE)
write.csv(res_HP$coefs,  file.path(wdir, "results", "coefficients_HealthcareProxy.csv"),  row.names = FALSE)
write.csv(res_ANY$coefs, file.path(wdir, "results", "coefficients_AnyDirective.csv"),     row.names = FALSE)

ggplot2::ggsave(file.path(wdir, "results", "ROC_AdvanceDirective.png"), res_AD$plot,  width = 6, height = 5)
ggplot2::ggsave(file.path(wdir, "results", "ROC_HealthcareProxy.png"), res_HP$plot,  width = 6, height = 5)
ggplot2::ggsave(file.path(wdir, "results", "ROC_AnyDirective.png"),    res_ANY$plot, width = 6, height = 5)

# ================================================================================================
# 15. Export bootstrap coefficient stability / 95% CI (ElasticNet only)
# ================================================================================================
write.csv(res_AD$boot_ci,  file.path(wdir, "results", "bootstrapCI_AdvanceDirective.csv"), row.names = FALSE)
write.csv(res_HP$boot_ci,  file.path(wdir, "results", "bootstrapCI_HealthcareProxy.csv"),  row.names = FALSE)
write.csv(res_ANY$boot_ci, file.path(wdir, "results", "bootstrapCI_AnyDirective.csv"),     row.names = FALSE)

# ================================================================================================
# 15b. Export Random Forest variable importance (accuracy-benchmark model - see fit_three_models())
# ================================================================================================
# Importance here means "how much this predictor helps prediction accuracy", NOT a signed
# coefficient - there is no "increases/decreases the outcome" direction to report, unlike coefs.csv.
write.csv(res_AD$rf_importance,  file.path(wdir, "results", "rfImportance_AdvanceDirective.csv"), row.names = FALSE)
write.csv(res_HP$rf_importance,  file.path(wdir, "results", "rfImportance_HealthcareProxy.csv"),  row.names = FALSE)
write.csv(res_ANY$rf_importance, file.path(wdir, "results", "rfImportance_AnyDirective.csv"),     row.names = FALSE)

# ================================================================================================
# 16. Print readable results summary to console for each outcome
# ================================================================================================
# Shows: test-set performance, significant Full predictors, ElasticNet non-zero coefficients,
# and bootstrap stability - all in one place per outcome.
# Remember: Full results are EXPLORATORY given N=94 vs ~21 predictors (see section 12 caveat).
# ElasticNet + bootstrap stability are the primary, trustworthy results.

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

# ================================================================================================
# 17. One plot per significant/retained predictor: raw relationship to its outcome
# ================================================================================================
# For every term flagged either by Full GLM (p < .05) or retained (non-zero) by ElasticNet, this
# produces ONE plot showing that predictor's relationship to the outcome directly - not the
# model's coefficient, but the actual data: predictor on the x-axis, outcome on the y-axis.
#   - Factor predictor (e.g. power_of_attorney_dich: No/Yes): bar chart of the PROPORTION with
#     outcome = "Yes" at each level of the predictor, with a normal-approximation 95% CI and the n
#     at each bar. E.g. "power_of_attorney_dich" (x: No/Yes) -> proportion with
#     "advance_directive_dich" = Yes (y), directly answering "does having power of attorney go
#     with also having an advance directive?"
#   - Numeric predictor (e.g. age, PDQ_score): the outcome (0/1) jittered against the predictor,
#     with a fitted logistic curve (geom_smooth, family = binomial) showing the predicted
#     probability of the outcome across the predictor's range.
# IMPORTANT INTERPRETIVE NOTE: these are RAW, UNADJUSTED (bivariate) relationships straight from
# the data - they do not hold other predictors constant the way the Full/ElasticNet coefficients
# do. That is deliberate: the goal here is graphical intuition for what a flagged term means in
# plain data terms, not a restatement of the adjusted model. A predictor can look weaker/stronger
# here than its adjusted coefficient suggests if it is correlated with other predictors in the
# model (see the correlation matrix from the collinearity check earlier in this project).
# Multi-level factors (e.g. cohabitation_dich: No/Yes/Institutional) show ALL levels together even
# if only one level's dummy term was flagged, since that is the more interpretable/complete view.

# Maps a (possibly dummy-coded, e.g. "power_of_attorney_dichYes") term name back to the original
# predictor column it came from (e.g. "power_of_attorney_dich"), by prefix-matching against the
# known predictor list for that outcome. Numeric predictors match themselves exactly.
term_to_predictor <- function(term, predictors) {
  matches <- predictors[vapply(predictors, function(p) startsWith(term, p), logical(1))]
  if (length(matches) == 0) return(NA_character_)
  matches[which.max(nchar(matches))]   # longest match, in case of any shared-prefix ambiguity
}

# Builds one plot for a single predictor -> outcome relationship, using the ORIGINAL (non-dummy)
# data so factors show all their levels together.
plot_predictor_effect <- function(data, predictor, outcome, outcome_label, predictor_label = predictor) {
  df <- data %>%
    dplyr::select(dplyr::all_of(c(predictor, outcome))) %>%
    tidyr::drop_na()
  df$.outcome_num <- as.numeric(df[[outcome]] == "Yes")

  if (is.factor(df[[predictor]]) || is.character(df[[predictor]])) {
    df[[predictor]] <- as.factor(df[[predictor]])
    summ <- df %>%
      dplyr::group_by(.data[[predictor]]) %>%
      dplyr::summarise(
        n     = dplyr::n(),
        n_yes = sum(.outcome_num),
        prop  = n_yes / n,
        se    = sqrt(prop * (1 - prop) / n),
        .groups = "drop"
      ) %>%
      dplyr::mutate(ci_lower = pmax(0, prop - 1.96 * se),
                     ci_upper = pmin(1, prop + 1.96 * se))

    ggplot2::ggplot(summ, ggplot2::aes(x = .data[[predictor]], y = prop)) +
      ggplot2::geom_col(fill = "steelblue", width = 0.6) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), width = 0.15) +
      ggplot2::geom_text(ggplot2::aes(label = paste0("n=", n)), vjust = -1, size = 3.2) +
      ggplot2::scale_y_continuous(limits = c(0, 1.08)) +
      ggplot2::labs(
        title = paste0(predictor_label, "  ->  ", outcome_label),
        x = predictor_label,
        y = paste0("Proportion with ", outcome_label, " = Yes")
      ) +
      ggplot2::theme_minimal()

  } else {
    ggplot2::ggplot(df, ggplot2::aes(x = .data[[predictor]], y = .outcome_num)) +
      ggplot2::geom_jitter(height = 0.04, width = 0, alpha = 0.4, colour = "steelblue") +
      ggplot2::geom_smooth(method = "glm", method.args = list(family = "binomial"),
                            colour = "firebrick", fill = "firebrick", alpha = 0.15, na.rm = TRUE) +
      ggplot2::scale_y_continuous(limits = c(-0.08, 1.08), breaks = c(0, 1)) +
      ggplot2::labs(
        title = paste0(predictor_label, "  ->  ", outcome_label),
        x = predictor_label,
        y = paste0("P(", outcome_label, " = Yes)")
      ) +
      ggplot2::theme_minimal()
  }
}

# For one outcome: finds every Full-significant (p < .05) or ElasticNet-retained (non-zero) term,
# maps each back to its original predictor column, and builds one plot per unique predictor
# (de-duplicating multi-level factors that contributed more than one flagged dummy term).
plot_significant_effects <- function(res, predictors, outcome, outcome_label) {
  full_terms <- res$coefs %>%
    dplyr::filter(model == "Full", significant == TRUE, term != "(Intercept)") %>%
    dplyr::pull(term)
  pen_terms <- res$coefs %>%
    dplyr::filter(model == "ElasticNet", term != "(Intercept)") %>%
    dplyr::pull(term)

  sig_terms <- union(full_terms, pen_terms)
  if (length(sig_terms) == 0) {
    message(sprintf(
      "plot_significant_effects: no Full-significant or ElasticNet-retained terms for %s - nothing to plot.",
      outcome_label
    ))
    return(list())
  }

  predictor_cols <- unique(vapply(sig_terms, term_to_predictor, character(1), predictors = predictors))
  predictor_cols <- predictor_cols[!is.na(predictor_cols)]
  if (length(predictor_cols) < length(unique(sig_terms))) {
    unmatched <- setdiff(sig_terms, predictors[vapply(predictors, function(p)
      any(startsWith(sig_terms, p)), logical(1))])
    if (length(unmatched) > 0) {
      warning(sprintf(
        "plot_significant_effects (%s): could not map term(s) back to a source predictor column, skipped: %s",
        outcome_label, paste(unmatched, collapse = ", ")
      ))
    }
  }

  model_data <- res$fits$model_data   # original (non-dummy, all-levels) data used for this outcome
  plots <- list()
  for (p in predictor_cols) {
    plots[[p]] <- plot_predictor_effect(model_data, p, outcome, outcome_label, predictor_label = p)
  }
  plots
}

# Prints and saves each plot in a named list, e.g. results/Effect_AdvanceDirective_power_of_attorney_dich.png
save_effect_plots <- function(plots, outcome_tag) {
  for (nm in names(plots)) {
    print(plots[[nm]])
    ggplot2::ggsave(
      file.path(wdir, "results", sprintf("Effect_%s_%s.png", outcome_tag, nm)),
      plots[[nm]], width = 6, height = 4.5
    )
  }
}

effect_plots_AD  <- plot_significant_effects(res_AD,  predictors_AD,  "advance_directive_dich", "Advance Directive")
effect_plots_HP  <- plot_significant_effects(res_HP,  predictors_HP,  "power_of_attorney_dich", "Healthcare Proxy")
effect_plots_ANY <- plot_significant_effects(res_ANY, predictors_ANY, "both_directives_dich",   "Both Directives")

save_effect_plots(effect_plots_AD,  "AdvanceDirective")
save_effect_plots(effect_plots_HP,  "HealthcareProxy")
save_effect_plots(effect_plots_ANY, "AnyDirective")

# ================================================================================================
# 18. Repeated-holdout stability check: is a single train/test AUC reproducible, or a lucky split?
# ================================================================================================
# Motivation ("Both Directives" specifically): a single 80/20 split gave Full GLM (Firth) test-set
# AUC = 0.838, but ElasticNet's CV curve showed logLoss improving MONOTONICALLY all the way to full
# shrinkage for every alpha (no interior minimum - see res_ANY$fits$pen$results), and RandomForest
# also underperformed Full on that same split (AUC 0.675 vs 0.838). Both are hints that 0.838 might
# be a property of that one particular split rather than a reproducible effect - "Both Directives"
# is the smallest/most restricted outcome (requires AD AND HP both = Yes), so its test partition is
# probably only ~15-18 people, small enough for a 22-parameter Firth fit to land on a
# good-looking-by-chance AUC on any single split.
#
# This refits Full GLM (Firth) on many independent random 80/20 splits of the same data/outcome and
# reports the DISTRIBUTION of test-set AUCs, not just the one run_pipeline() happened to draw. If
# AUC swings widely across splits (regularly dropping toward/below 0.5), that confirms the
# single-split number wasn't reliable - consistent with what the ElasticNet CV curve already
# suggested. If it stays consistently high, that's a genuinely more surprising finding worth
# investigating further (and would argue AGAINST what the CV curve suggested, so worth flagging
# either way it comes out).
#
# Applied to all three outcomes (not just "Both Directives") for comparison - is this instability
# specific to "Both Directives", or does Full GLM swing similarly for AD/HP too?

repeated_holdout_auc <- function(outcome, predictors, data, n_reps = 50, p_split = 0.8, seed = 1) {
  set.seed(seed)
  model_data <- data %>%
    dplyr::select(dplyr::all_of(c(outcome, predictors))) %>%
    tidyr::drop_na(dplyr::all_of(outcome))

  results    <- vector("list", n_reps)
  error_msgs <- character(0)

  for (i in seq_len(n_reps)) {
    results[[i]] <- tryCatch({
      index         <- caret::createDataPartition(model_data[[outcome]], p = p_split, list = FALSE)
      train_data    <- model_data[index, ]
      test_data     <- model_data[-index, ]
      train_data_cc <- train_data[stats::complete.cases(train_data[, c(outcome, predictors)]), ]

      # pl = FALSE: only point predictions are needed here, not confidence intervals - keeps this
      # fast across many repeats (this is the same simplification bootstrap_glmnet_coefs()/the old
      # bootstrap_firth_coefs() used for their resampling loops, for the same reason).
      dd_train <- make_dummy_design(train_data_cc, outcome, predictors)
      fit_i <- logistf::logistf(
        dd_train$formula, data = dd_train$data, pl = FALSE, na.action = na.omit,
        control = logistf::logistf.control(maxit = 1000)
      )
      m_i <- evaluate_firth_metrics(fit_i, test_data, outcome, predictors)
      data.frame(rep = i, AUC = m_i$AUC, Accuracy = m_i$Accuracy, logLoss = m_i$logLoss)
    }, error = function(e) {
      error_msgs[[length(error_msgs) + 1]] <<- conditionMessage(e)
      NULL
    })
  }

  out <- dplyr::bind_rows(results)
  n_success <- nrow(out)
  message(sprintf("repeated_holdout_auc (%s): %d/%d splits succeeded.", outcome, n_success, n_reps))
  if (n_success < n_reps) {
    message(sprintf(
      "  %d split(s) failed (e.g. separation/rank-deficiency specific to that split - the same",
      n_reps - n_success
    ))
    message("  kind of failure fit_three_models()'s tryCatch reports for the single 'real' split) - excluded from summary.")
  }
  out
}

plot_repeated_holdout_auc <- function(auc_df, label, single_split_auc = NULL) {
  p <- ggplot2::ggplot(auc_df, ggplot2::aes(x = AUC)) +
    ggplot2::geom_histogram(binwidth = 0.05, boundary = 0, fill = "steelblue", colour = "white") +
    ggplot2::geom_vline(xintercept = 0.5, colour = "grey40", linetype = "dotted") +
    ggplot2::labs(
      title = paste0("Full GLM (Firth) test-set AUC across ", nrow(auc_df), " random 80/20 splits"),
      subtitle = label,
      x = "Test-set AUC", y = "Count"
    ) +
    ggplot2::theme_minimal()

  if (!is.null(single_split_auc)) {
    p <- p +
      ggplot2::geom_vline(xintercept = single_split_auc, colour = "firebrick",
                           linetype = "dashed", linewidth = 1) +
      ggplot2::labs(subtitle = paste0(label, sprintf(
        " - dashed red line = original single-split AUC (%.3f)", single_split_auc
      )))
  }
  p
}

N_HOLDOUT_REPS <- 50   # lower for a quick test run; 50+ recommended for a stable-looking distribution

holdout_AD  <- repeated_holdout_auc("advance_directive_dich", predictors_AD,  eol, n_reps = N_HOLDOUT_REPS)
holdout_HP  <- repeated_holdout_auc("power_of_attorney_dich", predictors_HP,  eol, n_reps = N_HOLDOUT_REPS)
holdout_ANY <- repeated_holdout_auc("both_directives_dich",   predictors_ANY, eol, n_reps = N_HOLDOUT_REPS)

cat("\n--- Repeated-holdout Full GLM (Firth) AUC: summary across", N_HOLDOUT_REPS, "splits ---\n")
cat("\nAdvance Directive:\n");  print(summary(holdout_AD$AUC))
cat("\nHealthcare Proxy:\n");   print(summary(holdout_HP$AUC))
cat("\nBoth Directives:\n");    print(summary(holdout_ANY$AUC))

fig_holdout_AD  <- plot_repeated_holdout_auc(holdout_AD,  "Advance Directive",
                                              single_split_auc = res_AD$metrics$AUC[res_AD$metrics$model == "Full"])
fig_holdout_HP  <- plot_repeated_holdout_auc(holdout_HP,  "Healthcare Proxy",
                                              single_split_auc = res_HP$metrics$AUC[res_HP$metrics$model == "Full"])
fig_holdout_ANY <- plot_repeated_holdout_auc(holdout_ANY, "Both Directives",
                                              single_split_auc = res_ANY$metrics$AUC[res_ANY$metrics$model == "Full"])

print(fig_holdout_AD)
print(fig_holdout_HP)
print(fig_holdout_ANY)

ggplot2::ggsave(file.path(wdir, "results", "RepeatedHoldoutAUC_AdvanceDirective.png"), fig_holdout_AD,  width = 7, height = 5)
ggplot2::ggsave(file.path(wdir, "results", "RepeatedHoldoutAUC_HealthcareProxy.png"),  fig_holdout_HP,  width = 7, height = 5)
ggplot2::ggsave(file.path(wdir, "results", "RepeatedHoldoutAUC_BothDirectives.png"),   fig_holdout_ANY, width = 7, height = 5)

write.csv(holdout_AD,  file.path(wdir, "results", "RepeatedHoldoutAUC_AdvanceDirective.csv"), row.names = FALSE)
write.csv(holdout_HP,  file.path(wdir, "results", "RepeatedHoldoutAUC_HealthcareProxy.csv"),  row.names = FALSE)
write.csv(holdout_ANY, file.path(wdir, "results", "RepeatedHoldoutAUC_BothDirectives.csv"),   row.names = FALSE)

# One combined boxplot showing the summary() quartiles (min/Q1/median/Q3/max, i.e. exactly what
# summary(holdout_AD$AUC) etc. printed above as text) for all three outcomes side by side, instead
# of reading three separate histograms. Each outcome's original single-split AUC (the number
# run_pipeline() reported once) is overlaid as a red diamond, so you can see at a glance whether
# that single number sat near the middle of its distribution (reproducible) or out in a tail
# (a lucky/unlucky split).
holdout_summary_df <- dplyr::bind_rows(
  dplyr::mutate(holdout_AD,  outcome = "Advance Directive"),
  dplyr::mutate(holdout_HP,  outcome = "Healthcare Proxy"),
  dplyr::mutate(holdout_ANY, outcome = "Both Directives")
)
single_split_df <- data.frame(
  outcome = c("Advance Directive", "Healthcare Proxy", "Both Directives"),
  AUC     = c(res_AD$metrics$AUC[res_AD$metrics$model == "Full"],
              res_HP$metrics$AUC[res_HP$metrics$model == "Full"],
              res_ANY$metrics$AUC[res_ANY$metrics$model == "Full"])
)

fig_holdout_summary <- ggplot2::ggplot(holdout_summary_df, ggplot2::aes(x = outcome, y = AUC)) +
  ggplot2::geom_boxplot(fill = "steelblue", alpha = 0.4, outlier.shape = NA, width = 0.5) +
  ggplot2::geom_jitter(width = 0.08, height = 0, alpha = 0.35, colour = "steelblue4") +
  ggplot2::geom_point(data = single_split_df, ggplot2::aes(x = outcome, y = AUC),
                       colour = "firebrick", shape = 18, size = 5) +
  ggplot2::geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey40") +
  ggplot2::labs(
    title = paste0("Full GLM (Firth) test-set AUC: ", N_HOLDOUT_REPS, " random 80/20 splits per outcome"),
    subtitle = "Red diamond = original single-split AUC. Dotted line = chance level (AUC = 0.5).",
    x = NULL, y = "Test-set AUC"
  ) +
  ggplot2::theme_minimal()

print(fig_holdout_summary)
ggplot2::ggsave(file.path(wdir, "results", "RepeatedHoldoutAUC_Summary_AllOutcomes.png"),
                 fig_holdout_summary, width = 8, height = 5.5)

