################################################################################
# Project: Parkinson's EOLC Planning - predictors of advance care planning
# Script:  AnalyseVorsorge_v3_2.R
# Authors: Anna and David Pedrosa
# Version: 3.2
################################################################################
#
# CONTENTS
#   1. PREAMBLE       1.1 Version history          1.2 User settings
#                     1.3 Packages                 1.4 Directories
#   2. PREPROCESSING  2.1 Data import              2.2 PDQ-39 subscales
#                     2.3 Codebook                 2.4 Dichotomisation / recoding
#                     2.5 Combined outcome (BOTH)  2.6 Rescaling & comorbidity
#                     2.7 Noise calibration columns
#   3. ANALYSES       3.1 Descriptive tables       3.2 Model helper functions
#                     3.3 Optimism-corrected AUC   3.4 ElasticNet stability
#                     3.5 Univariable screen       3.6 Per-outcome pipeline
#                     3.7 Predictor correlations   3.8 Fit all three outcomes
#                     3.9 Sensitivity analysis
#   4. REPORT         4.1 Figure functions         4.2 Save figures
#                     4.3 Performance summary      4.4 Interpretation guide
################################################################################


################################################################################
# 1. PREAMBLE                                                               ####
################################################################################


## 1.1 Version history =========================================================
# WHAT CHANGED vs. AnalyseVorsorge_09_07_26_Firth.R (v2.2)
# -----------------------------------------------------------
#
# 3. ElasticNet is NO LONGER reported as coefficients. Penalized coefficients are
#    shrunk toward zero by construction and have no valid SEs or p-values; percentile
#    bootstrap CIs on them do not achieve nominal coverage. It is used ONLY as a
#    SELECTION STABILITY screen: how often does each predictor survive with a
#    non-zero coefficient across bootstrap resamples? alpha and lambda are re-tuned
#    INSIDE each resample (v2.2 reused a single CV-selected lambda, which
#    understated exactly the selection uncertainty the bootstrap was meant to show).
#
# 4. RandomForest removed. At N=94 with 21 predictors it cannot be evaluated
#    reliably and contributed nothing interpretable.
#
# 5. The other directive (AD as predictor of PoA and vice versa) is EXCLUDED from
#    the primary models. Only 10 of 94 people have one document without the other
#    (OR ~ 80); conditioning on it makes every other coefficient mean "adjusted for
#    having already done the other half of the same thing". It is available as an
#    explicit sensitivity analysis (see SENSITIVITY_INCLUDE_OTHER_DIRECTIVE below).
#
# 6. Naming fixed: the third outcome is BOTH directives, not ANY. All objects and
#    output files now say BOTH.
#
# 7. All model fitting goes through a pre-built design matrix (build_design()), so
#    factor contrasts are identical across the main fit, every bootstrap resample,
#    and every prediction. No dummy-name sanitising / make_dummy_design() needed.
#
# Version 3.1 (this edit)
# -----------------------
# 8. NOISE CALIBRATION. Four pure-random columns (2 continuous, 2 binary at differing
#    prevalence) are carried through the ElasticNet stability screen. The highest
#    percentage any of them achieves is the "noise floor" - an empirical, dataset-specific
#    answer to "what selection percentage counts as meaningful?", replacing an arbitrary
#    borrowed threshold. Drawn as a dashed reference line on every stability plot.
#
# 9. REPORTING UNITS for continuous predictors (section 1.2.2, UNITS block). Odds ratios are
#    now expressed per 5 years of age, per 100 mg LEDD, per 10 UPDRS points, per 10 PDQ
#    points, per 5 BDI points, instead of per single unit. This is a change of units
#    only - identical model, identical p-values and AUC.
#
# 10. SEVERITY PROXY switched from UPDRS-III to LEDD, and Charlson dichotomised. See the
#     rationale in section 1.2.5; briefly, updrs_III correlates 0.08 with disease duration in
#     this cohort (52/94 are stimulated, so treated motor scores understate disease
#     stage), whereas LEDD correlates 0.54 and is balanced across DBS groups. Charlson is
#     0 in 69/94 people, so a continuous slope would rest on a handful of individuals.
#
# 11. CORRELATION HEATMAP (section 3.7) of the candidate predictors, upper triangle,
#     clustered. Needed to read the stability screen correctly: correlated predictors
#     split their selection probability between them, so a mediocre percentage can mean
#     "weak" OR "sharing signal with a partner", and only the correlation structure
#     distinguishes the two.


## 1.2 User settings - edit these, then source the whole file ==================

### 1.2.1 General --------------------------------------------------------------
PPR_INSTALL_PACKAGES <- TRUE

# Number of bootstrap resamples. 200 for a quick look, 1000 for the final analysis.
B_BOOT <- 1000

# Set TRUE to additionally fit the AD and PoA models WITH the other directive as a
# covariate (sensitivity analysis - see note 5 above). Does not affect primary output.
SENSITIVITY_INCLUDE_OTHER_DIRECTIVE <- TRUE

### 1.2.2 Reporting units for continuous predictors ----------------------------
# Rescaling sets the units in which the OR is expressed:
#
#   age        SD  8.7  -> 5 years  (~0.6 SD). Set to 10 if you prefer per-decade.
#   LEDD       SD 487   -> 100 mg   (~0.2 SD), the conventional unit in PD papers.
#   UPDRS_sum  SD 31    -> 5 points
#   PDQ_score  SD 13.7  -> 5 points
#   bdi_score  SD  4.6  -> 5 points
#   MOCA_score SD  3.2  -> 1 point (left unscaled)
#   duration   SD  5.9  -> 1 year  (left unscaled)
#
AGE_UNIT   <- 5     # years per unit of the age term
LEDD_UNIT  <- 100   # mg per unit of the LEDD term
UPDRS_UNIT <- 5    # UPDRS points per unit (approx. MCID)
PDQ_UNIT   <- 5    # PDQ points per unit (approx. MCID)
BDI_UNIT   <- 5     # BDI points per unit (approx. MCID)

### 1.2.3 Stability screen settings --------------------------------------------
# STABILITY_ALPHA: the elastic net mixing parameter. Only the LASSO component (alpha
#   towards 1) sets coefficients to exactly zero; the ridge component (alpha towards 0)
#   never does. Fixed at 1 = pure lasso, because this screen exists to RANK variables,
#   not to predict. Letting cross-validation choose alpha optimises prediction, and
#   prediction prefers ridge - which is how the first run ended up at alpha = 0.25 with
#   a 50% noise floor. Supply a vector here (e.g. c(0.5, 0.75, 1)) to go back to
#   CV-tuned alpha; the function will then report the chosen-alpha distribution.
STABILITY_ALPHA <- 1
#
# STABILITY_RESAMPLE: "subsample" draws n/2 WITHOUT replacement (Meinshausen &
#   Buhlmann's stability selection); "bootstrap" draws n WITH replacement. The bootstrap
#   duplicates ~63% of people, weakening the effective penalty and inflating every
#   selection percentage - which is also why published 0.6-0.9 thresholds do not
#   transfer to bootstrap numbers. Subsampling is the default.
STABILITY_RESAMPLE <- "subsample"
#
# STABILITY_LAMBDA: "lambda.1se" is the sparser of the two CV rules and the conventional
#   choice for stability selection; "lambda.min" keeps more variables and will raise the
#   noise floor.
STABILITY_LAMBDA <- "lambda.min"# "lambda.1se"

### 1.2.4 Noise calibration variables ------------------------------------------
# Rather than borrowing an arbitrary stability threshold (0.6? 0.7?) from a literature
# whose assumptions do not hold here, this generates pure-noise columns with no
# relationship to anything, carries them through the ElasticNet stability screen, and
# reports what selection percentage NOISE achieves in this dataset, at this N, with this
# correlation structure. That is the empirical floor: a real predictor has to beat it.
#
# Two continuous and two binary, because selection probability depends on the marginal
# distribution of a column - a rare binary and a normal continuous do not behave the
# same way under the lasso penalty, so a single noise type would be a misleading
# reference. The binary prevalences (0.40 / 0.15) bracket the real binary predictors,
# from nursing_support_dich (0.49) down to care_preference_dich (0.10).
#
# Generated ONCE with a fixed seed so all three outcomes see the SAME noise columns and
# their stability results are directly comparable.
NOISE_SEED <- 4711
noise_vars <- c("noise_normal_1", "noise_normal_2", "noise_binary_40", "noise_binary_15")

### 1.2.5 Pre-specified PRIMARY predictor set ----------------------------------
#
# Rationale for each:
#   age_u                      - advance care planning rises steeply with age in every
#                                published cohort; the single most established predictor.
#   gender_dich                - consistently reported sex difference in ACP uptake.
#   marital_status_dich        - a spouse is the usual proxy decision-maker; strong
#                                theoretical link to holding a power of attorney.
#                                (chosen over cohabitation_dich: same construct, but
#                                cohabitation's third level has only n=5.)
#   higher_education_dich      - health literacy / access to legal advice.
#   end_of_life_thoughts_dich  - the ACP-engagement variable; having thought about
#                                end-of-life wishes is the proximal step before
#                                documenting them.
#                                (chosen over sharing_thoughts_dich: the two agree on
#                                90% of people - 14/6/3/71 - so including both would
#                                split one effect across two terms.)
#   LEDD_u                     - disease severity proxy. Chosen over UPDRS-III after
#                                checking both against disease duration in this cohort:
#                                LEDD correlates 0.54 with duration overall and 0.68 /
#                                0.49 within the no-DBS / DBS strata, whereas updrs_III
#                                correlates 0.08 - i.e. the motor exam does NOT track
#                                progression here, because 52/94 patients are stimulated
#                                and their treated motor scores are lower despite ~3.5
#                                years longer disease. LEDD means are near-identical
#                                across DBS groups (734 vs 724 mg), so it is not badly
#                                confounded by stimulation the way UPDRS-III is.
#                                CAVEAT for the discussion: LEDD is still a TREATMENT
#                                variable and reflects prescribing practice and
#                                tolerability as well as disease burden.
#   charlson_dich              - comorbidity. Dichotomised rather than continuous
#                                because the raw score is 0 in 69 of 94 people, with
#                                counts of 12/8/3/1/1 at scores 1-6; a continuous slope
#                                would be leveraged largely by the single person scoring
#                                6, and assumes a constant log-odds change per point
#                                across a range with single-digit cell counts.
predictors_primary <- c(
  "age_u",
  "gender_dich",
  "marital_status_dich",
  "higher_education_dich",
  "death_preference_dich", #replaced care preference '
  # "end_of_life_thoughts_dich",
  "LEDD_u",
  "charlson_dich"
)

### 1.2.6 EXPLORATORY predictor set (sensitivity Firth fit) --------------------
# Deliberately NOT used for the primary p-values. housing_dich is excluded: it is
# numerically identical to cohabitation_dich == "Institutional" in this dataset
# (both isolate exactly the same 5 people), which makes the design matrix rank-deficient.
# Charlson_withoutage rather than Charlson_withage, since the latter bakes in age and
# would double-count it alongside the standalone age term.
# updrs_III is NOT included alongside UPDRS_u: they correlate at 0.88.
predictors_exploratory <- c(
  "age_u", "gender_dich", "nationality_dich", "marital_status_dich",
  "spiritual_orientation_dich", "cohabitation_dich",
  "nursing_support_dich", "rural_dich", "higher_education_dich",
  "end_of_life_thoughts_dich", "sharing_thoughts_dich",
  "care_preference_dich", "death_preference_dich",
  "duration", "dbs_dich", "LEDD_u", "UPDRS_u", "PDQ_u", "MOCA_score",
  "bdi_u", "charlson_dich"
)

### 1.2.7 STABILITY-SCREEN predictor set ---------------------------------------
# The exploratory set PLUS the noise columns. Kept separate so the noise variables
# calibrate the stability screen without also polluting the sensitivity Firth fit,
# which is already at EPP < 2 and does not need four more parameters.
predictors_stability <- c(predictors_exploratory, noise_vars)


## 1.3 Packages ================================================================
required_pkgs <- c("dplyr", "tidyr", "ggplot2", "readxl", "tableone",
                   "logistf", "glmnet", "pROC")

install_missing <- isTRUE(as.logical(Sys.getenv("PPR_INSTALL_PACKAGES", "FALSE")))
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  if (!install_missing) {
    stop(paste0("Missing package(s): ", paste(missing_pkgs, collapse = ", "),
                "\nInstall them, or set PPR_INSTALL_PACKAGES=TRUE."), call. = FALSE)
  }
  install.packages(missing_pkgs, dependencies = TRUE)
}
suppressPackageStartupMessages(for (p in required_pkgs) library(p, character.only = TRUE))


## 1.4 Directories =============================================================
username <- Sys.info()[["login"]]
if (username == "dpedr") {
  wdir <- "D:/EOL_predictors/";           data_dir <- file.path(wdir, "data")
} else if (username == "david") {
  wdir <- "/media/storage/EOL_predictors/"; data_dir <- file.path(wdir, "data")
} else {
  wdir <- "/Users/annapedrosa/Documents/dokumente/Arbeit_Anna/Aktuelle Projekte/StudiePoD_PD"
  data_dir <- wdir
}
setwd(wdir)
results_dir <- file.path(wdir, "results_v3_2")
dir.create(results_dir, showWarnings = FALSE)


################################################################################
# 2. PREPROCESSING                                                          ####
################################################################################

## 2.1 Data import =============================================================
eol <- readxl::read_excel(file.path(data_dir, "Datenmatrix_end-of-life_final_version11.10.23.xlsx")) # raw data

names(eol) <- trimws(names(eol))
eol <- eol %>% dplyr::rename(Thoughts_discussed_with = Thoughts_dicussed_with)

stopifnot(nrow(eol) == 94) # double_check integrity


## 2.2 PDQ-39 subscales (recoded from raw values) ==============================
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


## 2.3 Codebook and labelled display copy ======================================
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
  # CONFIRMED with the paper questionnaire: codes 3 and 4 are both rural categories.
  Residential_location = c("0" = "City", "1" = "Middle town", "2" = "Small town",
                           "3" = "Rural town/municipality", "4" = "Rural town/municipality"),
  Professional_education = c("0" = "None", "1" = "Recognised vocational training",
                             "2" = "Fachhochschule", "3" = "University", "4" = "Other"),
  Advance_directive = c("0" = "No", "1" = "Yes"),
  Power_of_attorney = c("0" = "No", "1" = "Yes"),
  thoughts_about_end_of_life_wishes = c("0" = "Never", "1" = "Rarely", "2" = "Now and then", "3" = "Frequently"),
  Sharing_of_thoughts = c("0" = "Never", "1" = "Rarely", "2" = "Now and then", "3" = "Frequently", "4" = "I don't know"),
  Thoughts_discussed_with = c("0" = "Partner", "1" = "Children", "2" = "GP", "3" = "Neurologist",
                              "4" = "Friend", "5" = "Other family member", "6" = "Other"),
  preferred_place_of_care = c("0" = "Own household", "1" = "Nursing home", "2" = "Hospital", "3" = "Hospice",
                              "4" = "Relative's household", "5" = "Other", "6" = "Not important", "7" = "Palliative care unit"),
  preferred_place_of_death = c("0" = "Own household", "1" = "Nursing home", "2" = "Hospital", "3" = "Hospice",
                               "4" = "Relative's household", "5" = "Other", "6" = "Not important", "7" = "Palliative care unit")
)

eol_display <- eol
for (var in names(codebook)) {
  if (var %in% names(eol_display)) {
    # duplicated labels (Residential_location 3 and 4 -> same label) are collapsed by factor()
    eol_display[[var]] <- factor(eol_display[[var]],
                                 levels = names(codebook[[var]]),
                                 labels = codebook[[var]])
  }
}


## 2.4 Dichotomisation / recoding ==============================================
# All splits below were verified against the raw data. Observed counts are provided.

eol <- eol %>%
  mutate(
    # --- Socio-demographic predictors -----------------------------------------
    gender_dich = factor(dplyr::case_when(gender == 0 ~ 0, gender == 1 ~ 1, TRUE ~ NA_real_),
                         levels = c(0, 1), labels = c("Female", "Male")),

    marital_status_dich = factor(dplyr::case_when(
      marital_status %in% c(1, 2) ~ 1,
      marital_status %in% c(0, 3, 4, 5) ~ 0, TRUE ~ NA_real_),
      levels = c(0, 1), labels = c("Not married", "Married")),

    nationality_dich = factor(dplyr::case_when(
      nationality == 0 ~ 0, nationality %in% 1:11 ~ 1, TRUE ~ NA_real_),
      levels = c(0, 1), labels = c("German", "Non-German")),

    spiritual_orientation_dich = factor(dplyr::case_when(
      religion_belief_worldview %in% 0:7 ~ 0,
      religion_belief_worldview == 8 ~ 1, TRUE ~ NA_real_),
      levels = c(0, 1), labels = c("Affiliated", "None")),

    # Cohabitation without 0 in the data: 1 = alone (n=19), 2-6 = with someone (n=70).
    # 5 people in assisted living / nursing home (Housing_situation 2 or 3) get
    # explicit "Institutional" level rather than NA

    cohabitation_dich = factor(dplyr::case_when(
      Housing_situation %in% c(2, 3) ~ "Institutional",
      Cohabitation == 1 ~ "Alone",
      Cohabitation %in% 2:6 ~ "With someone",
      TRUE ~ NA_character_),
      levels = c("With someone", "Alone", "Institutional")),

    # No support (n=48) vs. any support (n=46).
    nursing_support_dich = factor(dplyr::case_when(
      Nursing_support == 0 ~ 0, Nursing_support %in% 1:4 ~ 1, TRUE ~ NA_real_),
      levels = c(0, 1), labels = c("No", "Yes")),

    # Codes 3 AND 4 are both rural -> rural n=31, non-rural n=63.
    rural_dich = factor(dplyr::case_when(
      Residential_location %in% c(3, 4) ~ 1,
      Residential_location %in% 0:2 ~ 0, TRUE ~ NA_real_),
      levels = c(0, 1), labels = c("No", "Yes")),

    # Fachhochschule / university (n=39) vs. none or vocational (n=55). Code 4 does not occur.
    higher_education_dich = factor(dplyr::case_when(
      Professional_education %in% c(2, 3) ~ 1,
      Professional_education %in% c(0, 1, 4) ~ 0, TRUE ~ NA_real_),
      levels = c(0, 1), labels = c("No", "Yes")),

    # --- Outcomes -------------------------------------------------------------
    # OUTCOMES. Only codes 0/1 occur - there is no "Maybe" category in this dataset.
    advance_directive_dich = factor(dplyr::case_when(
      Advance_directive == 1 ~ 1, Advance_directive == 0 ~ 0, TRUE ~ NA_real_),
      levels = c(0, 1), labels = c("No", "Yes")),      # Yes n=58, No n=36

    power_of_attorney_dich = factor(dplyr::case_when(
      Power_of_attorney == 1 ~ 1, Power_of_attorney == 0 ~ 0, TRUE ~ NA_real_),
      levels = c(0, 1), labels = c("No", "Yes")),      # Yes n=54, No n=40

    # --- End-of-life attitudes / ACP engagement -------------------------------
    # Never (n=20) vs. ever thought about it (n=74).
    end_of_life_thoughts_dich = factor(dplyr::case_when(
      thoughts_about_end_of_life_wishes == 0 ~ 0,
      thoughts_about_end_of_life_wishes %in% 1:3 ~ 1, TRUE ~ NA_real_),
      levels = c(0, 1), labels = c("No", "Yes")),

    # Ever shared (n=77) vs. never (n=17). Code 4 ("I don't know") does not occur in
    # this dataset; if it occurs, should be NA, NOT folded in with "Never".
    sharing_thoughts_dich = factor(dplyr::case_when(
      Sharing_of_thoughts %in% 1:3 ~ 1,
      Sharing_of_thoughts == 0 ~ 0,
      Sharing_of_thoughts == 4 ~ NA_real_, TRUE ~ NA_real_),
      levels = c(0, 1), labels = c("No", "Yes")),

    # --- Clinical -------------------------------------------------------------
    dbs_dich = factor(dplyr::case_when(dbs == 0 ~ 0, dbs == 1 ~ 1, TRUE ~ NA_real_),
                      levels = c(0, 1), labels = c("No", "Yes")),

    # --- Preferred place of care / death --------------------------------------
    # Whether ANY specific place was named (not which one). "Not important" = code 6.
    # care: Yes n=84 / No n=9 / NA n=1  exploratory set only!
    care_preference_dich = factor(dplyr::case_when(
      preferred_place_of_care %in% c(0:5, 7) ~ 1,
      preferred_place_of_care == 6 ~ 0, TRUE ~ NA_real_),
      levels = c(0, 1), labels = c("No", "Yes")),
    
    # death: Yes n=63 / No n=31
    death_preference_dich = factor(dplyr::case_when(
      preferred_place_of_death %in% c(0:5, 7) ~ 1,
      preferred_place_of_death == 6 ~ 0, TRUE ~ NA_real_),
      levels = c(0, 1), labels = c("No", "Yes"))
  )


## 2.5 Combined outcome: BOTH directives =======================================
# Third outcome: BOTH documents present (n=51) vs. not both (n=43).
eol <- eol %>%
  mutate(both_directives_dich = factor(dplyr::case_when(
    advance_directive_dich == "Yes" & power_of_attorney_dich == "Yes" ~ 1,
    !is.na(advance_directive_dich) & !is.na(power_of_attorney_dich) ~ 0,
    TRUE ~ NA_real_),
    levels = c(0, 1), labels = c("No", "Yes")))

cat("\n--- Outcome cross-tabulation (AD x PoA) ---\n")
print(table(AD = eol$advance_directive_dich, PoA = eol$power_of_attorney_dich))
cat("Only", sum(eol$advance_directive_dich != eol$power_of_attorney_dich, na.rm = TRUE),
    "of 94 people have one document without the other.\n")


## 2.6 Rescaled continuous predictors and dichotomised comorbidity =============
# The "_u" suffix means "expressed in the reporting unit set at the top of the script", 
# so an odds ratio for age_u reads "per AGE_UNIT years" rather than per year.

eol <- eol %>%
  mutate(
    age_u   = age       / AGE_UNIT,
    LEDD_u  = LEDD      / LEDD_UNIT,
    UPDRS_u = UPDRS_sum / UPDRS_UNIT,
    PDQ_u   = PDQ_score / PDQ_UNIT,
    bdi_u   = bdi_score / BDI_UNIT,

    # Any comorbidity (n=25) vs. none (n=69). See the rationale in section 1.2.5.
    charlson_dich = factor(dplyr::case_when(
      Charlson_withoutage == 0 ~ 0,
      Charlson_withoutage >= 1 ~ 1,
      TRUE ~ NA_real_),
      levels = c(0, 1), labels = c("None", "Any"))
  )


## 2.7 Noise calibration columns ===============================================
# Pure random columns with no relationship to any outcome; i.e. fixed seed, so that all three outcomes 
# can be compared to each other

set.seed(NOISE_SEED)
eol <- eol %>%
  mutate(
    noise_normal_1  = rnorm(dplyr::n()),
    noise_normal_2  = rnorm(dplyr::n()),
    noise_binary_40 = factor(rbinom(dplyr::n(), 1, 0.40), levels = c(0, 1), labels = c("No", "Yes")),
    noise_binary_15 = factor(rbinom(dplyr::n(), 1, 0.15), levels = c(0, 1), labels = c("No", "Yes"))
  )

cat("\n--- Noise columns vs. outcomes (should all be near zero) ---\n")
for (oc in c("advance_directive_dich", "power_of_attorney_dich", "both_directives_dich")) {
  y <- as.integer(eol[[oc]] == "Yes")
  rr <- vapply(noise_vars, function(v) {
    x <- eol[[v]]
    suppressWarnings(cor(as.numeric(if (is.factor(x)) as.integer(x) else x), y,
                         use = "complete.obs"))
  }, numeric(1))
  cat(sprintf("  %-24s %s\n", oc, paste(sprintf("%s=%+.2f", noise_vars, rr), collapse = "  ")))
}


################################################################################
# 3. ANALYSES                                                               ####
################################################################################

## 3.1 Descriptive tables (Table 1) ============================================

### 3.1.1 Raw codes, labelled --------------------------------------------------
vars_raw <- c("gender", "age", "age_at_diagnosis", "duration", "nationality",
              "marital_status", "religion_belief_worldview", "Housing_situation",
              "Cohabitation", "Nursing_support", "Residential_location",
              "Professional_education", "Advance_directive", "Power_of_attorney",
              "thoughts_about_end_of_life_wishes", "Sharing_of_thoughts",
              "Thoughts_discussed_with", "PDQ_score", "UPDRS_sum", "updrs_III",
              "bdi_score", "MOCA_score", "LEDD",
              "Charlson_withoutage", "Charlson_withage", "dbs",
              "preferred_place_of_care", "preferred_place_of_death")
cat_raw <- setdiff(vars_raw, c("age", "age_at_diagnosis", "duration", "PDQ_score",
                               "UPDRS_sum", "updrs_III", "bdi_score", "MOCA_score",
                               "LEDD", "Charlson_withoutage", "Charlson_withage"))

tab_raw <- tableone::CreateTableOne(vars = vars_raw, data = eol_display, factorVars = cat_raw)
print(tab_raw, showAllLevels = TRUE)
write.csv(print(tab_raw, showAllLevels = TRUE, printToggle = FALSE, quote = FALSE),
          file.path(results_dir, "Table1_raw.csv"))

### 3.1.2 Dichotomised / recoded variables -------------------------------------
vars_dich <- c("gender_dich", "age", "age_at_diagnosis", "duration", "nationality_dich",
               "marital_status_dich", "spiritual_orientation_dich", "cohabitation_dich",
               "nursing_support_dich", "rural_dich", "higher_education_dich",
               "advance_directive_dich", "power_of_attorney_dich", "both_directives_dich",
               "end_of_life_thoughts_dich", "sharing_thoughts_dich",
               "care_preference_dich", "death_preference_dich",
               "PDQ_score", "UPDRS_sum", "updrs_III", "bdi_score", "MOCA_score",
               "LEDD", "Charlson_withoutage", "charlson_dich", "dbs_dich")
cat_dich <- grep("_dich$", vars_dich, value = TRUE)

tab_dich <- tableone::CreateTableOne(vars = vars_dich, data = eol, factorVars = cat_dich)
print(tab_dich, showAllLevels = TRUE)

### 3.1.3 Stratified by outcome ------------------------------------------------
# Stratified by each outcome - this is usually the actual Table 1 of the paper.
for (oc in c("advance_directive_dich", "power_of_attorney_dich", "both_directives_dich")) {
  cat("\n\n#### Stratified by", oc, "####\n")
  tb <- tableone::CreateTableOne(vars = setdiff(vars_dich, oc), strata = oc,
                                 data = eol, factorVars = setdiff(cat_dich, oc))
  print(tb, showAllLevels = TRUE)
  write.csv(print(tb, showAllLevels = TRUE, printToggle = FALSE, quote = FALSE),
            file.path(results_dir, paste0("Table1_by_", oc, ".csv")))
}


## 3.2 Model helper functions ==================================================

### 3.2.1 Design matrix, rank check, events per predictor ----------------------
# Builds a single design matrix ONCE per (outcome, predictor set). Everything
# downstream - the main fit, every bootstrap resample, every prediction - reuses these
# exact columns, so factor contrasts can never drift between fit and predict. This
# replaces v2.2's make_dummy_design()/name-sanitising machinery entirely.
build_design <- function(data, outcome, predictors) {
  cc <- data[stats::complete.cases(data[, c(outcome, predictors)]), , drop = FALSE]
  form_rhs <- as.formula(paste("~", paste(predictors, collapse = " + ")))
  X <- model.matrix(form_rhs, data = cc)[, -1, drop = FALSE]   # drop intercept column
  y <- as.integer(cc[[outcome]] == "Yes")                       # 1 = Yes, 0 = No

  n_dropped <- nrow(data) - nrow(cc)
  if (n_dropped > 0) {
    message(sprintf("build_design(%s): %d of %d rows dropped as incomplete (N = %d used).",
                    outcome, n_dropped, nrow(data), nrow(cc)))
  }
  list(X = X, y = y, terms = colnames(X), n = nrow(cc),
       n_events = sum(y), outcome = outcome, predictors = predictors)
}

# Rank check. If this reports anything, two predictors are (near-)linearly dependent
# and logistf will fail with "LAPACK dpotrf ... not positive definite" rather than
# with a helpful message. Catch it here instead.
check_rank <- function(dsn, label = "") {
  qrX <- qr(cbind(1, dsn$X))
  if (qrX$rank < ncol(dsn$X) + 1) {
    warning(sprintf("%s: design matrix is rank-deficient (rank %d < %d columns). %s",
                    label, qrX$rank, ncol(dsn$X) + 1,
                    "Inspect cor(dsn$X) for a pair at |r| = 1."), call. = FALSE)
    varying <- apply(dsn$X, 2, function(z) length(unique(z)) > 1)
    if (sum(varying) < 2) return(invisible(NULL))
    cm <- cor(dsn$X[, varying, drop = FALSE])
    diag(cm) <- 0
    bad <- which(abs(cm) > 0.98, arr.ind = TRUE)
    if (nrow(bad) > 0) print(unique(t(apply(bad, 1, function(i) sort(colnames(cm)[i])))))
  }
  invisible(NULL)
}

# Events per predictor (parameter, strictly). Peduzzi et al. suggest >= 10.
epp <- function(dsn) {
  min(dsn$n_events, dsn$n - dsn$n_events) / ncol(dsn$X)
}

### 3.2.2 Firth fit, prediction, coefficient table -----------------------------
# Firth-penalized logistic fit on a design built by build_design().
fit_firth <- function(dsn, pl = TRUE, maxit = 500) {
  df <- as.data.frame(dsn$X)
  names(df) <- make.names(dsn$terms, unique = TRUE)  # syntactic names for the formula
  df$.y <- dsn$y
  form <- as.formula(paste(".y ~", paste(names(df)[names(df) != ".y"], collapse = " + ")))
  logistf::logistf(form, data = df, pl = pl,
                   control   = logistf::logistf.control(maxit = maxit),
                   plcontrol = logistf::logistpl.control(maxit = maxit))
}

# Predicted probabilities from a logistf fit for an arbitrary design matrix.
# Computed directly from the coefficients rather than via predict(), so it can never
# be confused by newdata whose factor levels differ from the training data.
predict_firth <- function(fit, X) {
  as.numeric(plogis(cbind(1, X) %*% as.numeric(coef(fit))))
}

# Coefficient table with odds ratios and profile-likelihood (not Wald) CIs.
# Firth's penalized likelihood keeps these finite and valid under separation, which is
# exactly why the original plain glm() version produced coefficients in the hundreds
# with p-values pinned near 0.999.
firth_table <- function(fit, term_labels, alpha = 0.05) {
  data.frame(
    term       = c("(Intercept)", term_labels),
    estimate   = as.numeric(coef(fit)),
    OR         = exp(as.numeric(coef(fit))),
    OR_low     = exp(as.numeric(fit$ci.lower)),
    OR_high    = exp(as.numeric(fit$ci.upper)),
    p_value    = as.numeric(fit$prob),
    row.names  = NULL
  ) %>%
    mutate(significant = p_value < alpha)
}


## 3.3 Optimism-corrected performance (replaces the 80/20 split) ===============
# Efron/Harrell bootstrap optimism correction:
#   1. Fit on the full sample     -> apparent AUC (optimistic, the model saw these data)
#   2. For b in 1..B: resample with replacement, refit,
#        - AUC of the resample model ON the resample   (boot-apparent)
#        - AUC of the resample model ON the ORIGINAL data (boot-test)
#        optimism_b = boot-apparent - boot-test
#   3. corrected AUC = apparent AUC - mean(optimism)
#
# This uses all N for both fitting and validation, and averages over B resamples rather
# than resting on one arbitrary 18-person test partition. It is the standard approach
# for exactly this situation (small N, model development and validation from one sample).
optimism_corrected_auc <- function(dsn, B = 500, seed = 2026) {
  set.seed(seed)

  fit_app <- fit_firth(dsn, pl = FALSE)
  p_app   <- predict_firth(fit_app, dsn$X)
  auc_app <- as.numeric(pROC::auc(pROC::roc(dsn$y, p_app, quiet = TRUE,
                                            levels = c(0, 1), direction = "<")))

  optimism <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    idx <- sample.int(dsn$n, dsn$n, replace = TRUE)
    dsn_b <- list(X = dsn$X[idx, , drop = FALSE], y = dsn$y[idx],
                  terms = dsn$terms, n = dsn$n, n_events = sum(dsn$y[idx]))
    # A resample can be degenerate (all one outcome class, or a dummy column all zero).
    if (length(unique(dsn_b$y)) < 2) next
    res <- tryCatch({
      fit_b <- fit_firth(dsn_b, pl = FALSE)
      a_in  <- as.numeric(pROC::auc(pROC::roc(dsn_b$y, predict_firth(fit_b, dsn_b$X),
                                              quiet = TRUE, levels = c(0, 1), direction = "<")))
      a_out <- as.numeric(pROC::auc(pROC::roc(dsn$y, predict_firth(fit_b, dsn$X),
                                              quiet = TRUE, levels = c(0, 1), direction = "<")))
      a_in - a_out
    }, error = function(e) NA_real_)
    optimism[b] <- res
  }

  n_ok <- sum(!is.na(optimism))
  message(sprintf("optimism_corrected_auc: %d/%d resamples usable.", n_ok, B))
  list(apparent  = auc_app,
       optimism  = mean(optimism, na.rm = TRUE),
       corrected = auc_app - mean(optimism, na.rm = TRUE),
       n_ok      = n_ok,
       roc       = pROC::roc(dsn$y, p_app, quiet = TRUE, levels = c(0, 1), direction = "<"))
}


## 3.4 ElasticNet SELECTION STABILITY (not coefficients) =======================
# WHY THIS IS HERE AND WHAT IT IS NOT
# -----------------------------------
# ElasticNet coefficients are shrunk toward zero by design, carry no standard errors and
# no p-values, and percentile bootstrap intervals around them do not achieve their
# nominal coverage (the estimator is non-differentiable at zero). They must NOT appear in
# a results table as odds ratios. v2.2 reported them as the primary result; that is the
# main thing this version changes.
#
# What ElasticNet IS good for at N=94 with 20 candidate predictors is SELECTION
# STABILITY. Refit it on many bootstrap resamples and count how often each predictor
# survives with a non-zero coefficient. A predictor retained in 85% of resamples is
# telling you something Firth alone cannot; one retained in 12% is noise, and reporting
# that fact is genuinely informative for a paper of this size.
#
# Crucially, lambda (and alpha, if you supply a grid) is RE-TUNED by cross-validation
# inside every resample. v2.2 reused one CV-selected lambda across all 1000 resamples,
# which held the tuning fixed and therefore understated precisely the selection
# uncertainty being measured.
#
# selection_pct is reported; a median coefficient is included for DIRECTION ONLY
# (is the association positive or negative), explicitly not as an effect size.
#
# ------------------------------------------------------------------------------------
# WHY alpha AND THE RESAMPLING SCHEME ARE NOW FIXED RATHER THAN CV-TUNED (v3.2)
# ------------------------------------------------------------------------------------
# The first run of this screen produced a noise floor near 50% - i.e. pure random columns
# survived half the resamples - with cross-validation choosing alpha = 0.25. Two causes,
# both now addressed by the settings above:
#
# 1. LOW ALPHA DOES NOT SELECT. Only the lasso part of the penalty sets coefficients to
#    exactly zero; the ridge part never does. At alpha = 0.25 just a quarter of the
#    penalty is doing selection work, so fitted models stay large and almost everything
#    survives - noise included. CV picked a low alpha because the deviance curve is flat
#    at N=94 with 21 candidates, and a flat curve means CV is nearly indifferent. But CV
#    optimises PREDICTION, and prediction likes ridge. This screen is not predicting
#    anything; it is ranking variables. Optimising for the wrong objective is what
#    produced the uninformative floor, so alpha is now FIXED at 1 (pure lasso), which is
#    what the stability-selection literature assumes.
#
# 2. THE BOOTSTRAP INFLATES SELECTION FREQUENCIES. Resampling n from n WITH replacement
#    duplicates roughly 63% of people, which weakens the effective penalty relative to
#    the information actually present and pushes every selection percentage up.
#    Meinshausen & Buhlmann's stability selection deliberately SUBSAMPLES at n/2 without
#    replacement instead. That is now the default. It is also the reason the published
#    0.6-0.9 thresholds were never comparable to the earlier bootstrap numbers.
#
# Expect the noise floor to drop into roughly the 10-25% range under these settings, and
# the real predictors to spread out enough to be distinguishable. Both toggles are kept
# so the earlier behaviour can be reproduced if you want to show the contrast.
#
# DIAGNOSTIC TO WATCH: mean model size q out of p candidate columns. If there is little
# real signal, a null variable's expected selection rate is roughly q/p - so if the
# average model keeps 12 of 25 columns, noise lands near 48% by arithmetic alone, not
# because anything is broken. The function prints observed noise floor next to this q/p
# prediction so the two can be compared directly.
elasticnet_stability <- function(dsn, B = 500,
                                 alphas    = STABILITY_ALPHA,
                                 scheme    = STABILITY_RESAMPLE,
                                 lambda_rule = STABILITY_LAMBDA,
                                 nfolds = NULL, seed = 2026) {
  set.seed(seed)
  scheme <- match.arg(scheme, c("subsample", "bootstrap"))
  p   <- ncol(dsn$X)
  n_s <- if (scheme == "subsample") floor(dsn$n / 2) else dsn$n

  # Subsamples are half the size, so 10-fold CV can leave too few events per fold.
  # Scale folds to the smaller class rather than hard-coding 10.
  if (is.null(nfolds)) {
    min_class <- min(table(dsn$y)) * (n_s / dsn$n)
    nfolds    <- max(3, min(10, floor(min_class / 3)))
  }

  sel <- matrix(0,          nrow = B, ncol = p, dimnames = list(NULL, dsn$terms))
  est <- matrix(NA_real_,   nrow = B, ncol = p, dimnames = list(NULL, dsn$terms))
  alpha_chosen <- rep(NA_real_, B)
  model_size   <- rep(NA_real_, B)

  for (b in seq_len(B)) {
    idx <- if (scheme == "subsample") sample.int(dsn$n, n_s, replace = FALSE)
           else                        sample.int(dsn$n, n_s, replace = TRUE)
    Xb <- dsn$X[idx, , drop = FALSE]; yb <- dsn$y[idx]
    if (length(unique(yb)) < 2 || min(table(yb)) < nfolds) next
    # Constant columns can appear in a resample (e.g. a rare dummy sampled all-zero);
    # glmnet handles these, but drop them explicitly so the fit is well defined.
    keep <- apply(Xb, 2, function(z) length(unique(z)) > 1)
    if (sum(keep) < 2) next

    # NOTE: the result is returned from tryCatch and assigned OUTSIDE it. Assigning with
    # `<<-` from inside a tryCatch expression would write to the global environment, not
    # to this function's frame, silently leaving sel/est/alpha_chosen empty.
    out_b <- tryCatch({
      # Shared foldid so the CV error is comparable across the alpha grid (a no-op when
      # alphas is a single fixed value, which is now the default).
      foldid <- sample(rep(seq_len(nfolds), length.out = length(yb)))
      best <- NULL
      for (a in alphas) {
        cvf <- glmnet::cv.glmnet(Xb[, keep, drop = FALSE], yb, family = "binomial",
                                 alpha = a, foldid = foldid, type.measure = "deviance")
        lam <- if (lambda_rule == "lambda.min") cvf$lambda.min else cvf$lambda.1se
        cur <- list(alpha = a, lambda = lam,
                    cvm = cvf$cvm[which(cvf$lambda == lam)][1], fit = cvf)
        if (is.null(best) || cur$cvm < best$cvm) best <- cur
      }
      cf <- as.matrix(coef(best$fit, s = best$lambda))[-1, 1]  # drop intercept
      list(cf = cf, alpha = best$alpha)
    }, error = function(e) NULL)

    if (!is.null(out_b)) {
      sel[b, names(out_b$cf)] <- as.numeric(out_b$cf != 0)
      est[b, names(out_b$cf)] <- out_b$cf
      alpha_chosen[b] <- out_b$alpha
      model_size[b]   <- sum(out_b$cf != 0)
    }
  }

  n_ok      <- sum(!is.na(alpha_chosen))
  keep_rows <- !is.na(alpha_chosen)
  q_mean    <- mean(model_size, na.rm = TRUE)

  # ---- Diagnostics ----
  cat(sprintf("\n  [stability diagnostics] scheme = %s (n per resample = %d of %d), %d/%d usable, CV folds = %d\n",
              scheme, n_s, dsn$n, n_ok, B, nfolds))
  if (length(alphas) > 1) {
    cat("  [stability diagnostics] alpha chosen by CV: ")
    print(round(prop.table(table(alpha_chosen[keep_rows])), 2))
    cat(sprintf("  [stability diagnostics] median alpha = %.2f%s\n",
                median(alpha_chosen, na.rm = TRUE),
                if (median(alpha_chosen, na.rm = TRUE) < 0.5)
                  "  <-- WARNING: below 0.5, little is being zeroed; percentages will be inflated" else ""))
  } else {
    cat(sprintf("  [stability diagnostics] alpha FIXED at %.2f (%s)\n",
                alphas[1], if (alphas[1] == 1) "pure lasso" else "elastic net"))
  }
  cat(sprintf("  [stability diagnostics] lambda rule = %s\n", lambda_rule))
  cat(sprintf("  [stability diagnostics] model size q: mean %.1f, median %.0f, range %.0f-%.0f, of p = %d columns\n",
              q_mean, median(model_size, na.rm = TRUE),
              min(model_size, na.rm = TRUE), max(model_size, na.rm = TRUE), p))
  cat(sprintf("  [stability diagnostics] fair-share prediction for a NULL column, q/p = %.0f%%\n",
              100 * q_mean / p))
  cat("  [stability diagnostics] compare that with the observed NOISE FLOOR printed below.\n")
  cat("      close together -> the penalty is not discriminating; nothing is being selected on merit\n")
  cat("      floor well BELOW q/p -> real signal is crowding the noise columns out. Good.\n")

  med <- apply(est[keep_rows, , drop = FALSE], 2, function(z) median(z, na.rm = TRUE))
  med[is.nan(med)] <- NA_real_   # a term never selected in any resample

  out <- data.frame(
    term           = dsn$terms,
    selection_pct  = 100 * colMeans(sel[keep_rows, , drop = FALSE]),
    median_coef    = med,
    # "0" = never selected, so no direction exists
    direction      = ifelse(is.na(med), "0", ifelse(med > 0, "+", "-")),
    row.names = NULL
  ) %>% arrange(desc(selection_pct))

  attr(out, "q_mean")      <- q_mean
  attr(out, "fair_share")  <- 100 * q_mean / p
  attr(out, "scheme")      <- scheme
  attr(out, "median_alpha") <- median(alpha_chosen, na.rm = TRUE)
  out
}


## 3.5 Unadjusted (univariable) Firth screen ===================================
# One predictor at a time, no adjustment. This is the honest way to show every candidate
# variable in a paper this size, and it is what most readers actually want alongside the
# adjusted primary model. Label it clearly as UNADJUSTED - it does not control for
# anything, and with 20 predictors some p < .05 results are expected by chance alone.
univariable_firth <- function(data, outcome, predictors) {
  do.call(rbind, lapply(predictors, function(p) {
    dsn <- build_design(data, outcome, p)
    fit <- fit_firth(dsn)
    tb  <- firth_table(fit, dsn$terms)
    tb  <- tb[tb$term != "(Intercept)", , drop = FALSE]
    tb$predictor <- p
    tb$n <- dsn$n
    tb
  })) %>%
    dplyr::select(predictor, term, OR, OR_low, OR_high, p_value, n) %>%
    arrange(p_value)
}


## 3.6 Per-outcome pipeline: run_outcome() =====================================
run_outcome <- function(outcome, label, data = eol, B = B_BOOT) {
  cat("\n\n============================================================\n")
  cat(label, "\n")
  cat("============================================================\n")

  # ---- PRIMARY: Firth on the pre-specified small set, full N ----
  dsn_p <- build_design(data, outcome, predictors_primary)
  check_rank(dsn_p, paste(label, "primary"))
  cat(sprintf("\nPRIMARY MODEL: N = %d, events = %d, parameters = %d, EPP = %.1f\n",
              dsn_p$n, dsn_p$n_events, ncol(dsn_p$X), epp(dsn_p)))

  fit_p <- fit_firth(dsn_p)
  tab_p <- firth_table(fit_p, dsn_p$terms)
  cat("\n--- Firth-penalized logistic regression (adjusted, profile-likelihood CIs) ---\n")
  print(tab_p %>% mutate(across(c(OR, OR_low, OR_high), ~round(.x, 3)),
                         p_value = round(p_value, 4)))

  perf_p <- optimism_corrected_auc(dsn_p, B = min(B, 500))
  cat(sprintf("\nAUC: apparent %.3f | optimism %.3f | OPTIMISM-CORRECTED %.3f\n",
              perf_p$apparent, perf_p$optimism, perf_p$corrected))

  # ---- UNADJUSTED screen over all candidates ----
  cat("\n--- Unadjusted (univariable) Firth screen, all candidate predictors ---\n")
  uni <- univariable_firth(data, outcome, predictors_exploratory)
  print(uni %>% mutate(across(c(OR, OR_low, OR_high), ~round(.x, 3)),
                       p_value = round(p_value, 4)))

  # ---- EXPLORATORY: ElasticNet selection stability, with noise calibration ----
  dsn_s <- build_design(data, outcome, predictors_stability)
  check_rank(dsn_s, paste(label, "stability"))
  cat(sprintf("\nSTABILITY set: N = %d, events = %d, parameters = %d (incl. %d noise), EPP = %.1f\n",
              dsn_s$n, dsn_s$n_events, ncol(dsn_s$X), length(noise_vars), epp(dsn_s)))
  cat("\n--- ElasticNet selection stability (NOT effect sizes - see section 3.4) ---\n")
  stab <- elasticnet_stability(dsn_s, B = min(B, 500))

  # Identify which rows are the noise columns. build_design() may append a factor level
  # to the column name (e.g. "noise_binary_40Yes"), so match by prefix.
  stab$is_noise <- vapply(stab$term,
                          function(t) any(startsWith(t, noise_vars)), logical(1))
  noise_floor <- if (any(stab$is_noise)) max(stab$selection_pct[stab$is_noise]) else NA_real_

  print(stab %>% mutate(selection_pct = round(selection_pct, 1),
                        median_coef   = round(median_coef, 3)))

  fair <- attr(stab, "fair_share")
  cat(sprintf("\nNOISE FLOOR: the highest-scoring pure-noise column reached %.1f%%.\n", noise_floor))
  cat(sprintf("Fair-share prediction for a null column (q/p): %.1f%%.\n", fair))
  if (!is.na(noise_floor) && !is.na(fair)) {
    if (noise_floor >= 0.8 * fair) {
      cat("=> Noise is at or near its fair share. The penalty is NOT discriminating between\n",
          "   signal and noise at this N. Read the whole table as inconclusive rather than\n",
          "   picking the top rows off it.\n", sep = "")
    } else {
      cat("=> Noise is clearly below its fair share, so real predictors are crowding the noise\n",
          "   columns out. Terms well above the floor are worth discussing; the SIZE of the gap\n",
          "   matters more than the absolute percentage.\n", sep = "")
    }
  }
  cat("Terms at or below the floor are indistinguishable from a column known to contain no signal.\n")

  # ---- SENSITIVITY: same primary model, full predictor set ----
  cat("\n--- SENSITIVITY: Firth on the full exploratory set (EXPLORATORY, low EPP) ---\n")
  # Built WITHOUT the noise columns - they belong in the stability screen only, and this
  # fit is already at EPP < 2 without four extra parameters.
  dsn_e <- build_design(data, outcome, predictors_exploratory)
  check_rank(dsn_e, paste(label, "exploratory"))
  cat(sprintf("(exploratory set: %d parameters, EPP = %.1f)\n", ncol(dsn_e$X), epp(dsn_e)))
  fit_e <- tryCatch(fit_firth(dsn_e), error = function(e) {
    warning("Full-set Firth fit failed: ", conditionMessage(e)); NULL })
  tab_e <- if (!is.null(fit_e)) firth_table(fit_e, dsn_e$terms) else NULL
  if (!is.null(tab_e)) print(tab_e %>% mutate(across(c(OR, OR_low, OR_high), ~round(.x, 3)),
                                              p_value = round(p_value, 4)))

  # ---- Exports ----
  short <- gsub("_dich$", "", outcome)
  write.csv(tab_p, file.path(results_dir, paste0("primaryFirth_", short, ".csv")), row.names = FALSE)
  write.csv(uni,   file.path(results_dir, paste0("univariable_", short, ".csv")),   row.names = FALSE)
  write.csv(stab,  file.path(results_dir, paste0("enetStability_", short, ".csv")), row.names = FALSE)
  if (!is.null(tab_e)) write.csv(tab_e, file.path(results_dir, paste0("fullFirth_", short, ".csv")), row.names = FALSE)

  list(design_primary = dsn_p, fit_primary = fit_p, table_primary = tab_p,
       perf = perf_p, univariable = uni, stability = stab, noise_floor = noise_floor,
       design_full = dsn_e, design_stability = dsn_s, table_full = tab_e, label = label)
}


## 3.7 Correlation structure of the candidate predictors =======================

### 3.7.1 Heatmap function -----------------------------------------------------
# WHY THIS MATTERS FOR READING THE STABILITY SCREEN
# -------------------------------------------------
# ElasticNet cannot tell two correlated predictors apart. When two variables carry
# overlapping signal it retains one in some resamples and the other in the rest, so BOTH
# show a mediocre selection percentage even when the underlying signal is strong. A term
# at 40% can therefore mean either "weak" or "sharing its signal with a correlated
# partner", and the only way to tell which is to look at the correlation structure.
#
# Known pairs in this dataset to check against the stability output:
#   end_of_life_thoughts_dich / sharing_thoughts_dich   agree on 90% of people
#   UPDRS_u / PDQ_u                                     r = 0.53
#   PDQ_u / bdi_u                                       r = 0.57
#   age_u / MOCA_score                                  r = -0.43
#   LEDD_u / duration                                   r = 0.54
#
# Computed on the DUMMY-CODED design matrix, so multi-level factors appear as one row per
# level. Note that a Pearson correlation between two binary columns is the phi
# coefficient and between a binary and a continuous column is the point-biserial
# correlation; both are legitimate here, but a phi coefficient is bounded well below 1
# when the two variables have very different marginal prevalences, so a small value
# between two rare binaries does not prove independence.
plot_correlation_heatmap <- function(data, predictors, label_size = 2.4, label_cut = 0.30) {
  dsn <- build_design(data, "advance_directive_dich", predictors)
  varying <- apply(dsn$X, 2, function(z) length(unique(z)) > 1)
  X <- dsn$X[, varying, drop = FALSE]

  cm <- cor(X, use = "pairwise.complete.obs")

  # Order by hierarchical clustering so blocks of mutually correlated variables sit
  # together and are visible as blocks rather than scattered across the grid.
  ord <- hclust(as.dist(1 - abs(cm)), method = "average")$order
  cm  <- cm[ord, ord]
  vn  <- colnames(cm)

  # Keep the UPPER triangle only (excluding the diagonal) - the matrix is symmetric, so
  # showing both halves doubles the ink for no extra information.
  cm[lower.tri(cm, diag = TRUE)] <- NA

  long <- as.data.frame(as.table(cm))
  names(long) <- c("var1", "var2", "r")
  long <- long[!is.na(long$r), ]
  long$var1 <- factor(long$var1, levels = vn)
  long$var2 <- factor(long$var2, levels = vn)
  long$lab  <- ifelse(abs(long$r) >= label_cut, sprintf("%.2f", long$r), "")

  ggplot(long, aes(x = var1, y = var2, fill = r)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = lab), size = label_size, colour = "grey15") +
    scale_fill_gradient2(low = "firebrick", mid = "white", high = "steelblue4",
                         midpoint = 0, limits = c(-1, 1), name = "r") +
    coord_equal() +
    labs(title = "Correlation among candidate predictors",
         subtitle = sprintf(paste0("Upper triangle, dummy-coded, ordered by hierarchical clustering. ",
                                   "Values shown where |r| >= %.2f.\n",
                                   "Correlated pairs split their ElasticNet selection probability - ",
                                   "read the stability plot with this in mind."), label_cut),
         x = NULL, y = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7),
          panel.grid  = element_blank())
}

### 3.7.2 Draw and save heatmap ------------------------------------------------
fig_corr <- plot_correlation_heatmap(eol, predictors_exploratory)
print(fig_corr)
ggsave(file.path(results_dir, "correlation_heatmap.png"), fig_corr, width = 9, height = 8, dpi = 200)

### 3.7.3 Numeric matrix and highly correlated pairs ---------------------------
# The numeric matrix, for anyone who wants to read exact values rather than squint at tiles.
{
  dsn_c <- build_design(eol, "advance_directive_dich", predictors_exploratory)
  cm_out <- round(cor(dsn_c$X[, apply(dsn_c$X, 2, function(z) length(unique(z)) > 1), drop = FALSE]), 3)
  write.csv(cm_out, file.path(results_dir, "correlation_matrix.csv"))
  # Anything above 0.6 is worth a second look before interpreting the stability screen.
  cm_flag <- cm_out; diag(cm_flag) <- 0
  hi <- which(abs(cm_flag) > 0.6, arr.ind = TRUE)
  if (nrow(hi) > 0) {
    cat("\n--- Predictor pairs with |r| > 0.6 (their selection percentages will be split) ---\n")
    pairs_hi <- unique(t(apply(hi, 1, function(i) sort(colnames(cm_flag)[i]))))
    for (k in seq_len(nrow(pairs_hi))) {
      cat(sprintf("  %-32s %-32s r = %+.2f\n", pairs_hi[k, 1], pairs_hi[k, 2],
                  cm_flag[pairs_hi[k, 1], pairs_hi[k, 2]]))
    }
  } else {
    cat("\nNo predictor pair exceeds |r| = 0.6.\n")
  }
}


## 3.8 Fit all three outcomes ==================================================
# Each call runs the full pipeline from section 3.6 and writes its CSVs.
res_AD   <- run_outcome("advance_directive_dich", "ADVANCE DIRECTIVE")
res_PoA  <- run_outcome("power_of_attorney_dich", "POWER OF ATTORNEY")
res_BOTH <- run_outcome("both_directives_dich",   "BOTH DIRECTIVES (AD and PoA)")


## 3.9 Sensitivity: the other directive as covariate ===========================
# Reported separately and explicitly, never as the primary model. See note 5 at the top:
# with only 10 of 94 people holding one document without the other, this term dominates
# the model and changes what every other coefficient means.
if (SENSITIVITY_INCLUDE_OTHER_DIRECTIVE) {
  cat("\n\n=== SENSITIVITY: primary model + the other directive as covariate ===\n")
  for (cfg in list(list(o = "advance_directive_dich", x = "power_of_attorney_dich", l = "AD | + PoA"),
                   list(o = "power_of_attorney_dich", x = "advance_directive_dich", l = "PoA | + AD"))) {
    dsn <- build_design(eol, cfg$o, c(predictors_primary, cfg$x))
    fit <- fit_firth(dsn)
    cat("\n---", cfg$l, "---\n")
    print(firth_table(fit, dsn$terms) %>%
            mutate(across(c(OR, OR_low, OR_high), ~round(.x, 3)), p_value = round(p_value, 4)))
    write.csv(firth_table(fit, dsn$terms),
              file.path(results_dir, paste0("sensitivity_", gsub("[^A-Za-z]", "", cfg$l), ".csv")),
              row.names = FALSE)
  }
}


################################################################################
# 4. REPORT                                                                 ####
################################################################################


## 4.1 Figure functions ========================================================

### 4.1.1 Forest plot (primary Firth model) ------------------------------------
# Forest plot of the primary Firth model (odds ratios, profile-likelihood CIs).
forest_plot <- function(res) {
  d <- res$table_primary %>%
    filter(term != "(Intercept)") %>%
    mutate(term = factor(term, levels = rev(term)))
  ggplot(d, aes(x = OR, y = term)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = OR_low, xmax = OR_high), height = 0.15) +
    geom_point(size = 2.6, colour = "steelblue4") +
    scale_x_log10() +
    labs(title = paste0("Adjusted odds ratios - ", res$label),
         subtitle = sprintf("Firth-penalized logistic regression, N = %d, %d events, EPP = %.1f",
                            res$design_primary$n, res$design_primary$n_events, epp(res$design_primary)),
         x = "Odds ratio (log scale, profile-likelihood 95% CI)", y = NULL) +
    theme_minimal()
}

### 4.1.2 Stability plot (ElasticNet selection frequencies) --------------------
# Selection-frequency plot for the ElasticNet stability screen.
# The noise columns are always shown regardless of where they rank, and the noise floor
# (highest-scoring noise column) is drawn as a vertical reference line. Read the plot as:
# anything to the LEFT of that line performed no better than a variable known to contain
# no signal. This replaces an arbitrary borrowed threshold with one calibrated on these
# data, at this N, with this correlation structure.
stability_plot <- function(res, top_n = 20) {
  d <- res$stability
  d <- rbind(head(d[!d$is_noise, ], top_n), d[d$is_noise, ])
  d <- d[order(d$selection_pct), ]
  d$term  <- factor(d$term, levels = d$term)
  d$fill_grp <- ifelse(d$is_noise, "noise", d$direction)

  p <- ggplot(d, aes(x = selection_pct, y = term, fill = fill_grp)) +
    geom_col() +
    scale_fill_manual(
      values = c("+" = "steelblue4", "-" = "firebrick", "0" = "grey80", "noise" = "grey45"),
      breaks = c("+", "-", "0", "noise"),
      labels = c("Positive", "Negative", "Never selected", "NOISE (reference)"),
      name = NULL) +
    labs(title = paste0("ElasticNet selection stability - ", res$label),
         x = "Selected in % of bootstrap resamples", y = NULL) +
    theme_minimal()

  if (!is.na(res$noise_floor)) {
    fair <- attr(res$stability, "fair_share")
    p <- p +
      geom_vline(xintercept = res$noise_floor, linetype = "dashed",
                 colour = "grey20", linewidth = 0.7)
    if (!is.null(fair) && !is.na(fair)) {
      p <- p + geom_vline(xintercept = fair, linetype = "dotted",
                          colour = "grey55", linewidth = 0.6)
    }
    p <- p + labs(subtitle = sprintf(
      paste0("Dashed = noise floor (%.0f%%): best pure-noise column. Dotted = fair share q/p (%.0f%%).\n",
             "Bars left of the dashed line are indistinguishable from noise. ",
             "Colour shows direction only - NOT effect size.\n",
             "Scheme: %s, alpha = %s."),
      res$noise_floor, if (is.null(fair)) NA_real_ else fair,
      attr(res$stability, "scheme"),
      paste(STABILITY_ALPHA, collapse = "/")))
  }
  p
}


## 4.2 Save figures for each outcome ===========================================
for (r in list(res_AD, res_PoA, res_BOTH)) {
  short <- gsub("[^A-Za-z]", "", r$label)
  ggsave(file.path(results_dir, paste0("forest_", short, ".png")), forest_plot(r), width = 7, height = 4.5)
  ggsave(file.path(results_dir, paste0("stability_", short, ".png")), stability_plot(r), width = 7.5, height = 6)
}


## 4.3 One-page performance summary across outcomes ============================
perf_summary <- data.frame(
  outcome        = c(res_AD$label, res_PoA$label, res_BOTH$label),
  N              = c(res_AD$design_primary$n, res_PoA$design_primary$n, res_BOTH$design_primary$n),
  events         = c(res_AD$design_primary$n_events, res_PoA$design_primary$n_events, res_BOTH$design_primary$n_events),
  EPP            = round(c(epp(res_AD$design_primary), epp(res_PoA$design_primary), epp(res_BOTH$design_primary)), 1),
  AUC_apparent   = round(c(res_AD$perf$apparent,  res_PoA$perf$apparent,  res_BOTH$perf$apparent), 3),
  AUC_corrected  = round(c(res_AD$perf$corrected, res_PoA$perf$corrected, res_BOTH$perf$corrected), 3)
)
cat("\n\n--- Primary model performance (optimism-corrected) ---\n")
print(perf_summary)
write.csv(perf_summary, file.path(results_dir, "performance_summary.csv"), row.names = FALSE)

cat("\nAll output written to:", results_dir, "\n")


## 4.4 Interpreting this output - a short guide ================================
# primaryFirth_*.csv     THE result table for the paper. Adjusted odds ratios with
#                        profile-likelihood 95% CIs and penalized-likelihood-ratio
#                        p-values, from the pre-specified 5-predictor model on all N.
#
# univariable_*.csv      Unadjusted associations for every candidate variable. Report
#                        alongside the primary table, labelled UNADJUSTED. Do not
#                        interpret individual p-values here as confirmatory: with 20
#                        predictors, one or two below .05 are expected by chance.
#
# enetStability_*.csv    selection_pct = how often the term survived across bootstrap
#                        resamples with alpha/lambda re-tuned each time. Read as
#                        evidence of robustness, NOT as effect size; median_coef gives
#                        direction only.
#                        Do NOT apply a fixed cutoff. Compare each term against the
#                        NOISE FLOOR printed under the table (the best-scoring pure-noise
#                        column, is_noise == TRUE). A term below the floor performed no
#                        better than a variable known to contain nothing. A term well
#                        above it is worth discussing. In between is genuinely ambiguous
#                        and should be described as such rather than dichotomised.
#                        ALSO CHECK the median alpha printed by elasticnet_stability():
#                        if cross-validation keeps choosing a low alpha (near-ridge),
#                        few coefficients are ever set exactly to zero and ALL
#                        percentages inflate toward 100, making the whole screen
#                        uninformative. Median alpha of 0.5 or above means the numbers
#                        are interpretable.
#
# correlation_heatmap    Read BEFORE the stability table. Correlated predictors split
# correlation_matrix.csv their selection probability, so two variables each at 40% may
#                        be one strong signal rather than two weak ones. Pairs above
#                        |r| = 0.6 are listed in the console.
#
# fullFirth_*.csv        The 20-predictor model. EXPLORATORY. Expect wide CIs. Useful
#                        to check that the primary estimates do not flip sign, not as
#                        a source of p-values.
#
# performance_summary    AUC_apparent is optimistic (the model saw these data).
#                        AUC_corrected is the honest estimate. If corrected AUC sits
#                        near 0.5-0.6, the defensible conclusion is that these
#                        variables do not predict advance care planning well in this
#                        sample - which, with N=94, is a legitimate and publishable
#                        finding rather than a failure.