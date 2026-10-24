# ===========================================================
# Project: Parkinson's EOLC Planning - predictors of advance care planning
# Script:  AnalyseVorsorge_v3_4.R
# Authors: Anna and David Pedrosa
# Version: 3.4
#
# Version history:       CHANGELOG.md
# Reading the output:    interpretation.R (separate, not part of the analysis)
# ===========================================================


# ===========================================================
# 0. USER SETTINGS - edit these, then source the whole file
# ===========================================================

PPR_INSTALL_PACKAGES <- TRUE
B_BOOT <- 1000 # Number of bootstrap resamples
# Set TRUE to additionally fit the AD and PoA models WITH the other directive as a
# covariate (sensitivity analysis, section 12). Does not affect primary output.
SENSITIVITY_INCLUDE_OTHER_DIRECTIVE <- TRUE

# --- UNITS FOR CONTINUOUS PREDICTORS ---------------------------------------------
#   age        -> 5 years
#   LEDD       -> 100 mg   (conventional unit in PD papers)
#   UPDRS_sum  -> 5 points (MCID)
#   PDQ_score  -> 5 points (MCID)
#   bdi_score  -> 5 points (MCID)
#   MOCA_score -> 1 point  (unscaled)
#   duration   -> 1 year   (unscaled)
#
AGE_UNIT   <- 5     # years per unit of the age term
LEDD_UNIT  <- 100   # mg per unit of the LEDD term
UPDRS_UNIT <- 5    # UPDRS points per unit
PDQ_UNIT   <- 5    # PDQ points per unit
BDI_UNIT   <- 5     # BDI points per unit

# --- STABILITY SCREEN SETTINGS ----------------------------------------------------
# STABILITY_ALPHA: the elastic net mixing parameter. Only the LASSO component (alpha
#   towards 1) sets coefficients to exactly zero; the ridge component (alpha towards 0)
#   never does. Fixed at 1 = pure lasso, because this screen exists to RANK variables,
#   not to predict. Letting cross-validation choose alpha optimises prediction, and
#   prediction prefers ridge. Supply a vector here (e.g. c(0.5, 0.75, 1)) to go back to
#   CV-tuned alpha; the function will then report the chosen-alpha distribution.
STABILITY_ALPHA <- 1
STABILITY_RESAMPLE <- "subsample" # "subsample" draws n/2 WITHOUT replacement
STABILITY_LAMBDA <- "lambda.min"# "lambda.1se"
# Smallest lambda tried, as a fraction of the largest. glmnet's default (1e-4) reaches
# almost unpenalised fits, which do not converge on half-samples (separation)
# and are never the CV optimum anyway. 0.01 keeps the path in the range that matters.
STABILITY_LMR <- 0.01

# --- NOISE CALIBRATION VARIABLES --------------------------------------------------
# This generates pure-noise columns (two continous and two binary) with no relationship 
# to anything, carries them through the ElasticNet stability screen, and reports what 
# selection percentage noise achieves in this dataset, at this N, with this correlation 
# structure.
#
NOISE_SEED <- 4711
noise_vars <- c("noise_normal_1", "noise_normal_2", "noise_binary_40", "noise_binary_15")

# --- THE PRE-SPECIFIED PRIMARY PREDICTOR SET -------------------------------------
# This is the single most important decision in the script and it should be made on
# clinical/theoretical grounds BEFORE looking at any model output. Change it now, once,
# and then leave it alone.
#
# Rationale for each:
#   age_u                      - advance care planning rises with increasing age 
#   gender_dich                - consistently reported sex difference
#   marital_status_dich        - a spouse is the usual proxy decision-maker
#   higher_education_dich      - health literacy / access to legal advice.
#   death_preference_dich      - a stated preferred place of death; marker of having
#                                engaged with end-of-life wishes (replaces
#                                end_of_life_thoughts_dich / care_preference_dich).
#   LEDD_u                     - disease severity proxy
#   charlson_dich              - comorbidity. Dichotomised rather than continuous
#                                because most participants score 0

predictors_primary <- c(
  "age_u",
  "gender_dich",
  "marital_status_dich",
  "higher_education_dich",
  "death_preference_dich",
  # "end_of_life_thoughts_dich",
  "LEDD_u",
  "charlson_dich"
)

# --- THE EXPLORATORY SET (sensitivity Firth fit) ----------------------------------
# Deliberately NOT used for the primary p-values. housing_dich is excluded: it is
# numerically identical to cohabitation_dich == "Institutional" in this dataset, which
# makes the design matrix rank-deficient.
# Charlson_withoutage rather than Charlson_withage, since the latter bakes in age and
# would double-count it alongside the standalone age term.
# updrs_III is NOT included alongside UPDRS_u: they are highly correlated.
predictors_exploratory <- c(
  "age_u", "gender_dich", "nationality_dich", "marital_status_dich",
  "spiritual_orientation_dich", "cohabitation_dich",
  "nursing_support_dich", "rural_dich", "higher_education_dich",
  "end_of_life_thoughts_dich", "sharing_thoughts_dich",
  "care_preference_dich", "death_preference_dich",
  "duration", "dbs_dich", "LEDD_u", "UPDRS_u", "PDQ_u", "MOCA_score",
  "bdi_u", "charlson_dich"
)

# --- THE STABILITY-SCREEN SET -----------------------------------------------------
predictors_stability <- c(predictors_exploratory, noise_vars) # The exploratory set PLUS the noise columns. 


# ===========================================================
# 1. Packages
# ===========================================================
required_pkgs <- c("dplyr", "tidyr", "ggplot2", "readxl", "tableone",
                   "logistf", "glmnet", "pROC")

install_missing <- isTRUE(PPR_INSTALL_PACKAGES)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  if (!install_missing) {
    stop(paste0("Missing package(s): ", paste(missing_pkgs, collapse = ", "),
                "\nInstall them, or set PPR_INSTALL_PACKAGES=TRUE."), call. = FALSE)
  }
  install.packages(missing_pkgs, dependencies = TRUE)
}
suppressPackageStartupMessages(for (p in required_pkgs) library(p, character.only = TRUE))

# ===========================================================
# 2. Directories and data
# ===========================================================
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
results_dir <- file.path(wdir, "results_final")
dir.create(results_dir, showWarnings = FALSE)

eol <- readxl::read_excel(file.path(data_dir, "Datenmatrix_end-of-life_final_version11.10.23.xlsx"))

names(eol) <- trimws(names(eol))
eol <- eol %>% dplyr::rename(Thoughts_discussed_with = Thoughts_dicussed_with)

stopifnot(nrow(eol) == 94)


# ===========================================================
# 3. PDQ-39 subscales
# ===========================================================
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
message(sprintf("PDQ_score: mean = %.2f, SD = %.2f",
                mean(eol$PDQ_score), sd(eol$PDQ_score)))


# ===========================================================
# 4. Codebook and labelled display copy
# ===========================================================
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

# ===========================================================
# 5. Dichotomisation / recoding
# ===========================================================
eol <- eol %>%
  mutate(
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

    cohabitation_dich = factor(dplyr::case_when(
      Housing_situation %in% c(2, 3) ~ "Institutional",
      Cohabitation == 1 ~ "Alone",
      Cohabitation %in% 2:6 ~ "With someone",
      TRUE ~ NA_character_),
      levels = c("With someone", "Alone", "Institutional")),

    nursing_support_dich = factor(dplyr::case_when(
      Nursing_support == 0 ~ 0, Nursing_support %in% 1:4 ~ 1, TRUE ~ NA_real_),
      levels = c(0, 1), labels = c("No", "Yes")),

    rural_dich = factor(dplyr::case_when(
      Residential_location %in% c(3, 4) ~ 1,
      Residential_location %in% 0:2 ~ 0, TRUE ~ NA_real_),
      levels = c(0, 1), labels = c("No", "Yes")),

    # Fachhochschule / university vs. none, vocational or other
    higher_education_dich = factor(dplyr::case_when(
      Professional_education %in% c(2, 3) ~ 1,
      Professional_education %in% c(0, 1, 4) ~ 0, TRUE ~ NA_real_),
      levels = c(0, 1), labels = c("No", "Yes")),

    # OUTCOMES.
    advance_directive_dich = factor(dplyr::case_when(
      Advance_directive == 1 ~ 1, Advance_directive == 0 ~ 0, TRUE ~ NA_real_),
      levels = c(0, 1), labels = c("No", "Yes")),

    power_of_attorney_dich = factor(dplyr::case_when(
      Power_of_attorney == 1 ~ 1, Power_of_attorney == 0 ~ 0, TRUE ~ NA_real_),
      levels = c(0, 1), labels = c("No", "Yes")),

    # Never vs. ever thought about it
    end_of_life_thoughts_dich = factor(dplyr::case_when(
      thoughts_about_end_of_life_wishes == 0 ~ 0,
      thoughts_about_end_of_life_wishes %in% 1:3 ~ 1, TRUE ~ NA_real_),
      levels = c(0, 1), labels = c("No", "Yes")),

    # Ever shared vs. never. Code 4 ("I don't know") is set to NA, not folded into "Never".
    sharing_thoughts_dich = factor(dplyr::case_when(
      Sharing_of_thoughts %in% 1:3 ~ 1,
      Sharing_of_thoughts == 0 ~ 0,
      Sharing_of_thoughts == 4 ~ NA_real_, TRUE ~ NA_real_),
      levels = c(0, 1), labels = c("No", "Yes")),

    dbs_dich = factor(dplyr::case_when(dbs == 0 ~ 0, dbs == 1 ~ 1, TRUE ~ NA_real_),
                      levels = c(0, 1), labels = c("No", "Yes")),

    # Whether ANY specific place was named (not which one). "Not important" = code 6.
    care_preference_dich = factor(dplyr::case_when(
      preferred_place_of_care %in% c(0:5, 7) ~ 1,
      preferred_place_of_care == 6 ~ 0, TRUE ~ NA_real_),
      levels = c(0, 1), labels = c("No", "Yes")),
    death_preference_dich = factor(dplyr::case_when(
      preferred_place_of_death %in% c(0:5, 7) ~ 1,
      preferred_place_of_death == 6 ~ 0, TRUE ~ NA_real_),
      levels = c(0, 1), labels = c("No", "Yes"))
  )

# Third outcome: BOTH documents present vs. not both.
eol <- eol %>%
  mutate(both_directives_dich = factor(dplyr::case_when(
    advance_directive_dich == "Yes" & power_of_attorney_dich == "Yes" ~ 1,
    !is.na(advance_directive_dich) & !is.na(power_of_attorney_dich) ~ 0,
    TRUE ~ NA_real_),
    levels = c(0, 1), labels = c("No", "Yes")))

cat("\n--- Outcome cross-tabulation (AD x PoA) ---\n")
print(table(AD = eol$advance_directive_dich, PoA = eol$power_of_attorney_dich))
cat("Only", sum(eol$advance_directive_dich != eol$power_of_attorney_dich, na.rm = TRUE),
    "of", nrow(eol), "people have one document without the other.\n")


# ===========================================================
# 5b. Rescaled continuous predictors and dichotomised comorbidity
# ===========================================================
# See the UNITS block at the beginningin section 0. The "_u" suffix means "expressed in the reporting
# unit set at the top of the script"

eol <- eol %>%
  mutate(
    age_u   = age       / AGE_UNIT,
    LEDD_u  = LEDD      / LEDD_UNIT,
    UPDRS_u = UPDRS_sum / UPDRS_UNIT,
    PDQ_u   = PDQ_score / PDQ_UNIT,
    bdi_u   = bdi_score / BDI_UNIT,

    # Any comorbidity vs. none
    charlson_dich = factor(dplyr::case_when(
      Charlson_withoutage == 0 ~ 0,
      Charlson_withoutage >= 1 ~ 1,
      TRUE ~ NA_real_),
      levels = c(0, 1), labels = c("None", "Any"))
  )

message(sprintf("LEDD: n = %d, missing = %d, zeros = %d, median = %.0f, max = %.0f",
                sum(!is.na(eol$LEDD)), sum(is.na(eol$LEDD)),
                sum(eol$LEDD == 0, na.rm = TRUE),
                median(eol$LEDD, na.rm = TRUE), max(eol$LEDD, na.rm = TRUE)))

cat("\n--- LEDD by DBS status (check this is what you expect) ---\n")
print(eol %>% group_by(dbs_dich) %>%
        summarise(n = n(), mean_LEDD = round(mean(LEDD, na.rm = TRUE)),
                  sd_LEDD = round(sd(LEDD, na.rm = TRUE)),
                  mean_duration = round(mean(duration), 1), .groups = "drop"))


# ===========================================================
# 5c. Noise calibration columns
# ===========================================================
# Pure random columns with no relationship to any outcome. The fixed seed matters: all 
# three outcomes must see the SAME noise columns, otherwise their stability tables are
# are calibrated against different reference variables and cannot be compared

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


# ===========================================================
# 6. Descriptive tables
# ===========================================================
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

tab_raw <- tableone::CreateTableOne(vars = vars_raw, data = as.data.frame(eol_display), factorVars = cat_raw)
print(tab_raw, showAllLevels = TRUE)
write.csv(print(tab_raw, showAllLevels = TRUE, printToggle = FALSE, quote = FALSE),
          file.path(results_dir, "Table1_raw.csv"))

vars_dich <- c("gender_dich", "age", "age_at_diagnosis", "duration", "nationality_dich",
               "marital_status_dich", "spiritual_orientation_dich", "cohabitation_dich",
               "nursing_support_dich", "rural_dich", "higher_education_dich",
               "advance_directive_dich", "power_of_attorney_dich", "both_directives_dich",
               "end_of_life_thoughts_dich", "sharing_thoughts_dich",
               "care_preference_dich", "death_preference_dich",
               "PDQ_score", "UPDRS_sum", "updrs_III", "bdi_score", "MOCA_score",
               "LEDD", "Charlson_withoutage", "charlson_dich", "dbs_dich")
cat_dich <- grep("_dich$", vars_dich, value = TRUE)

tab_dich <- tableone::CreateTableOne(vars = vars_dich, data = as.data.frame(eol), factorVars = cat_dich)
print(tab_dich, showAllLevels = TRUE)

# Stratified by each outcome - this is usually the actual Table 1 of the paper.
# as.data.frame(): with a tibble, tableone silently returns NA p-values for continuous variables.
for (oc in c("advance_directive_dich", "power_of_attorney_dich", "both_directives_dich")) {
  cat("\n\n#### Stratified by", oc, "####\n")
  tb <- tableone::CreateTableOne(vars = setdiff(vars_dich, oc), strata = oc,
                                 data = as.data.frame(eol), factorVars = setdiff(cat_dich, oc))
  print(tb, showAllLevels = TRUE)
  write.csv(print(tb, showAllLevels = TRUE, printToggle = FALSE, quote = FALSE),
            file.path(results_dir, paste0("Table1_by_", oc, ".csv")))
}


# ===========================================================
# 7. Core helpers
# ===========================================================

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


# ===========================================================
# 8. Optimism-corrected performance (replaces the 80/20 split)
# ===========================================================
# Efron/Harrell bootstrap optimism correction [ 10.1002/(SICI)1097-0258(19960229)15:4<361::AID-SIM168>3.0.CO;2-4]: 
#   1. Fit on the full sample     -> apparent AUC (optimistic, the model saw these data)
#   2. For b in 1..B: resample with replacement, refit,
#        - AUC of the resample model ON the resample   (boot-apparent)
#        - AUC of the resample model ON the ORIGINAL data (boot-test)
#        optimism_b = boot-apparent - boot-test
#   3. corrected AUC = apparent AUC - mean(optimism)

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


# ===========================================================
# 9. ElasticNet SELECTION STABILITY (not coefficients)
# ===========================================================
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
  at_edge      <- rep(NA, B)

  for (b in seq_len(B)) {
    idx <- if (scheme == "subsample") sample.int(dsn$n, n_s, replace = FALSE)
           else                        sample.int(dsn$n, n_s, replace = TRUE)
    Xb <- dsn$X[idx, , drop = FALSE]; yb <- dsn$y[idx]
    if (length(unique(yb)) < 2 || min(table(yb)) < nfolds) next
    keep <- apply(Xb, 2, function(z) length(unique(z)) > 1)
    if (sum(keep) < 2) next

    out_b <- tryCatch({
      # Shared foldid so the CV error is comparable across the alpha grid (a no-op when
      # alphas is a single fixed value, which is now the default).
      foldid <- sample(rep(seq_len(nfolds), length.out = length(yb)))
      best <- NULL
      for (a in alphas) {
        cvf <- glmnet::cv.glmnet(Xb[, keep, drop = FALSE], yb, family = "binomial",
                                 alpha = a, foldid = foldid, type.measure = "deviance",
                                 lambda.min.ratio = STABILITY_LMR)
        lam <- if (lambda_rule == "lambda.min") cvf$lambda.min else cvf$lambda.1se
        cur <- list(alpha = a, lambda = lam, at_edge = lam == min(cvf$lambda),
                    cvm = cvf$cvm[which(cvf$lambda == lam)][1], fit = cvf)
        if (is.null(best) || cur$cvm < best$cvm) best <- cur
      }
      cf <- as.matrix(coef(best$fit, s = best$lambda))[-1, 1]  # drop intercept
      list(cf = cf, alpha = best$alpha, at_edge = best$at_edge)
    }, error = function(e) NULL)

    if (!is.null(out_b)) {
      sel[b, names(out_b$cf)] <- as.numeric(out_b$cf != 0)
      est[b, names(out_b$cf)] <- out_b$cf
      alpha_chosen[b] <- out_b$alpha
      model_size[b]   <- sum(out_b$cf != 0)
      at_edge[b]      <- out_b$at_edge
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
  cat(sprintf("  [stability diagnostics] lambda rule = %s, lambda.min.ratio = %g\n", lambda_rule, STABILITY_LMR))
  cat(sprintf("  [stability diagnostics] chosen lambda at the small end of the path in %.0f%% of resamples%s\n",
              100 * mean(at_edge, na.rm = TRUE),
              if (mean(at_edge, na.rm = TRUE) > 0.10) "  <-- WARNING: lower STABILITY_LMR" else ""))
  cat(sprintf("  [stability diagnostics] model size q: mean %.1f, median %.0f, range %.0f-%.0f, of p = %d columns\n",
              q_mean, median(model_size, na.rm = TRUE),
              min(model_size, na.rm = TRUE), max(model_size, na.rm = TRUE), p))
  cat(sprintf("  [stability diagnostics] fair-share prediction for a NULL column, q/p = %.0f%%\n",
              100 * q_mean / p))

  # Median over the resamples in which the term was SELECTED. Including the zeros would
  # give a median of exactly 0 for any term selected in <50% of resamples, which the
  # direction rule below would then mislabel as negative.
  med <- apply(est[keep_rows, , drop = FALSE], 2, function(z) {
    z <- z[!is.na(z) & z != 0]
    if (length(z) == 0) NA_real_ else median(z)
  })

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


# ===========================================================
# 10. Unadjusted (univariable) Firth screen
# ===========================================================
# One predictor at a time, no adjustment. This is UNADJUSTED - it does not control for
# anything, so that with 20 predictors some p < .05 are expected by chance.
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


# ===========================================================
# 11. Run everything for one outcome (MAIN analysis)
# ===========================================================
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
  cat("\n--- ElasticNet selection stability (NOT effect sizes - see section 9) ---\n")
  stab <- elasticnet_stability(dsn_s, B = min(B, 500))

  stab$is_noise <- vapply(stab$term,
                          function(t) any(startsWith(t, noise_vars)), logical(1))
  noise_floor <- if (any(stab$is_noise)) max(stab$selection_pct[stab$is_noise]) else NA_real_

  print(stab %>% mutate(selection_pct = round(selection_pct, 1),
                        median_coef   = round(median_coef, 3)))

  fair <- attr(stab, "fair_share")
  cat(sprintf("\nNOISE FLOOR: the highest-scoring pure-noise column reached %.1f%%.\n", noise_floor))
  cat(sprintf("Fair-share prediction for a null column (q/p): %.1f%%.\n", fair))

  # ---- SENSITIVITY: same primary model, full predictor set ----
  cat("\n--- SENSITIVITY: Firth on the full exploratory set (EXPLORATORY, low EPP) ---\n")
  # Built WITHOUT the noise columns
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

# ===========================================================
# 11b. Correlation heatmap of the candidate predictors
# ===========================================================
# Computed on the dummy-coded design matrix (one row per factor level); Pearson r between
# two binary columns is the phi coefficient, between binary and continuous point-biserial.
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

fig_corr <- plot_correlation_heatmap(eol, predictors_exploratory)
print(fig_corr)
ggsave(file.path(results_dir, "correlation_heatmap.png"), fig_corr, width = 9, height = 8, dpi = 200)

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


res_AD   <- run_outcome("advance_directive_dich", "ADVANCE DIRECTIVE")
res_PoA  <- run_outcome("power_of_attorney_dich", "POWER OF ATTORNEY")
res_BOTH <- run_outcome("both_directives_dich",   "BOTH DIRECTIVES (AD and PoA)")


# ===========================================================
# 12. Sensitivity: the other directive as covariate
# ===========================================================
# Reported separately, never as the primary model: the two documents almost always co-occur,
# so this term dominates the model and changes what every other coefficient means.
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


# ===========================================================
# 13. Figures
# ===========================================================

# Forest plot of the primary Firth model (odds ratios, profile-likelihood CIs).
forest_plot <- function(res) {
  d <- res$table_primary %>%
    filter(term != "(Intercept)") %>%
    mutate(term = factor(term, levels = rev(term)))
  ggplot(d, aes(x = OR, y = term)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(xmin = OR_low, xmax = OR_high), width = 0.15, orientation = "y") +
    geom_point(size = 2.6, colour = "steelblue4") +
    scale_x_log10() +
    labs(title = paste0("Adjusted odds ratios - ", res$label),
         subtitle = sprintf("Firth-penalized logistic regression, N = %d, %d events, EPP = %.1f",
                            res$design_primary$n, res$design_primary$n_events, epp(res$design_primary)),
         x = "Odds ratio (log scale, profile-likelihood 95% CI)", y = NULL) +
    theme_minimal()
}

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
         x = "Selected in % of resamples", y = NULL) +
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
             "Colour shows direction only, not effect size.\n",
             "Scheme: %s, alpha = %s."),
      res$noise_floor, if (is.null(fair)) NA_real_ else fair,
      attr(res$stability, "scheme"),
      paste(STABILITY_ALPHA, collapse = "/")))
  }
  p
}

for (r in list(res_AD, res_PoA, res_BOTH)) {
  short <- gsub("[^A-Za-z]", "", r$label)
  ggsave(file.path(results_dir, paste0("forest_", short, ".png")), forest_plot(r), width = 7, height = 4.5)
  ggsave(file.path(results_dir, paste0("stability_", short, ".png")), stability_plot(r), width = 7.5, height = 6)
}


# ===========================================================
# 14. One-page summary across outcomes
# ===========================================================
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
