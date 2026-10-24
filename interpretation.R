# ===========================================================
# interpretation.R  -  writing aid, NOT part of the repository
# ===========================================================
# Reads the CSV output of AnalyseVorsorge_v3_4.R and walks through it in the order the
# results section would present it. The numbers in the comments are from the run of
# 2026-10-24 (B = 500 resamples); the code below recomputes everything from the CSVs, so
# rerun it after any new analysis run and check the comments still match.
#
# Order:
#   1. Sample and outcomes
#   2. How strongly AD and PoA overlap (sensitivity models)
#   3. Unadjusted associations
#   4. Primary model and discrimination
#   5. Exploratory 21-predictor model
#   6. Selection stability screen and predictor correlations
#   7. Synthesis and open points
# ===========================================================

library(dplyr)

res_dir <- "/Users/annapedrosa/Documents/dokumente/Arbeit_Anna/Aktuelle Projekte/StudiePoD_PD/results_final"
rd <- function(f) read.csv(file.path(res_dir, f), check.names = FALSE)
outcomes <- c(AD = "advance_directive", PoA = "power_of_attorney", BOTH = "both_directives")

show_or <- function(d) {
  d %>% filter(term != "(Intercept)") %>%
    mutate(across(c(OR, OR_low, OR_high), ~ round(.x, 2)), p_value = signif(p_value, 2)) %>%
    select(term, OR, OR_low, OR_high, p_value)
}


# ===========================================================
# 1. Sample and outcomes  (Table1_raw.csv)
# ===========================================================
# N = 94; 66 men (70.2%); age 65.1 (SD 8.7); disease duration 10.5 (5.9) years;
# DBS 52 (55.3%); LEDD 728 (487) mg; PDQ-39 SI 25.5 (13.7); MoCA 25.9 (3.2).
# Advance directive 58 (61.7%), power of attorney 54 (57.4%), both 51 (54.3%).
# 20 (21.3%) had never thought about end-of-life wishes; 31 (33.0%) said the place of
# death was not important, 41 (43.6%) preferred their own home.
print(rd("Table1_raw.csv"))


# ===========================================================
# 2. Overlap of the two documents  (sensitivity_*.csv)
# ===========================================================
# Cross-tab: neither 33, AD only 7, PoA only 3, both 51 -> only 10 of 94 discordant.
# Adding the other document to the primary model:
#   AD  | + PoA : OR 63.7 (16.4-420),  p < 1e-12
#   PoA | + AD  : OR 55.6 (15.7-280),  p < 1e-12
# With it included, nothing else is significant (age in the AD model drops to
# OR 1.48, 0.95-2.43). Reason for excluding it from the primary models: it explains
# almost everything and changes what the other coefficients mean.
# Writing point: AD and PoA are in practice one decision in this cohort, which is also
# why the results for the three outcomes look so similar.
print(show_or(rd("sensitivity_ADPoA.csv")))
print(show_or(rd("sensitivity_PoAAD.csv")))


# ===========================================================
# 3. Unadjusted associations  (Table1_by_*.csv, univariable_*.csv)
# ===========================================================
# Only nationality differs between groups in all three outcomes:
#   non-German: 1/58 with AD vs 7/36 without (p = .009); PoA 1/54 vs 7/40 (p = .021);
#   BOTH 0/51 vs 8/43 (p = .004).
#   Univariable Firth OR (non-German vs German): AD 0.10 (0.01-0.50), PoA 0.13 (0.01-0.61),
#   BOTH 0.04 (0.0003-0.34). Very wide CIs: only 8 non-German participants.
# Age: AD OR 1.28 per 5 years (1.00-1.66, p = .046); PoA 1.16 and BOTH 1.18, both ns.
# LEDD (AD p = .055) and disease duration (AD p = .063) show trends for AD only.
# All other candidates p > .1.
# With 21 candidates, one or two p < .05 would be expected by chance - but nationality's
# p-values (.001-.008) and its consistency across outcomes make chance less likely.
#
# Table1_by_*.csv: continuous p-values were NA in the first run (tibble issue, fixed in
# the script). Expected after rerun (ANOVA, AD / PoA / BOTH): age .044 / .201 / .161;
# duration .062 / .327 / .200; LEDD .055 / .111 / .103; all others > .25.
for (o in outcomes) {
  cat("\n####", o, "\n")
  print(rd(paste0("univariable_", o, ".csv")) %>%
          mutate(across(c(OR, OR_low, OR_high), ~ round(.x, 2)), p_value = signif(p_value, 2)) %>%
          head(6))
}


# ===========================================================
# 4. Primary model and discrimination  (primaryFirth_*.csv, performance_summary.csv)
# ===========================================================
# Seven predictors, N = 94, EPP 5.1 (AD), 5.7 (PoA), 6.1 (BOTH).
#   AD:   age OR 1.31 per 5 years (1.01-1.72, p = .038); LEDD OR 1.10 per 100 mg
#         (1.00-1.23, p = .052). Sex, marital status, education, preferred place of
#         death and comorbidity ns.
#   PoA:  no predictor significant (age 1.19, 0.93-1.54; LEDD 1.07, 0.99-1.18).
#   BOTH: no predictor significant (age 1.20, 0.94-1.55; LEDD 1.07, 0.99-1.18).
# Discrimination (apparent -> optimism-corrected AUC):
#   AD 0.689 -> 0.592;  PoA 0.642 -> 0.529;  BOTH 0.618 -> 0.505.
# => The pre-specified factors barely discriminate; for PoA and BOTH the corrected AUC
#    is at chance level. This is a finding in itself: having these documents is not
#    explained by the sociodemographic and clinical factors one would expect.
# Note: nationality, the strongest unadjusted correlate, is NOT in the primary set.
for (o in outcomes) { cat("\n####", o, "\n"); print(show_or(rd(paste0("primaryFirth_", o, ".csv")))) }
print(rd("performance_summary.csv"))


# ===========================================================
# 5. Exploratory 21-predictor model  (fullFirth_*.csv)  -  EPP below 2, wide CIs
# ===========================================================
# Consistent across all three outcomes:
#   age        OR 1.48-1.54 per 5 years, p = .010-.017
#   non-German OR 0.02-0.06, p < .005
#   care preference (any place named) OR 6.1-7.4, lower CI 1.2-1.4, p = .017-.030
# PoA additionally LEDD OR 1.14 (1.00-1.31, p = .046).
# Care preference: unadjusted ORs were only 2.1-2.7 and ns; it becomes "significant" only
# after adjustment for 20 others, and only 9 people said "not important". Treat as
# hypothesis-generating at most.
# Use of this model: check that primary estimates keep their direction. They do (age,
# LEDD positive; comorbidity and preferred place of death around 1 or below).
for (o in outcomes) {
  cat("\n####", o, "\n")
  print(show_or(rd(paste0("fullFirth_", o, ".csv"))) %>% filter(p_value < 0.1))
}


# ===========================================================
# 6. Selection stability  (enetStability_*.csv)  and correlations
# ===========================================================
# How to read it:
#   selection_pct = share of half-samples in which lasso kept the term. Robustness of
#   selection, NOT effect size; median_coef/direction only give the sign.
#   Reference values: the noise floor (best of the four noise columns) and the fair share
#   q/p (= mean selection_pct over all columns, what a null column gets on average).
#   The floor is a MAXIMUM over four null columns, so it lies above the fair share even
#   when nothing is going on. The useful comparison is the MEAN of the noise columns vs
#   the fair share, and each real term vs the floor.
#   Correlated predictors split their selection share between them (see below).
#
# Results (floor / fair share / terms above the floor):
#   AD:   21.6% / 11.4%  ->  non-German 51.4 (-), age 37.6 (+), duration 27.2 (+), LEDD 23.6 (+)
#   PoA:  16.2% /  8.7%  ->  non-German 38.4 (-), age 19.2 (+), care preference 19.0 (+)
#   BOTH: 25.2% / 12.0%  ->  non-German 65.4 (-), age 27.6 (+)
# Only nationality is clearly separated from noise (2-3 x the floor) in all three.
# Age is above the floor everywhere but by a small margin except for AD.
# Duration and LEDD are correlated (r = 0.54), so for AD they may be one signal
# ("longer / more advanced disease") split in two.
# Everything else is at or below the noise floor.
#
# Correlations worth mentioning (correlation_matrix.csv):
#   end-of-life thoughts / sharing thoughts r = 0.70;  married / living alone r = -0.65;
#   PDQ / BDI 0.58; nursing support / UPDRS 0.56; LEDD / duration 0.54; UPDRS / PDQ 0.53;
#   age / MoCA -0.45; MoCA / comorbidity -0.42.
for (o in outcomes) {
  d <- rd(paste0("enetStability_", o, ".csv"))
  floor_ <- max(d$selection_pct[d$is_noise])
  cat(sprintf("\n#### %s: noise floor %.1f%%, mean noise %.1f%%, fair share %.1f%%\n",
              o, floor_, mean(d$selection_pct[d$is_noise]), mean(d$selection_pct)))
  print(d %>% filter(!is_noise, selection_pct > floor_) %>% select(term, selection_pct, direction))
}


# ===========================================================
# 7. Synthesis and open points for the discussion
# ===========================================================
# - Main message: in this PD cohort, having an AD or PoA is poorly explained by the
#   measured sociodemographic, clinical and attitudinal factors (corrected AUC 0.51-0.59).
# - The two documents are almost always completed together.
# - The most consistent correlate is nationality (non-German participants rarely held a
#   document), but it rests on 8 people; possible explanations (language, legal
#   familiarity, cultural norms around family decision-making) are speculative.
# - Older age is modestly associated with AD (primary model) and consistently positive
#   elsewhere; disease burden (LEDD/duration) shows a weak tendency for AD only.
# - Attitudinal variables (thought about / shared end-of-life wishes, preferred place of
#   death) show no association - noteworthy, as ACP models assume they precede
#   documentation.
#
# Open points to settle before writing:
# - Primary predictor set: if death_preference_dich replaced care preference /
#   end-of-life thoughts after results were seen, do not call the set "pre-specified".
# - Nationality is the strongest signal but not in the primary model; say so explicitly
#   and justify (too few non-German participants for an adjusted estimate).
