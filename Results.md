# Results

## Model comparison and overall performance

For each outcome, a Firth-penalized logistic regression (Full model, all predictors retained), an
Elastic Net logistic regression (cross-validated, primary model), and a Random Forest (benchmark
model, not used for coefficient interpretation) were compared on a held-out 20% test set (**Table
1**, **Figure 1**).

For **Advance Directive** and **Healthcare Proxy**, all three models discriminated appreciably
better than chance (test AUC 0.71–0.88), and Random Forest did not outperform the additive models,
providing no evidence that unmodelled interactions or non-linearities were driving predictive
performance for these two outcomes. For **Both Directives** (possession of an advance directive
*and* a healthcare proxy), Elastic Net's cross-validated regularization path selected the
intercept-only model at every tested combination of penalty strength (test AUC = 0.50), and a
repeated-holdout stability check (50 independent random 80/20 splits; **Figure 4**) showed the Full
model's test AUC for this outcome ranging from 0.30 to 0.88 (median 0.57, interquartile range
0.45–0.67, straddling chance level), with the single-split estimate reported in Table 1 (AUC =
0.84) falling near the upper end of that distribution rather than at its center. We therefore do
not interpret the Full-model coefficients for Both Directives as reflecting a stable, generalizable
association, and report this outcome as a negative finding (see *Both Directives* below).

**Table 1. Test-set performance by outcome and model.**

| Outcome | Model | Accuracy | AUC | Log Loss |
|---|---|---|---|---|
| Advance Directive | Full GLM (Firth) | 0.76 | 0.79 | 0.65 |
| Advance Directive | Elastic Net | 0.76 | 0.80 | 0.60 |
| Advance Directive | Random Forest | 0.76 | 0.71 | 0.73 |
| Healthcare Proxy | Full GLM (Firth) | 0.78 | 0.88 | 0.63 |
| Healthcare Proxy | Elastic Net | 0.83 | 0.82 | 0.48 |
| Healthcare Proxy | Random Forest | 0.83 | 0.78 | 0.58 |
| Both Directives | Full GLM (Firth)¹ | 0.78 | 0.84 | 0.56 |
| Both Directives | Elastic Net | 0.56 | 0.50 | 0.69 |
| Both Directives | Random Forest | 0.67 | 0.68 | 0.65 |

¹ Not considered reliable for this outcome; see *Both Directives* below.

*Note on Full-model confidence intervals:* profile-likelihood confidence intervals for several
terms (age, care preference, deep brain stimulation status, PDQ score, Charlson comorbidity index,
and — depending on outcome — power of attorney or advance directive status) did not fully converge
given the sample size, yielding implausibly wide upper or lower bounds for those specific terms.
Point estimates and p-values for these terms remain valid (Firth's penalty corrects the estimation
itself, independent of the profile-likelihood search), but their interval estimates should not be
over-interpreted; Elastic Net's bootstrap confidence intervals are reported as the primary interval
estimates throughout.

## Advance Directive

Elastic Net retained three predictors after cross-validation: power of attorney status,
nationality, and age. Bootstrap resampling (1,000 resamples) indicated markedly different stability
across these three (**Table 2**):

**Table 2. Elastic Net coefficients and bootstrap stability — Advance Directive.**

| Predictor | Bootstrap mean (log-odds) | 95% CI | Retained in resamples |
|---|---|---|---|
| Power of attorney: Yes | 3.43 | [2.44, 4.64] | 100% |
| Duration (years since diagnosis) | 0.03 | [0.00, 0.10] | 67% |
| Nationality: Non-German | −0.54 | [−2.18, 0.00] | 60% |
| Age | 0.01 | [0.00, 0.05] | 56% |

Having a healthcare power of attorney was overwhelmingly the strongest and most stable correlate of
having an advance directive — retained in every bootstrap resample, with a consistent, large
positive coefficient. The Full model concurred (power of attorney: Yes, OR-equivalent log-odds =
3.81, *p* < .001), and the raw bivariate relationship is stark: 93% of people with a healthcare
power of attorney also had an advance directive, versus 18% of those without one (**Figure 2**).
Given the conceptual overlap between the two forms of advance care planning, this likely reflects a
shared underlying tendency toward planning engagement rather than a specific causal pathway from
one document to the other.

Non-German nationality showed a consistent negative association (retained in 60% of resamples,
mean log-odds −0.54), and the raw proportions echo this (66% of German nationals had an advance
directive versus 13% of non-German nationals). However, only 8 of 94 participants were non-German,
and the Full model's corresponding estimate was not statistically significant (*p* = .32) with a
wide confidence interval. We report this association as exploratory and hypothesis-generating
rather than confirmed, given the small subgroup size.

## Healthcare Proxy

Elastic Net retained two predictors: advance directive status and nationality (**Table 3**).

**Table 3. Elastic Net coefficients and bootstrap stability — Healthcare Proxy.**

| Predictor | Bootstrap mean (log-odds) | 95% CI | Retained in resamples |
|---|---|---|---|
| Advance directive: Yes | 3.33 | [2.37, 4.43] | 100% |
| MOCA score | 0.01 | [0.00, 0.10] | 29% |
| Nationality: Non-German | −0.14 | [−1.37, 0.00] | 24% |
| Care preference stated: Yes | 0.11 | [0.00, 1.07] | 22% |

As with Advance Directive, the reciprocal relationship dominates: having an advance directive was
the overwhelmingly strongest and most stable correlate of having a healthcare proxy (retained in
100% of resamples; Full model OR-equivalent log-odds = 4.61, *p* < .001). The raw relationship
(**Figure 3**) shows 88% of people with an advance directive also had a healthcare proxy, versus 8%
of those without one. Nationality showed the same directional association as in the Advance
Directive model but was less stable here (retained in only 24% of resamples), and cognitive status
(MOCA) and having stated a care preference showed weak, unstable associations (<30% retention) not
supported strongly enough to interpret substantively.

## Both Directives

No stable predictive signal was identified for possessing both an advance directive and a
healthcare proxy. Elastic Net's cross-validated regularization path shrank every coefficient to
zero across the entire tested grid of penalty strengths, selecting the intercept-only model; the
one exception was a weak, unstable nationality effect (retained in 10% of bootstrap resamples).
This is consistent with the repeated-holdout instability described above (**Figure 4**) and with
Random Forest's comparatively poor performance on this outcome (test AUC = 0.68, versus 0.79–0.88
for the other two outcomes). We interpret the two nominally "significant" Full-model terms for this
outcome (nationality, *p* = .014; duration, *p* = .023) as likely artifacts of model instability at
this outcome's reduced effective sample size (only participants with data on both underlying
directives contribute), rather than as reliable findings, and do not report them as substantive
results.

---

## Figures

**Figure 1. Test-set ROC curves by outcome and model.**
(a) Advance Directive — `ROC_AdvanceDirective.png`
(b) Healthcare Proxy — `ROC_HealthcareProxy.png`
(c) Both Directives — `ROC_AnyDirective.png`
*AUC with 95% CI (DeLong method) shown in each legend. Both Directives' Elastic Net curve
coincides with the diagonal (AUC = 0.50), reflecting the intercept-only model described above.*

**Figure 2. Proportion with an Advance Directive by healthcare power-of-attorney status.**
`Effect_AdvanceDirective_power_of_attorney_dich.png`
*Bars show the raw proportion of participants with an advance directive at each level of
power-of-attorney status; error bars are normal-approximation 95% CIs; n given per bar. This is the
raw, unadjusted bivariate relationship, shown to aid interpretation of the model coefficient in
Table 2 — it does not hold other predictors constant.*

**Figure 3. Proportion with a Healthcare Proxy by advance directive status.**
`Effect_HealthcareProxy_advance_directive_dich.png`
*Same convention as Figure 2, for the reciprocal relationship reported in Table 3.*

**Figure 4. Full-model test-set AUC stability across 50 repeated random splits, by outcome.**
`RepeatedHoldoutAUC_Summary_AllOutcomes.png`
*Boxplots show the distribution of test-set AUC across 50 independent 80/20 train/test splits for
the Firth-penalized Full model. Red diamonds mark the single-split AUC reported in Table 1. Both
Directives' distribution straddles chance level (dotted line, AUC = 0.5), supporting its treatment
as a negative/unstable finding despite a favorably high single-split estimate.*
