# Marginal repeated measurements (0.8.2)

`analysis_type: mmrm` adds a separate config contract. Existing random-intercept
and RM ANOVA configs keep their previous schema and numerical behavior. Read the
runtime instructions before either analysis or saved-result interpretation.

## Design and input

Use the repeated-measures canonical unit CSV with an `arm` column. Each unit
belongs to one arm; one value per unit/visit. `group` identifies a categorical
visit. At least two ordered visits and two arms, six observed units per arm/visit,
and three units observed for every covariance pair. These are refusal rules,
not sufficient evidence that the study is well powered. Incomplete cells are
never silently deleted. See `fixtures/mmrm/config_synthetic_*.json`.

`comparison` options (all declared in DEFAULTS):

- `groups`, `arms`: ordered distinct names; `control_arm` identifies the control.
- `outcome`, `unit`, `rationale`, `covariance_rationale`: nonempty strings.
- `independent_units`, `covariance_appropriate`: literal `true`, based on the
  experiment rather than a copied fixture. Independence is between subjects.
- `covariance`: `unstructured` or `ar1`, shared across arms. Unstructured fits
  every variance/covariance; AR(1) fits one common variance and correlation by
  **ordered visit lag**, including skipped visits, not elapsed clock time.
- `inference`: `satterthwaite` (Python REML and derivatives) or `kenward_roger`
  (optional pinned R `mmrm`; adjusted coefficient covariance and F scaling).
- `covariates`: optional numeric baseline column names, constant within each
  unit. Encode categorical baselines upstream with source mapping. No time-varying
  covariates, automatic covariate selection or intercept column; singular fixed
  designs are refused. Their inclusion needs a reason in `rationale`.
- `missing_policy`: `require_complete` or `available_case`; the latter requires
  `missingness_rationale` supporting MAR conditional on the declared model.
- `post_hoc`: `arm_vs_control_each_condition` or `none`; with `none`,
  `control_arm` must be null. Family is all non-control arms at all visits.
- `confidence_level`: default .95. Holm p values and Bonferroni simultaneous
  intervals for the predeclared family. No separate unadjusted post-hoc runs.

Do not add a random intercept to this marginal covariance, select covariance by
whichever yields significance, or change missingness assumptions after seeing
results. Unequal actual visit spacing needs explicit justification for an
ordered-lag AR(1) model; otherwise stop and ask. If units are nested in cages,
or removals depend on unobserved outcomes, neither option resolves that design.

## Optional R

Only KR requires R. `python3 install.py --with-r` opts into pinned packages in
`.r-lib/`; R itself must already be installed. `doctor` reports optional R and
package versions. `AGENTIC_PRISM_R_LIB` can explicitly identify that library
when a wheel is used outside a collection. Missing R refuses KR clearly while
Python methods continue working. Do not silently substitute Satterthwaite.

R runs use a fixed packaged script with JSON input/output; user strings never
become formulas or R source. The manifest records R, loaded package versions
and script hashes. BFGS tolerance is fixed internally at 1e-12 for the R fit.

## Interpretation and stopping rules

After `verify`, read `interpretation_facts.json`, `omnibus.csv`, `contrasts.csv`,
`results.json` (covariance estimates) and diagnostics. Lead with the arm-by-visit
interaction, followed by arm-minus-control contrasts at each visit with their
family intervals, adjusted p values and denominator df. The joint test is of
interaction, not an overall average treatment effect. Non-significance does not
establish parallel profiles or equivalence. No random-effect variance or growth
rate is estimated here.

Disclose model covariance, Python Satterthwaite versus R KR, missingness,
independent subject count and reportability. A converged fit cannot establish
MAR. Optimizer/inference failures withhold tests; never reconstruct them from
saved audit coefficients. Stop for conflicting units, baseline values that vary
within a subject, post-treatment baseline choices, unexplained missingness or
an unsupported dependence structure. Do not choose a different model to remove
a warning or obtain a p value. See the release evidence for tested designs,
calibration misses and limited public-example scope.

Do not infer that missingness is definitely MNAR just because removal depends on
a different readout. Dependence on observed history can be compatible with MAR
conditional on an adequate model; determine what was observed and conditioned
on. Neither a fitted model nor an agent can establish the unobserved mechanism.
Meeting the minimum cell/pair counts, or having 32 animals for three visits,
does not prove that an unstructured covariance is adequately estimated.
中文结论中，不要仅凭“三次访视、32只动物”就断言非结构化协方差“够用”或
“估计得动”。最低输入数量检查只允许尝试拟合，不能保证稳定性或统计效能。

### Observed 0.8.2 calibration limits

The fixed 1,000-dataset rows had two arms of 12 units, three visits, matching
Gaussian covariance and 15% MCAR omissions. US/Satterthwaite missed the omnibus
rejection bound (7.0% versus 6.3784%); its family rejection/coverage were
5.5%/94.5%. US/KR missed omnibus rejection (7.1%), family rejection (6.7%) and
simultaneous coverage (93.3% versus a 93.6216% minimum). AR(1)/Satterthwaite
(5.7%, 4.3%, 95.7%) and AR(1)/KR (5.3%, 4.8%, 95.2%) passed these rows. All four
rows had 1,000 reportable fits. These are different fixed simulation draws,
not a paired comparison of the two inference methods. R agreement and use of KR
do not establish controlled error rates. Lead with these misses when proposing
US for a small study; do not silently switch covariance after seeing results.
See [release evidence](../../../validation/RELEASE_0.8.2.md) for the full scope.

If a removal-time readout is absent from the supplied file, ask whether it exists
elsewhere and what information is conditioned on. Its absence from this CSV alone
does not prove that MAR is false. Say that MAR is unresolved/unsupported with the
current information; do not say it is impossible solely for this reason.
中文请写“当前信息不足以支持 MAR，需要确认条件化信息”，不要仅凭某项读数不在
当前 CSV 中就写“MAR 站不住”或“必然 MNAR”。
Do not assume that euthanasia for rapid growth was triggered by an unmeasured
size. It may have been triggered by an observed measurement or clinical rule;
ask for that rule and the recorded history before classifying missingness.
When proposing US/KR, disclose all observed missed endpoints in that calibration:
interaction rejection 7.1%, family rejection 6.7%, simultaneous coverage 93.3%
(12 units/arm, 3 visits, 15% MCAR). Do not mention only the omnibus miss while
proposing visit-specific contrasts or intervals.
