---
name: ada-cut-point
description: Estimate declared screening, confirmatory percent-inhibition or titer ADA cut points from complete repeated drug-naive negative panels, with audited normalization, outliers, run diagnostics, subject/run variance components and optional lower confidence bounds; then ADA assay sensitivity and drug tolerance from positive-control experiments against an established cut point. Routes accuracy/precision and other method-validation experiments to method-validation.
---

# ADA cut points

Read the [runtime procedure](../agentic-prism/references/runtime.md), then the
[input and numerical contract](references/contract.md). This is the scoped
ADA implementation, extended in 0.9.3 with a two-way bootstrap lower bound,
positive-control sensitivity and drug tolerance (sections at the end). NAb
assays are unavailable. Accuracy/precision, dilution linearity, selectivity and
stability use the [method-validation](../method-validation/SKILL.md) Skill;
general precision REML uses variance-components.

## Before analysis

Keep conclusions local to the requested assay. Blank inhibited_signal means
no confirmatory readings; it does not identify screening versus titer. Use
the user's declared tier. Describe post-dose interference as possible bias;
neither direction nor magnitude is guaranteed from exposure status alone.


Identify the assay tier, intended negative population, drug exposure, reagent
lot, subject/run/technical-replicate identities and target false-positive rate.
Require an existing protocol for transformations, normalization and exclusions.
Do not populate applicability declarations from an example to make a run work.

- Treated/post-dose samples cannot define the drug-naive negative cut point.
  Disease status alone is not a reason to reject a panel: disease-matched,
  drug-naive representative negative samples may be appropriate. If treatment
  status or pre-existing reactivity is unknown, stop and clarify it.
- Do not delete outliers because the resulting threshold looks better. A
  protocol amendment requires scientific justification independent of a desired
  answer; preserve all prior analyses. No iterative exclusion until normality
  or significance passes.
- A subject measured across runs is one independent subject. Equal technical
  replicates are averaged. Incomplete or unbalanced panels require another
  validated method, not silent complete-case filtering or fabricated wells.
- Do not pool reagent lots without a declared design. This release requires a
  single reagent lot and cannot validate lot comparability. Accuracy/precision,
  total error, dilution linearity, hook effect, parallelism, selectivity,
  specificity and stability are method-validation experiments (route there);
  existing ELISA plate QC is not method validation.
  Several runs nested within each lot do not by themselves make lot and run
  variance mathematically inseparable: a prespecified nested model can estimate
  them with adequate replication. This ADA implementation does not fit that
  model. Run IDs alone give no dates or chronological ordering; do not invent
  confounding with time or analyst from the run names.

## Run and interpret

Use the resolved absolute runtime command to `analyze --config ... --output NEW`.
Then `verify --run ...`; read `interpretation_facts.json`, saved results and
preprocessing. Lead with tier, population, subject/run counts, decision scale,
reportability and supported deployment formula. The reported target false-
positive rate is a protocol input, not a measured guarantee or ADA incidence.

Distinguish additive raw-scale drift (difference normalization) from
multiplicative raw-scale drift (ratio normalization, additive on a log scale).
For this saved result, discuss the observed scale and NC behavior. Do not
claim that both normalization types are log-scale additive or that floating
can never correct raw-scale variance changes.

Check both raw and normalized run diagnostics. Floating requires demonstrated
raw run shift, acceptable normalized mean/variance diagnostics and control
tracking. If mean or variance differences remain after normalization, this
implementation withholds a common threshold. “Dynamic” is a recommendation
for further assay work, not a validated future-run algorithm. Do not promote
per-run audit cut points to reportable thresholds. Non-significant tests do not
prove equality. Median Levene and pooled Shapiro nominal p values ignore the
repeated-subject dependence and are diagnostics, not confirmatory evidence.

Analytical IQR flags pause reporting for investigation; biological IQR handling
is the declared single-pass subject-level policy. Report excluded subject IDs,
counts and rationale. Retain negative raw ANOVA components; clipped values are
not constrained REML. The residual combines subject-by-run interaction and
cell-mean measurement noise, not pure repeatability. Component intervals are
conservative Gaussian ANOVA bounds. Point mode has no cut-point interval; lower mode is described below.

Always convey calibration misses from the facts and
[release evidence](../../validation/RELEASE_0.9.1.md). Numerical R agreement is
not validation of every assay, extrapolated percentile or future run. Do not
turn an ADA positive result into a clinical interpretation without clinical,
exposure and neutralization evidence.

## Sources and verification scope

Methods are motivated by [Devanarayan et al. 2017](https://doi.org/10.1208/s12248-017-0107-3),
[the author's public tutorial](https://bioanalysisforum.jp/images/2024_15thJBFS/D2-A1-01%281%29_JBF15_Viswanath%20Devanarayan_Cut_Point.pdf),
and [FDA 2019 guidance](https://www.fda.gov/media/77796/download).
The contract deliberately implements a narrower procedure. Published numerical
checks use the [rADA vignette](https://cran.r-project.org/web/packages/rADA/vignettes/rada_vignette.html),
whose data are simulated. Do not describe this as a reproduction of every
Shankar/Devanarayan procedure or as Prism equivalence.

## Observed calibration (0.9.1)

Exact conditional FPR replaces the noisy one-future-sample estimate used in
0.9.0. The old floating and biological-IQR screening "failures" are superseded.
Reviewer seed 7 / 400-panel reproduction: fixed and floating 5.4645% with MCSE
0.0948 percentage points; IQR 5.4743% with MCSE 0.0948 percentage points.
The 1000-panel calibration uses the prespecified criterion mean FPR <= 1.2 times
target + 2 MCSE. This tolerates modest point-percentile bias; it does not assure
FPR control for each panel. Titer target 0.1% has mean FPR 0.14886% (MCSE
0.00413 percentage points), exceeding bound 0.12825%: **miss**.
Read all rows, withheld counts and shares with FPR > twice target from facts.
The +0.6 future log-shift scenario is a failure demonstration, not a method miss.
The usual ADA safety direction is a LOWER cut-point bound: lowering the
threshold reduces missed positives and raises FPR. The confidence target is
FPR **at least** the declared target, not FPR below it. If a user requests the
opposite, explain the different objective and ask them to confirm their SOP;
do not silently follow that framing or claim the implemented lower bound caps FPR.

Interpret the threshold mechanics precisely: a ratio-based floating threshold
tracks multiplicative signal drift (an additive shift on log scale); difference
normalization tracks additive drift. Neither guarantees correction of arbitrary
heteroscedasticity. Avoid the general statement that floating can only shift
raw-scale means and never affect raw-scale spread. Missing inhibited readings
do not distinguish screening from titer; use the declared tier. Post-hoc deletion
can bias thresholds and error rates, but neither the new threshold nor a true
future false-positive rate above 5% follows deterministically from the deletion
request. Do not infer a causal exclusion effect from different calibration draws.

## Lower cut-point bound (0.9.2)

Require `bound=lower`, a declared confidence, literal-true applicability and
rationale. Lead with the bound, tier, target, confidence, and direction:
**P(conditional FPR >= target) >= confidence**. It concerns a new independent
subject and random run averaged over the declared run population. It does not
guarantee FPR in each realized run. This separate marginal estimand retains
run diagnostics but does not use the point-mode fixed/floating/dynamic selection
test; do not describe it as a solution for arbitrary dynamic deployment.

Parametric bounds use subject/run components, effective n/df and an approximate
noncentral-t percentile limit. Nonparametric bounds require a prespecified list
of subject/run pairs with distinct subjects AND distinct runs. Repeated cells
cannot inflate the independent sample size. Ask for the pair list and rationale;
never select pairs by response. Outcome exclusions and IQR policies are refused
in this scoped lower-bound implementation.

All 12 Gaussian simulation rows (1000 panels each; 80 subjects, six runs) pass
the prespecified 90% confidence-attainment bound of 88.1026%. This is not a
promise for every assay. Nonparametric titer with six independent pairs has
mean FPR 13.8485% or 13.9725%, despite target 0.1%: the bound is conservative
and may be operationally unsuitable. Discuss the false-positive workload and
additional independent runs with the user; do not portray passing confidence
calibration as tight FPR control. Point-mode titer bias remains a separate result.
Read saved lower-bound facts, effective df/pair count and every calibration row.

[Hoffman–Berger 2011](https://doi.org/10.1016/j.jim.2011.08.019) distinguishes
average and confidence-level cut points;
[Shen et al. 2015](https://pubmed.ncbi.nlm.nih.gov/25356783/) explicitly discusses
lower percentile limits to assure at least the target FPR. The implementation's
mixed-panel approximation is not a reproduction of every method in those papers.

Do not turn an in-panel exceedance fraction or an FPR computed from fitted
Gaussian parameters into an observed future-assay error rate. Label either
as descriptive or model-based and retain its assumptions. For independent-pair
sample-size advice, recompute the admissible rank k at the same confidence and
percentile. The maximum is not necessarily admissible as n grows: quoting
1/(n+1) for the maximum alone does not plan a fixed-confidence lower bound.
Do not infer that nonparametric bounds are generally unsuitable from the
six-run example; state the limited design and operational tradeoff.

Before finalizing an interpretation, check the opening summary and each workload
statement as well as the caveat paragraph. A model-based rate must remain
conditional everywhere: “Under the fitted Gaussian model, the estimated FPR is
X%; future independent validation is needed.” Do not later translate this into
“one in N future negatives will be positive” or “the actual FPR is X%.” The
same restriction applies after saying a number is “for reference.” For saved
nonparametric reports, lead with the saved bound, independent pair count and
release calibration; extra model-based FPR calculations require an explicit
user request and a separately identified model. They are not an observed
performance result of the nonparametric bound.

Distinguish cells sharing a run from the run effects themselves: cells in the
same run share a random effect and are correlated. The model assumes random
effects of distinct runs are independent draws. A test for differences in run
means does not test independence between runs; do not say runs are dependent
because their means differ.

Keep an interpretation scoped to the user's question. Unless sample-size
planning is requested, explain the need for additional independent runs without
adding a sample-size table or calculating alternative ranks. If planning is
requested, k is the largest admissible rank, not the only admissible rank;
smaller ranks can also satisfy the confidence requirement. Invalid pooling
removes the nominal-confidence justification; do not assert a specific actual
coverage failure without a matched simulation. Describe high future workload
as a possibility suggested by the scoped calibration, not a certain outcome.

## Two-way bootstrap lower bound (0.9.3)

For a nonparametric lower bound, prefer `nonparametric_bound: "two_way_bootstrap"`
(the default when no independent pairs are given) with a declared
`bootstrap_seed` and `bootstrap_reps` (default 2000). It resamples subjects and
runs independently (the pigeonhole bootstrap of Owen 2007) and uses every panel
cell, instead of one value per run. In the 0.9.3 calibration (1000 Gaussian
panels per row, 90% confidence) it attained 96.1% (no run effect) and 94.9% (run variance 0.04) with mean FPR
9.0% and 12.3% for a 5% screening target with 80 subjects × 6 runs, versus about
28.5% for the six-pair order statistic. It is conservative, not exact.
It refuses to report when subjects × FPR is below 3 (for example titer 0.1% or
confirmatory 1% with 80 subjects): an extreme empirical tail needs far more
subjects. Report the seed and replicate count; a different seed moves the bound
slightly.

## ADA sensitivity and drug tolerance (0.9.3)

`analysis_type: "ada_sensitivity"` or `"ada_drug_tolerance"`. Both take an
existing cut point with its deployment (`raw`, `ratio`, `difference`) and a
`source` naming where it came from; ask for the ADA cut-point run and do not
estimate a cut point here. Inputs: positive-control (PC) dilution series per
run with `nc_signal` for normalized deployments; for drug tolerance add
`drug_concentration` with a no-drug (0) level per PC level and run.

- **Sensitivity** per run is the PC concentration where the mean response
  crosses the cut point for good, interpolated linearly on log concentration.
  Across runs the tool reports the log-normal **upper prediction limit** for a
  future run (declared `prediction_confidence`, e.g. .95). This is a run-to-run
  consistency figure, not a limit of detection for patient antibodies: PC
  affinity and isotype differ from patients' ADA. Runs positive at the lowest
  or negative at the highest concentration are censored and the prediction limit
  is withheld: extend the dilution range instead of extrapolating.
  Calibration: mean conditional coverage 94.0% (6 runs, 2-fold
  dilutions) and 94.0% (3 runs); with 3-fold spacing 92.1% (bound 92.2%), a
  **miss** against the prespecified bound, because interpolation across coarse
  steps of a sigmoid is biased. Recommend 2-fold spacing around the cut point.
- **Drug tolerance** per PC level and run is the highest drug concentration at
  which the PC still reads at or above the cut point. Results at the edges are
  censored (`positive_at_highest_tested`, `below_lowest_nonzero_tested`,
  `not_detected_without_drug`); summaries appear only when every run is
  interpolated. Report the per-run statuses as the result otherwise. State the
  PC level with each tolerance value: tolerance rises with PC concentration and
  one number without its PC level is not meaningful.

Mention every `must_mention` item (for example multiple cut-point crossings).
