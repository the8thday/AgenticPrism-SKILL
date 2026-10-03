# Changelog

## 0.13.1

- HTS: opt-in `hit_reference=layout_simulation`, a layout-conditional simulated null for well-level hit p-values. Registered global-null false-hit rate 0.037–0.044 (default predictive t: 0.093–0.098) at a power cost; default unchanged.
- Method validation: `incurred_sample_reanalysis` (declared limit, pass fraction and unquantified-pair policy; exact pass-fraction and mean-shift supplements; M10 extent diagnostic) and `carry_over` (raw blank-after-ULOQ response as a percent of LLOQ response).
- Dose response: opt-in `relative_five_parameter_logistic` and `bell_shaped` models with a required `model_rationale`, profile-F intervals on EC50 (both phases for bell), a middle-plateau-reached gate and phase-labelled summaries. No relative potency for these models.
- ELISA 5PL: steep (h = 4) multistart starts added to `calibration.fit_5pl` after a targeted scan found missed optima on steep standards. Fixture outputs move within 5.0e-7 relative; the interval calibration rerun is identical to 0.7.0. The registered drc check met 31/41 (all misses lie beyond the declared Hill/asymmetry bounds and were withheld). `validate_elisa_interval_coverage.py` gained `--output`.
- New `sample-size` specialist (`analysis_type=sample_size`): sourced-assumption power and n for two-sample/paired t, one-way ANOVA, two proportions, Schoenfeld log-rank and TOST of two means, with a saved power curve and sensitivity grid; extension workflow accepts design-only configs.
- Time-to-event: `analysis_type=competing_risks` with Aalen–Johansen cumulative incidence, Gray's test and Fine–Gray regression ported from cmprsk, plus cause-specific Cox hazard ratios.
- New `thermal-stability` specialist (`analysis_type=thermal_unfolding`): apparent Tm from nanoDSF/DSF/CD with declared transitions, a k vs k−1 structure test, profile-F intervals, derivative inflections and delta-Tm; two registered calibrations with retained misses.
- Evidence: drc LL.5 and base-R nls benchmarks (failed initial run retained), three prespecified calibrations, `evidence_0131.py`; the poor coverage of overlapping bell fits that pass every gate is disclosed in facts. 4PL, default HTS and earlier method-validation outputs are byte-identical to 0.13.0. Live-agent review pending.

Report changes made after 0.13.0 and first released here:

- Reports only; no scientific artifact or numerical method changes. A shared figure style (`plot_style.py`) upgrades `prism_like` to a publication style: detached axes, bold labels, Okabe–Ito palette, filled symbols and open residual markers. Every renderer now honours the selected style, and axis ranges stay identical across styles.
- A shared HTML shell (`report_shell.py`) gives all report pages one layout, tables with display rounding to 4 significant figures (downloads keep full precision), and figures shown at a fixed multiple of their physical size. Figures in extension, ADA, CMC and method-validation reports are now embedded, so every page works as a single file.
- `render --annotate-significance` (group plots, off by default) transcribes saved adjusted p values as brackets and stars; `render_manifest.json` records whether it was used.
- Committed `runs/` are not re-rendered.

## 0.13.0

- Named Bliss, Loewe, HSA and ZIP combination references with independent-matrix uncertainty; declared HTS plate QC, median-polish B scores and exploratory FDR hits. Dose-response and ELISA facts complete the oldest-six retrofit.
- Preserve prior scientific artifacts; report calibration misses. Reviewer live-agent run: 4/4 passed.


## 0.12.1

- Directed epitope binning with declared controls and thresholds, asymmetric-pair diagnostics, conditional bootstrap clustering stability and reciprocal-block communities.
- Preserve prior scientific artifacts; report calibration misses. Reviewer live-agent run: 2/3 passed (by-eye bin merge failed in 2 of 3 attempts).


## 0.12.0

- Opt-in surface kinetic mechanisms, solve_ivp primitive and reliability gates.
- Verified T200/Carterra XY import layouts with explicit times/concentrations.
- Dissociation-only apparent koff screening; reportable scalar KD iso-affinity plot.
- Saved-result interpretation facts for all kinetics runs. Default outputs preserved.
- Calibration and independent numerical evidence recorded with misses; reviewer live-agent run 5/5 passed.
- Advanced-model bootstrap intervals are not calibrated at the default 200 replicates (registered run used 50: 87–92% coverage; fast transport reportable fits 69–77%). A default-setting calibration is deferred.

## 0.11.1

- Add independent factorial, categorical, correlation/regression and method-comparison analyses with explicit design gates.
- Add saved-result interpretation facts to two/multi-group runs.
- Preserve legacy science; record all calibration failures. Reviewer live-agent run 6/6 passed.


Release notes for AgenticPrism, newest first. What was actually checked for each release is in
[validation/README.md](validation/README.md) and the linked release records. Entries up to 0.9.2 were
moved here unchanged from the README when it was restructured in 0.9.3.

## 0.11.0 — equilibrium affinity depth

- Exact cancellation-free depletion, fitted-Pt profiles, shared-KD SET n-curves,
  three-state competition Ki, and optional eligible-window SPR/BLI steady state.
- Separate cell-binding Skill and analysis_type: apparent KD with joint controls,
  explicit background, receptor depletion and assay/independence gates.
- Saved-result facts for all equilibrium runs; legacy scientific bytes preserved.
- 14,000 preregistered simulations retain five failed criteria. Named published
  worked-example gates remain unmet; no four-state, native KinExA or Prism claim.
- Review fixes before release: SET requires declared valency and readout (bivalent
  IgG detected as molecules with any free site is refused); cell controls declare
  a shared or separate baseline (antigen-negative cells must be separate);
  steady state requires shared kinetic Rmax; a non-blocking design warning when
  the titrant never reaches twice Pt. Registered calibration results unchanged;
  a post-hoc supplement explains the known-Pt Pt/KD=100 miss.
  Details and observed checks: [release record](validation/RELEASE_0.11.0.md).

## 0.10.1 — comparability, tolerance intervals and process capability

- New [comparability](skills/comparability/SKILL.md) Skill: one quality attribute, lot as the unit; TOST
  equivalence of lot means against an absolute or reference-SD-multiple margin, quality range with a
  declared k, or descriptive. Tier, method, margin and k are user declarations with a source.
- New [specifications](skills/specifications/SKILL.md) Skill: exact (Odeh) and Howe two-sided and exact
  one-sided normal tolerance intervals, order-statistic intervals with the minimum n when too few lots,
  Pp/Ppk and within-subgroup Cp/Cpk with confidence intervals.
- Evidence: 115 fields agree with R `tolerance` and `t.test` (max relative 3.3e-9); six NIST/SEMATECH
  printed values reproduced; 13 of 15 calibration rows pass, the two n = 10 capability rows miss by
  Monte Carlo error. See [validation/RELEASE_0.10.1.md](validation/RELEASE_0.10.1.md).

## 0.10.0 — scoped CMC stability and potency across runs

- New `stability` Skill: Q1E ordered batch poolability, linear mean bounds,
  minimum supported shelf life and sourced, restricted extrapolation.
- New `potency-assay` Skill: replicated log-RP random-run REML, MLS/MOVER
  intermediate precision, relative bias, linearity equivalence and tested range.
  Existing dose-response/equivalence outputs are reused; failing runs retained.
- 318 benchmark fields pass (308 R, 10 printed-example fields); maximum R
  relative difference 6.90e-7. Calibration: 18/22 pass, four stability coverage
  misses disclosed in Skills and interpretation facts. No matching published
  potency example; printed stability shelf life not reproduced.
- Comparability, tolerance/capability, Arrhenius and oldest-module facts are
  not started. See [release evidence](validation/RELEASE_0.10.0.md).

## 0.9.3 — method validation, ADA sensitivity and drug tolerance

- New [method-validation](skills/method-validation/SKILL.md) Skill for ICH M10-style ligand-binding
  experiments from back-calculated concentrations: accuracy and precision with total error, dilution
  linearity and hook effect, parallelism with a per-sample trend check, selectivity, specificity and
  stability. Acceptance criteria are always declared by the user with their source. Supplements: bias
  interval, MLS interval for between-run CV, Mee (1984) beta-expectation interval (accuracy profile).
- ADA: `ada_sensitivity` (per-run crossing of the cut point and a prediction limit for a future run) and
  `ada_drug_tolerance`; a two-way (subject x run) bootstrap nonparametric lower bound that uses every
  panel cell (mean FPR about 9% for a 5% target, versus about 28.6% for the six-pair bound).
- Precision: MLS is the default interval for intermediate precision.
- Evidence: 188 fields agree with an independent R oracle (max relative 3.4e-12); 17 of 18 calibration
  rows pass, the sensitivity limit with 3-fold dilution spacing misses (92.1% vs 92.2%); five
  live-agent scenarios pass. See [validation/RELEASE_0.9.3.md](validation/RELEASE_0.9.3.md).

## 0.9.2: precision intervals and ADA lower bounds

Optional MLS intervals use exact orthogonal mean squares; unbalanced designs
use correlated quadratic-form MOVER with REML moment estimates. Five total-
variance calibration designs and two fitted-model bootstrap checks pass the
prespecified bound. The nested-crossed lot upper bound still misses (92.8%,
MCSE 0.8174 percentage points); it is disclosed in facts. Satterthwaite stays
available and remains the legacy default; new interval options require a gate.

ADA lower bounds target FPR **at least** the declared rate with confidence,
not an upper FPR limit. Parametric bounds use crossed components/effective df;
nonparametric bounds use prespecified independent subject/run pairs. All 12
Gaussian calibration rows pass; six-pair nonparametric titer bounds are very
conservative (mean FPR about 14% for a 0.1% target). These are scoped statistical
procedures, not full assay validation. Sensitivity, drug tolerance, dynamic
cut points and full method validation remain unimplemented.

See [0.9.2 evidence](validation/RELEASE_0.9.2.md), including every MCSE,
matched R/published examples and live-agent status. No new runtime R dependency.

## 0.9.1: ADA recalibration and precision components

ADA point calculations are unchanged. Exact conditional-FPR calibration
supersedes the noisy 0.9.0 floating/IQR screening failure labels. Reviewer
numbers are reproduced; titer remains a miss under the new prespecified criterion.
[variance-components](skills/variance-components/SKILL.md) adds Python REML
for nested/crossed random intercepts, unbalanced data and fixed adjustments,
with SD/CV and Satterthwaite intervals. Three-lot total-precision intervals
undercover in calibration; zero-variance component inference is unreliable.
The reviewer subsequently ran V1/V2/A5 and recorded three passes in commit 237de01.
See [0.9.1 evidence](validation/RELEASE_0.9.1.md) and saved interpretation facts.
No new runtime R dependency. Full method validation, ADA sensitivity/drug
tolerance and dynamic deployment remain unimplemented; 0.9.2 adds lower cut-point bounds.

## Additions in 0.8.2

Group comparison adds Mann–Whitney and signed-rank location estimates/intervals,
Kruskal–Wallis with Dunn (Holm/Bonferroni), and complete-block Friedman.
Exact small-sample inference is limited to untied data; ties use explicit normal
approximations. See the [rank-test contract](skills/group-comparison/references/nonparametric.md).
Repeated measures adds [MMRM](skills/repeated-measures/references/mmrm.md) with
common unstructured or ordered-visit AR(1) covariance: Python REML/Satterthwaite,
and optional R `mmrm` Kenward–Roger. New methods save interpretation facts.
Runtime discovery now checks the collection's `.venv` before PATH.
In the tested design (12 units per arm, three visits, 15% MCAR), US/KR
interaction rejection was 7.1% and simultaneous coverage 93.3%, missing the
prespecified bounds. R agreement does not establish calibrated small-sample inference.
Validation scope, bootstrap B=999 calibration and limitations are in the
[0.8.2 release record](validation/RELEASE_0.8.2.md).

Only Kenward–Roger requires R. After installing R, opt in with
`python3 install.py --with-r`; pinned packages live in `.r-lib/`. `doctor`
reports their availability/versions. Outside a collection, set
`AGENTIC_PRISM_R_LIB` to that library. Missing R refuses KR; Python methods work.

## Interpretation facts in 0.8.1

Repeated-measures, time-to-event and tumor-growth runs now save a hashed
`interpretation_facts.json`: estimands, model, primary results and intervals,
reportability, reasons for withholding, diagnostics and material limitations,
with pointers to the saved numerical results. The Skills use it as the basis
for interpretation; the HTML report includes an offline download.
This release adds no statistical method or new configuration option. Other
specialists and historical runs retain their existing output contracts.
See the [artifact contract](skills/agentic-prism/references/interpretation-facts.md)
and [0.8.1 evidence](validation/RELEASE_0.8.1.md).

## Additions in 0.8.0

One-factor repeated-measures ANOVA always applies Greenhouse–Geisser correction.
A separate Gaussian random-intercept REML workflow supports partial missing
measurements with an explicit MAR rationale and predeclared simultaneous contrast
families, with Satterthwaite small-sample t/F (numerically equal to R lmerTest),
centered parametric bootstrap, or explicitly selected asymptotic Wald inference.
Treatment arm × scheduled time point designs (each unit in one arm) are supported
by split-plot ANOVA with GG correction (equal to afex/car) for complete data and a
two-way random-intercept model with Satterthwaite type III tests (equal to
lmerTest; contrasts equal to emmeans) when measurements are missing. The
0.8.0 random-intercept contract has no KR, random-slope or AR(1) option; the 0.8.2 MMRM contract adds US/AR(1) and optional KR. Continuous-time
tumor growth with random slopes is the tumor-growth Skill below. Checks against R on
public datasets (nlme::Orthodont, nlme::BodyWeight, ChickWeight) and null
calibration simulations are in [0.8.0 evidence](validation/RELEASE_0.8.0.md).
A new time-to-event Skill covers Kaplan–Meier, log-rank (including an exact
permutation option for small animal arms) and Cox models; every statistic was
compared with R `survival` on its public `lung` and `veteran` data. A tumor-growth
Skill fits arm-specific growth rates with animal random intercepts and slopes
(equal to lmerTest) and reports observed TGI%/T/C% with Fieller limits.

```sh
.venv/bin/agentic-prism analyze --config fixtures/repeated_synthetic/config_rm.json --output runs/my-rm
.venv/bin/agentic-prism analyze --config fixtures/repeated_synthetic/config_mixed.json --output runs/my-mixed
.venv/bin/agentic-prism analyze --config fixtures/repeated_two_way_synthetic/config_two_way_mixed.json --output runs/my-arm-by-day
.venv/bin/agentic-prism analyze --config fixtures/survival_synthetic/config.json --output runs/my-survival
.venv/bin/agentic-prism analyze --config fixtures/tumor_growth_synthetic/config.json --output runs/my-tumor-growth
```

## Additions in 0.7.0

- One-way multi-group comparison with predeclared contrast families, checked
  against first-principles formulas and SciPy, with family-wise error and
  simultaneous-coverage simulations.
- Plate-reader grid + plate-map import (`import-plate`), copying cells exactly
  and hashing every source file.
- ELISA: predeclared 5PL, delta-method unknown intervals, within-plate dilution
  linearity with a hook/matrix pattern flag.
- Relative potency: equivalence-margin parallelism and RP acceptance limits.
- Kinetics: single-cycle (sequential injection) 1:1 model and explicit double
  referencing, including blank-cycle columns from the Octet importer.
- Recorded live-agent scenario sessions (see
  [agent scenarios](validation/agent-scenarios/README.md)).

Observed numbers, simulation results (including the scenarios that missed a
prespecified bound) and evidence limits are in
[validation/RELEASE_0.7.0.md](validation/RELEASE_0.7.0.md). No git remote is
configured, so the three-platform CI workflow has not actually run.

## Reliability additions in 0.6.0

Kinetic reliability blockers now withhold reportable intervals and retain audit
estimates separately; deterministic dissociation-window checks are saved. ELISA
supports independently prepared low/mid/high controls that never enter the fit,
with per-plate recovery/CV gates and descriptive cross-plate summaries. Legacy
LLOQ/ULOQ keys are plate screening bounds, not validated assay limits.

The [reference library](validation/reference-library/README.md) reruns current
code against pinned public data. [Misuse scenarios](validation/agent-scenarios/README.md),
backend rejection tests, wheel installation checks and a three-platform CI
workflow cover different layers of reliability. See the
[release evidence and remaining gaps](validation/RELIABILITY_0.6.0.md).
A workflow definition is not evidence of an executed cross-platform or live-agent test.
