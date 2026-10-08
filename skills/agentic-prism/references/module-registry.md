# Module registry — development 0.13.4

| Purpose | Skill | Capability state |
|---|---|---|
| Two-agent combination reference scores | [drug-combination](../../drug-combination/SKILL.md) | Experimental code, gates incomplete; 0.13.0: Bliss/Loewe/HSA/ZIP, existing 4PL fitter, declared scales and independent-matrix uncertainty; unsupported fits and unreplicated claims withheld; see observed evidence limits |
| HTS plate QC and exploratory hits | [import-plate](../../import-plate/SKILL.md) | Experimental code, gates incomplete; 0.13.0: explicit grid maps, Z-prime, SSMD, median-polish B scores and nominal predictive-t/BH screens; global-null false-hit calibration0.100 exceeds0.063784, so validated FDR control is UNMET for this default; failed QC withholds hits; no confirmed biological activity claim. 0.13.1 opt-in `hit_reference=layout_simulation`: false-hit 0.037–0.044 in registered Gaussian additive-plate rows, at a power cost; published example gate unmet; live-agent scenarios passed (0.13.1) |
| Directed epitope competition and bins | [epitope-binning](../../epitope-binning/SKILL.md) | Experimental code, gates incomplete; 0.12.1: explicit controls/thresholds, directed asymmetry, average linkage with conditional BP stability, reciprocal-block components; no structural epitope or AU claim; published controlled-response example UNMET |
| Single-site equilibrium binding, KD and diagnostics | [equilibrium-binding](../../equilibrium-binding/SKILL.md) | Implemented; see release validation record |
| Re-render an existing equilibrium run | Same specialist, render command | Implemented; no refitting |
| BLI/SPR 1:1 kon/koff and kinetic KD | [binding-kinetics](../../binding-kinetics/SKILL.md) | Implemented for independent cycles and (0.7.0) single-cycle series in canonical CSV, explicit reference-channel or double referencing, and the scoped Octet Results.txt layout; reliability gates with audit-only intervals for blocked fits. 0.12.0 experimental code (evidence gates incomplete) adds opt-in drift, heterogeneous ligand, bivalent analyte, transport and off-rate screening, scoped T200/Carterra XY imports and saved-result facts. Native binary formats remain unsupported |
| ELISA 4PL/5PL calibration and unknown concentration interpolation | [elisa-quantification](../../elisa-quantification/SKILL.md) | Implemented for per-plate standards with back-calculation QC, independent-control gates, descriptive cross-plate QC, range-checked inverse and dilution, (0.7.0) predeclared 5PL, calibration-conditional delta-method unknown intervals and within-plate dilution linearity. No mixed-model intermediate precision or matrix model |
| Relative EC50/IC50 from symmetric 4PL curves; replicate summaries; reference-vs-test relative potency | [dose-response](../../dose-response/SKILL.md) | Implemented for long-form CSV with explicit endpoint/direction/units; optional fixed plateaus and relative weighting; parallel-line RP with F-test or (0.7.0) predeclared-margin equivalence parallelism and optional RP acceptance limits. 0.13.1 opt-in asymmetric 5PL and bell-shaped models (declared rationale, no RP; drc/nls checks and calibration passed; overlapping-phase bell fits that pass the gates cover poorly; published example gate unmet; bell-shaped live-agent scenarios passed, 5PL not scenario-tested) |
| Single-sample and positive ratio tests | [group-comparison](../../group-comparison/SKILL.md) | 0.13.4 `location_test`: ordinary one-sample t; one-sample/paired log ratio and independent Welch log ratio; four base-R t.test comparisons pass; one-sample log-ratio coverage miss retained. No pseudocounts or imaginary controls. See [evidence](../../../validation/RELEASE_0.13.4.md) |
| Relative ECx/ICx and 4PL adequacy | [dose-response](../../dose-response/SKILL.md) | 0.13.4: declared percentages including EC10/20/80/90 and IC90; endpoint-specific profile-F intervals (pointwise), range and parent-fit withholding; symmetric 4PL only. Approximate replicated lack-of-fit F plus residual-sign diagnostic; passing does not prove adequacy. Valid endpoint base-R checks pass; audit-only hook discrepancies and narrow-range C80 undercoverage retained. [Evidence](../../../validation/RELEASE_0.13.4.md) |
| Two-group independent or paired comparisons | [group-comparison](../../group-comparison/SKILL.md) | Implemented: Welch or paired two-sided t test on independent-unit data |
| One-way comparison of ≥3 independent groups | [group-comparison](../../group-comparison/SKILL.md) | Implemented (0.7.0): classic or Welch ANOVA with predeclared Dunnett, Tukey-Kramer, Games-Howell or Holm-Welch family. Parametric independent-group contract; rank tests use the separate entry below |
| Repeated measurements: one within-unit factor, or between-unit arm × within-unit condition | [repeated-measures](../../repeated-measures/SKILL.md) | Implemented (0.8.0): complete RM ANOVA and split-plot ANOVA with GG correction; REML random-intercept models (one-factor and arm × condition) with explicit available-case policy and Satterthwaite (lmerTest-equivalent), parametric bootstrap (one-factor only) or Wald inference. No random slopes (see tumor-growth for continuous-time growth), KR, nested/crossed effects or serial correlation. See 0.8.0 evidence |
| Install / locate / check the runtime | `install.py`, `agentic-prism doctor` ([runtime instructions](runtime.md)) | Implemented (0.7.1): per-clone `.venv` with lock-constrained dependencies and self-check; macOS evidence only |
| Competing risks | [time-to-event](../../time-to-event/SKILL.md) | 0.13.1 (`analysis_type=competing_risks`): Aalen–Johansen CIF with log(−log) intervals, Gray's test, Fine–Gray (cmprsk crr sandwich) and cause-specific Cox; equal to cmprsk/survival within 7e-14 (incl. mgus2); 6/6 calibration rows passed; one censoring group, no strata; live-agent scenarios passed (0.13.1) |
| Time-to-event (survival, time to humane endpoint) | [time-to-event](../../time-to-event/SKILL.md) | Implemented (0.8.0): Kaplan–Meier with log/log-log limits, medians, landmark survival, log-rank (asymptotic or permutation) with pairwise Holm families, Cox (Efron) hazard ratios with optional covariates, cox.zph PH test; equal to R survival. No competing risks, recurrent events, time-varying covariates, frailty or interval censoring |
| Longitudinal tumor volumes (in-vivo efficacy) | [tumor-growth](../../tumor-growth/SKILL.md) | Implemented (0.8.0): log-volume random-slope mixed model with Satterthwaite (equal to lmerTest) growth rates, rate differences, doubling times and model T/C; observed TGI%/T/C% with Fieller limits and humane-endpoint dropout diagnostics. No Gompertz/logistic curves, cage effects or AR(1) |
| Plate-reader grid + plate map to long table | `agentic-prism import-plate` (used by the ELISA and dose-response Skills) | Implemented (0.7.0) for a generic 96/384-well grid CSV layout; no vendor-specific export parser |

| Rank-based unit-level comparisons | [group-comparison](../../group-comparison/SKILL.md) | 0.8.2: MW/signed-rank with HL and exact untied or normal approximate intervals; KW/Dunn Holm/Bonferroni; complete-block Friedman. See nonparametric contract and 0.8.2 evidence |
| Marginal repeated measurements | [repeated-measures](../../repeated-measures/SKILL.md) | 0.8.2: common US or ordered-visit AR(1), optional numeric baseline covariates; Python REML/Satterthwaite or optional R mmrm KR. Separate `analysis_type=mmrm` contract; interpretation facts included. US small-sample calibration misses are recorded in 0.8.2 evidence |

| ADA cut points | [ada-cut-point](../../ada-cut-point/SKILL.md) | 0.9.0 scoped implementation: complete balanced repeated negative panels, screening/confirmatory/titer point percentiles, normalization, prespecified transformations/outliers, subject/run ANOVA components. Fixed/floating with gates; dynamic recommendation with withheld deployment. Current calibration is superseded by 0.9.1 exact conditional FPR; titer miss disclosed. No sensitivity, drug tolerance or full method validation |

| ADA lower bound, sensitivity, drug tolerance | [ada-cut-point](../../ada-cut-point/SKILL.md) | 0.9.3: two-way (subject × run) bootstrap nonparametric lower bound with a tail-support gate; `ada_sensitivity` (per-run log-linear crossing, log-normal prediction limit) and `ada_drug_tolerance` (per PC level and run, censored at the edges). See 0.9.3 evidence |
| Bioanalytical method validation (LBA) | [method-validation](../../method-validation/SKILL.md) | 0.9.3: accuracy/precision with total error, dilution linearity and hook effect, parallelism with trend, selectivity, specificity, stability from back-calculated concentrations against user-declared criteria; bias, MLS and beta-expectation supplements. 0.13.1 adds incurred sample reanalysis (declared limits and unquantified-pair policy, calibrated supplements) and carry-over from raw responses; ISR live-agent scenario passed, carry-over not scenario-tested. No calibration-curve fitting |
| Thermal unfolding (apparent Tm) | [thermal-stability](../../thermal-stability/SKILL.md) | 0.13.1: two-state transitions with sloping baselines (1–3 declared), k vs k−1 structure test, profile-F Tm, derivative inflections, Welch delta-Tm; nls agreement 2.4e-13; two registered calibrations (run 1: 3 misses; run 2: 1 miss and a 0.9% overlap leak, all retained); DSC and light scattering refused; live-agent scenarios passed (0.13.1) |
| Nested replicates within independent units | [group-comparison](../../group-comparison/SKILL.md) | 0.13.2 (`analysis_type=nested_comparison`): random-intercept REML with Satterthwaite (lmerTest within 1.7e-7), unit-means analysis always reported and primary at a boundary fit, ICC and design effect; 5/5 calibration rows passed; one nesting level; live-agent scenarios passed (0.13.2) |
| Baseline-adjusted comparison (ANCOVA) | [group-comparison](../../group-comparison/SKILL.md) | 0.13.2 (`analysis_type=ancova`): pre-treatment covariates only, common-slope OLS, slope-homogeneity gate, adjusted means and contrasts (lm/emmeans/car within 4.2e-12); 4/4 calibration rows passed; live-agent scenarios passed (0.13.2) |
| Area under time curves | [curve-auc](../../curve-auc/SKILL.md) | 0.13.2: linear trapezoid per unit on a declared interval, declared baseline and early-end policy, unequal-withholding flag, Welch comparisons; 3/3 calibration rows passed; not PK/NCA; live-agent scenarios passed (0.13.2). 0.13.4 adds paired t contrasts with explicit pair policy, Holm p and Bonferroni family intervals; 4/4 paired calibration rows passed; see [new evidence](../../../validation/RELEASE_0.13.4.md) |
| Linear/quadratic standard curves | [standard-curve](../../standard-curve/SKILL.md) | 0.13.2: declared form and weighting, recovery acceptance, lack of fit, monotonicity, inversion intervals (investr within 5.1e-8), no extrapolation; 2/2 calibration rows passed; live-agent scenarios passed (0.13.2) |
| qPCR relative quantification | [qpcr](../../qpcr/SKILL.md) | 0.13.2: efficiency-corrected ΔΔCq with geometric-mean references, biological-replicate Welch/paired t, NTC/undetermined/technical-SD/reference-shift flags (base-R within 1.9e-13); 3/3 calibration rows passed; published example unmet; live-agent scenarios passed (0.13.2) |
| Common nonlinear models | [nonlinear-models](../../nonlinear-models/SKILL.md) | 0.13.3 (`analysis_type=nonlinear_fit`): one/two-phase decay, association, exponential and logistic growth, Michaelis-Menten; profile-F intervals, derived half-lives, design-support gates; base-R nls/confint agreement 18/18; 11/11 calibration rows passed; live-agent gate not run |
| Prospective sample size and power | [sample-size](../../sample-size/SKILL.md) | 0.13.1: two-sample/paired t, one-way ANOVA, two proportions (approximate), log-rank events (Schoenfeld), TOST for two means; sourced assumptions; pwr/PowerTOST agreement 9.5e-10, simulated power within bounds for 10/10 rows; no repeated-measures, clustered or multiplicity designs; live-agent scenarios passed (0.13.1) |
| Precision variance components | [variance-components](../../variance-components/SKILL.md) | 0.9.1–0.9.3: REML nested/crossed components; MLS/MOVER intervals for sums are the default from 0.9.3 |

Published validation evidence: [validation/README.md](../../../validation/README.md).

0.8.1 adds [interpretation facts](interpretation-facts.md) to repeated-measures,
time-to-event and tumor-growth runs. Other specialists retain their existing
contracts. No new statistical method was added in 0.8.1; at that release, nonparametric tests,
MMRM, KR, Gompertz, cage effects, MNAR sensitivity, competing risks and RMST
remained follow-up work. The 0.8.2 entries above supersede that earlier scope. See [0.8.1 evidence](../../../validation/RELEASE_0.8.1.md).

“Implemented” is not Prism numerical equivalence. No Prism project benchmark has
been run. New modules must supply their own data contract, estimand, assumptions,
tests and release status. A requested mixture of available and unavailable tasks
may complete the available work while explicitly identifying the remaining gap.

Route sensorgrams to the kinetic specialist; its Skill records current limits. Route
functional EC50/IC50 curves to the 4PL specialist, not the equilibrium KD fitter.

Routing examples:

- “Here is a titration table, calculate KD”: specialist checks equilibrium,
  concentration basis, signal meaning and valency before choosing the model.
- “Here is an ELISA plate, infer unknown concentrations”: calibration specialist;
  do not send its OD values through the KD workflow.
- “Compare three independent KD experiments”: equilibrium specialist summarizes
  log KD descriptively. For a predeclared two-group contrast, use the group
  specialist on one independent-unit value per group/experiment.
- “Compare four candidate antibodies with the isotype control”: group specialist,
  multi-group design with a vs-control family chosen before testing.
- “Same donors under three treatments”: repeated-measures specialist; establish
  complete one-factor design or justified random-intercept model. Do not use
  independent-group ANOVA or unadjusted paired t families.
- “Body weight of three dose groups over five weeks”: repeated-measures
  specialist, arm × condition design; report the interaction first.
- “Did the antibody slow tumor growth? Give TGI.”: tumor-growth specialist;
  lead with the growth-rate model, treat observed TGI as secondary when animals
  were removed early, and pair with time-to-event for the humane endpoint.
- “Which dose prolonged survival / delayed the humane endpoint?”: time-to-event
  specialist; establish event definition, time zero and censoring first.
- “Compute kon/koff from an endpoint table”: insufficient time-course data;
  request sensorgrams rather than manufacturing rate constants.
- “Single-cycle kinetics from Biacore”: kinetic specialist with
  `injection_design=single_cycle` after mapping the export to the canonical CSV.
- “Estimate IC50 from viability versus inhibitor concentration”: dose-response
  specialist; determine whether relative plateau midpoint is the intended IC50.
- “Is our candidate more potent than the reference antibody?” / relative potency:
  dose-response comparisons, pairing reference and test curves from the same
  plate, one comparison per independent experiment. Equivalence margins and RP
  acceptance limits must come from the user's historical data or protocol.
- “Make existing plots Prism-like”: render-only command on a saved run.

Current reliability scope: [0.6.0 evidence](../../../validation/RELIABILITY_0.6.0.md)
and [0.7.0 additions](../../../validation/RELEASE_0.7.0.md).

Repeated-measures addition: [0.8.0 evidence](../../../validation/RELEASE_0.8.0.md).

0.8.2 implements the rank-test and MMRM entries above; the 0.8.1 paragraph describes that earlier release. See [0.8.2 evidence](../../../validation/RELEASE_0.8.2.md).

## Precision components (0.9.1)

| Purpose | Skill | Capability state |
|---|---|---|
| General precision components | [variance-components](../../variance-components/SKILL.md) | 0.9.2: optional MLS/MOVER sums and one-sided component upper bounds; 0.9.1 reviewer scenarios passed. Python REML: nested/crossed random intercepts, unbalanced data, fixed adjustments; Satterthwaite component/total intervals, SD/CV. Small three-lot total coverage and true-zero component calibration misses disclosed. Not full method validation |

## 0.9.2 additions

ADA adds opt-in parametric or independent-pair order-statistic lower percentile
bounds: FPR at least target with declared confidence, marginal over subject/run.
Precision adds opt-in mean-square MLS/MOVER for total/positive sums and separate
one-sided component upper bounds. See [0.9.2 evidence](../../../validation/RELEASE_0.9.2.md)
for approximation boundaries, the nested-crossed lot upper-bound miss and the
conservative six-pair nonparametric titer result. No full method validation,
sensitivity/drug tolerance or dynamic deployment. New methods run in Python.

## CMC (0.10.0–0.10.1)

| Purpose | Skill | Capability state |
|---|---|---|
| Quantitative long-term stability and shelf life | [stability](../../stability/SKILL.md) | Implemented, scoped linear fixed-batch Q1E regression, ordered ANCOVA alpha .25, one/two-sided 95% mean bounds; sourced specification/direction/storage; restricted declared extrapolation. Four selection-related calibration misses. Published regression/ANCOVA reproduced, printed shelf-life gate not met |
| Replicated RP combination and nominal-level validation | [potency-assay](../../potency-assay/SKILL.md) | Implemented, Python random-run REML and MLS/MOVER; sourced criteria, per-run suitability, retained outlier flags, bias/precision/linearity/tested range. Requires at least two independent determinations per run per level. Published-example gate not met |
| Comparability / biosimilarity | [comparability](../../comparability/SKILL.md) | 0.10.1: one attribute, lot as unit; Welch or pooled TOST against an absolute or reference-SD-multiple margin, quality range with declared k and optional required fraction, or descriptive; tier recorded, never inferred. Not a totality-of-evidence judgement |
| Tolerance intervals, Cp/Cpk/Pp/Ppk | [specifications](../../specifications/SKILL.md) | 0.10.1: normal exact (Odeh) or Howe two-sided and exact one-sided tolerance intervals; order-statistic intervals with minimum n; Pp/Ppk overall and Cp/Cpk within subgroups with chi-square / normal-approximation intervals. n = 10 capability rows missed by Monte Carlo error (disclosed) |
| Remaining old-module facts; Arrhenius | None | All six oldest modules write facts as of 0.13.0. Arrhenius not started |

Evidence and limitations: [0.10.0 release](../../../validation/RELEASE_0.10.0.md), [0.10.1 release](../../../validation/RELEASE_0.10.1.md).

## 0.11.0 affinity additions (implemented with evidence limitations)

- `cell_binding`: separate [cell-binding](../../cell-binding/SKILL.md) specialist,
  apparent KD with joint nonspecific controls (declared shared or separate control
  baseline; antigen-negative cells separate) and receptor-depletion declaration.
  Published worked-example gate UNMET. See [release evidence](../../../validation/RELEASE_0.11.0.md).

- `equilibrium_binding` opts into `one_site_depletion`,
  `solution_equilibrium_titration`, or `competition_exact`; shared mass balance,
  profile-F, active-site/equilibration gates, SET valency/readout gate (bivalent
  any-free-site capture refused), non-blocking titration-span design warnings,
  and source-linked facts.
- `binding_kinetics` opts into `steady_state.enabled=true`: declared eligible
  plateau windows with shared kinetic Rmax, existing hyperbola, KD_ss/KD_kin diagnostic. Default kinetics
  scientific artifacts are unchanged.
- All new published worked-example gates UNMET; five registered calibration
  criteria failed. Known-Pt coverage 93.1% at Pt/KD=1 and 80.0% at 100; fitted-Pt
  titration diagnostic withholding 48.5% at 10 and 91.5% at 100; unweighted cell
  coverage 93.4% (bound 93.62%). Keep all failures and distinguish any point
  withholding (95.7% at fitted ratio 100) from the registered diagnostic metric.
  Relative weighting/quadratic cell coverage and standalone steady-state coverage
  are uncalibrated. Four-state competition is unavailable.

## 0.11.1 routine statistics

| Purpose | Skill | State |
|---|---|---|
| Independent crossed two-factor designs | group-comparison | Experimental code; gates incomplete: SS II/III, interaction-first EMMs, Tukey/Dunnett/Sidak/Holm; shared residual variance, no repeated units |
| Independent-subject categorical outcomes | group-comparison | Fisher 2x2 conditional OR, seeded RxC Monte Carlo, chi-square, Newcombe RD/Katz RR, trend, McNemar, Wilson/Clopper-Pearson |
| Association, linear prediction and method agreement | [correlation-regression](../../correlation-regression/SKILL.md) | Pearson/Spearman, OLS/WLS, declared-ratio Deming jackknife, Passing-Bablok with CUSUM, exact normal LoA quantile intervals |

Facts available for all new types and retrofitted two/multi-group runs.
Calibration misses retained; reviewer live-agent run 6/6 passed; published evidence
is scoped to the public dataset actually reproduced. See
[release record](../../../validation/RELEASE_0.11.1.md).

## 0.12.0 surface kinetic mechanisms

Code is available for local evaluation; the five-gate validation definition of
done is not met. Numerical/public example limits and coverage misses prevent a fully validated
capability claim; advanced-model intervals are uncalibrated at the default 200
bootstrap replicates (deferred). Reviewer live-agent run 5/5 passed.

Opt-in drift, heterogeneous-ligand, bivalent-analyte and mass-transport models
use `advanced_kinetics.py` and the shared `binding_ode.py`. Default 1:1 scientific
outputs remain unchanged. Advanced fits require predeclared mechanisms,
readout/valency/response-scale declarations, profile support, residual and
window checks, and segment-wise bootstrap. Calibration misses remain visible;
do not claim uniform 95% coverage. Bivalent binding has no single KD; surface
heterogeneity does not establish two epitopes. Off-rate screening reports only
apparent koff. `surface_import.py` supports the pinned T200 XY text and Carterra
XY workbook layouts with explicit metadata and source hashes. Other layouts
remain unsupported. Kinetics now saves extraction-only interpretation facts.
See [0.12.0 evidence](../../../validation/RELEASE_0.12.0.md). Reviewer live-agent run 5/5 passed;
PK/PD remains deferred.
