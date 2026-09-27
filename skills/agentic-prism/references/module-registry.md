# Module registry — development 0.9.3

| Purpose | Skill | Capability state |
|---|---|---|
| Single-site equilibrium binding, KD and diagnostics | [equilibrium-binding](../../equilibrium-binding/SKILL.md) | Implemented; see release validation record |
| Re-render an existing equilibrium run | Same specialist, render command | Implemented; no refitting |
| BLI/SPR 1:1 kon/koff and kinetic KD | [binding-kinetics](../../binding-kinetics/SKILL.md) | Implemented for independent cycles and (0.7.0) single-cycle series in canonical CSV, explicit reference-channel or double referencing, and the scoped Octet Results.txt layout; reliability gates with audit-only intervals for blocked fits. No Biacore/native parser, mass-transport or heterogeneous models |
| ELISA 4PL/5PL calibration and unknown concentration interpolation | [elisa-quantification](../../elisa-quantification/SKILL.md) | Implemented for per-plate standards with back-calculation QC, independent-control gates, descriptive cross-plate QC, range-checked inverse and dilution, (0.7.0) predeclared 5PL, calibration-conditional delta-method unknown intervals and within-plate dilution linearity. No mixed-model intermediate precision or matrix model |
| Relative EC50/IC50 from symmetric 4PL curves; replicate summaries; reference-vs-test relative potency | [dose-response](../../dose-response/SKILL.md) | Implemented for long-form CSV with explicit endpoint/direction/units; optional fixed plateaus and relative weighting; parallel-line RP with F-test or (0.7.0) predeclared-margin equivalence parallelism and optional RP acceptance limits |
| Two-group independent or paired comparisons | [group-comparison](../../group-comparison/SKILL.md) | Implemented: Welch or paired two-sided t test on independent-unit data |
| One-way comparison of ≥3 independent groups | [group-comparison](../../group-comparison/SKILL.md) | Implemented (0.7.0): classic or Welch ANOVA with predeclared Dunnett, Tukey-Kramer, Games-Howell or Holm-Welch family. Parametric independent-group contract; rank tests use the separate entry below |
| Repeated measurements: one within-unit factor, or between-unit arm × within-unit condition | [repeated-measures](../../repeated-measures/SKILL.md) | Implemented (0.8.0): complete RM ANOVA and split-plot ANOVA with GG correction; REML random-intercept models (one-factor and arm × condition) with explicit available-case policy and Satterthwaite (lmerTest-equivalent), parametric bootstrap (one-factor only) or Wald inference. No random slopes (see tumor-growth for continuous-time growth), KR, nested/crossed effects or serial correlation. See 0.8.0 evidence |
| Install / locate / check the runtime | `install.py`, `agentic-prism doctor` ([runtime instructions](runtime.md)) | Implemented (0.7.1): per-clone `.venv` with lock-constrained dependencies and self-check; macOS evidence only |
| Time-to-event (survival, time to humane endpoint) | [time-to-event](../../time-to-event/SKILL.md) | Implemented (0.8.0): Kaplan–Meier with log/log-log limits, medians, landmark survival, log-rank (asymptotic or permutation) with pairwise Holm families, Cox (Efron) hazard ratios with optional covariates, cox.zph PH test; equal to R survival. No competing risks, recurrent events, time-varying covariates, frailty or interval censoring |
| Longitudinal tumor volumes (in-vivo efficacy) | [tumor-growth](../../tumor-growth/SKILL.md) | Implemented (0.8.0): log-volume random-slope mixed model with Satterthwaite (equal to lmerTest) growth rates, rate differences, doubling times and model T/C; observed TGI%/T/C% with Fieller limits and humane-endpoint dropout diagnostics. No Gompertz/logistic curves, cage effects or AR(1) |
| Plate-reader grid + plate map to long table | `agentic-prism import-plate` (used by the ELISA and dose-response Skills) | Implemented (0.7.0) for a generic 96/384-well grid CSV layout; no vendor-specific export parser |

| Rank-based unit-level comparisons | [group-comparison](../../group-comparison/SKILL.md) | 0.8.2: MW/signed-rank with HL and exact untied or normal approximate intervals; KW/Dunn Holm/Bonferroni; complete-block Friedman. See nonparametric contract and 0.8.2 evidence |
| Marginal repeated measurements | [repeated-measures](../../repeated-measures/SKILL.md) | 0.8.2: common US or ordered-visit AR(1), optional numeric baseline covariates; Python REML/Satterthwaite or optional R mmrm KR. Separate `analysis_type=mmrm` contract; interpretation facts included. US small-sample calibration misses are recorded in 0.8.2 evidence |

| ADA cut points | [ada-cut-point](../../ada-cut-point/SKILL.md) | 0.9.0 scoped implementation: complete balanced repeated negative panels, screening/confirmatory/titer point percentiles, normalization, prespecified transformations/outliers, subject/run ANOVA components. Fixed/floating with gates; dynamic recommendation with withheld deployment. Current calibration is superseded by 0.9.1 exact conditional FPR; titer miss disclosed. No sensitivity, drug tolerance or full method validation |

| ADA lower bound, sensitivity, drug tolerance | [ada-cut-point](../../ada-cut-point/SKILL.md) | 0.9.3: two-way (subject × run) bootstrap nonparametric lower bound with a tail-support gate; `ada_sensitivity` (per-run log-linear crossing, log-normal prediction limit) and `ada_drug_tolerance` (per PC level and run, censored at the edges). See 0.9.3 evidence |
| Bioanalytical method validation (LBA) | [method-validation](../../method-validation/SKILL.md) | 0.9.3: accuracy/precision with total error, dilution linearity and hook effect, parallelism with trend, selectivity, specificity, stability from back-calculated concentrations against user-declared criteria; bias, MLS and beta-expectation supplements. No calibration-curve fitting, carry-over or incurred-sample reanalysis |
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
