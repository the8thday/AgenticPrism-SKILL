# AgenticPrism — modular scientific analysis skills

**English** | [简体中文](README.zh-CN.md) | [Changelog](CHANGELOG.md)

Agent Skills for the statistical analyses of antibody and biologics research
and development, from binding characterization to in vivo efficacy,
bioanalytical validation and immunogenicity. An AI agent reads a Skill, checks
that the analysis fits the experiment, writes an explicit JSON config and runs
one versioned Python package. Every run saves its inputs, configuration,
hashes, results, diagnostics, an `interpretation_facts.json` for the agent's
narrative (in most modules) and an offline HTML report.

Local development version **0.9.3**.

> **AI agents:** to install this collection for a user, follow
> [INSTALL_FOR_AGENTS.md](INSTALL_FOR_AGENTS.md).

## What it covers

| Stage | Question | Skill | Methods |
|---|---|---|---|
| Binding characterization | Equilibrium affinity (KD) | [equilibrium-binding](skills/equilibrium-binding/SKILL.md) | Single-site model with profile-F interval; summary over independent experiments |
| | Association and dissociation rates | [binding-kinetics](skills/binding-kinetics/SKILL.md) | Global 1:1 fit of independent-cycle or single-cycle BLI/SPR data, reference or double referencing, block-bootstrap intervals with reliability gates, Octet `Results.txt` import |
| In vitro function | EC50/IC50 and relative potency | [dose-response](skills/dose-response/SKILL.md) | Symmetric 4PL, relative midpoint with profile interval, parallel-line relative potency with F-test or equivalence-margin parallelism |
| | Concentration from a standard curve | [elisa-quantification](skills/elisa-quantification/SKILL.md) | Per-plate 4PL/5PL, back-calculation QC, independent QC gates, delta-method unknown intervals, dilution linearity; plate-grid import |
| In vivo efficacy | Tumor volumes over time | [tumor-growth](skills/tumor-growth/SKILL.md) | Log-volume random-slope model: growth rates, doubling times, rate differences, model T/C; observed TGI% and T/C% with Fieller limits and dropout diagnostics |
| | Survival, time to humane endpoint | [time-to-event](skills/time-to-event/SKILL.md) | Kaplan–Meier, log-rank (asymptotic or permutation), Cox hazard ratios, proportional-hazards test |
| | Body weight or other scheduled measurements | [repeated-measures](skills/repeated-measures/SKILL.md) | RM and split-plot ANOVA with GG correction, random-intercept models with Satterthwaite, marginal US/AR(1) MMRM with optional Kenward–Roger |
| Bioanalysis and immunogenicity | Method validation (ICH M10-style, ligand-binding assays) | [method-validation](skills/method-validation/SKILL.md) | Accuracy/precision with total error, dilution linearity and hook effect, parallelism with trend check, selectivity, specificity, stability; accuracy profile |
| | ADA cut points, sensitivity, drug tolerance | [ada-cut-point](skills/ada-cut-point/SKILL.md) | Screening/confirmatory/titer cut points, fixed or floating, confidence lower bounds; positive-control sensitivity with a run-to-run prediction limit; drug tolerance |
| | Repeatability and intermediate precision | [variance-components](skills/variance-components/SKILL.md) | REML for nested/crossed factors, unbalanced data, MLS/MOVER intervals |
| Any stage | Comparing groups | [group-comparison](skills/group-comparison/SKILL.md) | Welch or paired t, one-way ANOVA with Dunnett/Tukey/Games-Howell/Holm families, Mann–Whitney and signed-rank with Hodges–Lehmann, Kruskal–Wallis with Dunn, Friedman |

The [router Skill](skills/agentic-prism/SKILL.md) picks the specialist from the
experiment; the [module registry](skills/agentic-prism/references/module-registry.md)
lists exact capabilities and boundaries.

## How correctness is checked

Each method is compared with an established implementation on public or
fixture data, and its intervals and error rates are checked by simulation
with pass/fail bounds fixed before running. Misses are reported as misses.
Details, scope and every number: [validation/README.md](validation/README.md).

| Module | Compared with | Largest difference | Calibration by simulation |
|---|---|---|---|
| Equilibrium binding | BindCurve 0.2.0 (25 curves) | KD 1.3e-4 relative | CI coverage check |
| Kinetics | Published TitrationAnalysis fits of Octet RED384 data | rates/KD 2e-4 relative | Bootstrap coverage; single-cycle koff 91.3% vs 91.4% bound (miss) |
| Dose response | Independent SciPy four-parameter refit | 5.8e-7 relative | 4PL and relative-potency coverage |
| ELISA | Independent full-parameter fit with root finding | 3.9e-10 | Unknown-interval coverage |
| Group comparison | SciPy; R `wilcox.test`, `kruskal.test`, `friedman.test`, `dunn.test` | F 1e-10 relative; rank statistics exact or 1e-14; HL interval endpoints within 3e-4 absolute (root tolerance) | Family-wise error; Welch procedures miss with a 4-unit high-variance group |
| Repeated measures, MMRM | R nlme, lmerTest, afex, emmeans, mmrm, pbkrtest | 1.8e-6 relative (df, MMRM) | Satterthwaite passes; unstructured MMRM with 12 units per arm misses (7.0–7.1%) |
| Time to event | R survival (lung, veteran) | 1.6e-14 relative | Permutation log-rank, KM log-log, Cox and PH test pass |
| Tumor growth | R lmerTest | 1.9e-7 relative | Rate-difference coverage under dropout passes |
| ADA cut points | R base/car/lme4; rADA vignette | 2.1e-12 relative | Lower bounds pass; point titer cut point misses |
| Precision components | R VCA and lme4; CLSI EP05 example | 9.3e-7 absolute | MLS total passes; one component upper bound misses (92.8%) |
| Method validation, ADA sensitivity | R VCA `anovaVCA`, `lm`, `binom.test`, `t.test`, `approx` | 3.4e-12 relative | 17/18 pass; sensitivity limit with 3-fold dilutions misses |

Live-agent misuse scenarios (for example asking for a verdict without
acceptance criteria, omitting a hook effect, or quoting positive-control
sensitivity as a patient detection limit) are recorded per release in
[validation/agent-scenarios](validation/agent-scenarios/README.md).
**No numerical equivalence to GraphPad Prism is claimed.**

## Installation

Requirements: git, internet access for the first install, and either
[uv](https://docs.astral.sh/uv/) (recommended; it fetches the validated Python
3.13 by itself) or Python 3.12+.

```sh
git clone https://github.com/the8thday/AgenticPrism-SKILL.git AgenticPrism
cd AgenticPrism
python3 install.py          # Windows: py install.py
```

`install.py` creates `.venv` inside the clone, installs the package with the
dependency versions from `requirements-lock.txt`, and runs a self-check (an
analysis, a render and a hash verification). Rerun it after `git pull`; the
Skills check that the runtime and Skill versions match. Remove `.venv` to
uninstall.

**Using the Skills with an agent.** Point the agent at the router Skill:

> Use `<clone>/skills/agentic-prism/SKILL.md` to analyze my experiment files.

For agents that discover Skills from a folder (for example Claude Code's
`~/.claude/skills`), link them there:

```sh
python3 install.py --link-skills ~/.claude/skills
```

This creates symlinks and never overwrites existing entries. Do not copy Skill
folders: a copied Skill cannot find its runtime. Each Skill locates the runtime
from its real location and checks it with `agentic-prism doctor --collection
<clone>` ([runtime instructions](skills/agentic-prism/references/runtime.md)).
The twelve Skills and the package are one version and must stay together.

Platform evidence: install and results were checked on macOS (uv with Python
3.13, and pip with Python 3.14; identical results). Linux and Windows are
expected to work but the CI workflow has not yet run. See [validation/RELEASE_0.7.1.md](validation/RELEASE_0.7.1.md).

**Optional R.** Only MMRM Kenward–Roger needs R. After installing R, opt in with
`python3 install.py --with-r`; pinned packages live in `.r-lib/`, and `doctor`
reports them. Without R, only that method refuses; everything else runs.

## Usage

```sh
.venv/bin/agentic-prism analyze --config fixtures/public/config.json --output runs/my-public-run
.venv/bin/agentic-prism render --run runs/my-public-run --style standard
.venv/bin/agentic-prism verify --run runs/my-public-run
.venv/bin/agentic-prism analyze --config fixtures/kinetics_octet/config_300s.json --output runs/my-octet-300s
.venv/bin/agentic-prism import-octet --manifest fixtures/kinetics_octet/octet_manifest.json --output runs/my-octet-import
.venv/bin/agentic-prism analyze --config fixtures/dose_synthetic/config_ic50.json --output runs/my-ic50-run
.venv/bin/agentic-prism analyze --config fixtures/dose_synthetic/config_potency.json --output runs/my-potency-run
.venv/bin/agentic-prism analyze --config fixtures/elisa_synthetic/config.json --output runs/my-elisa-run
.venv/bin/agentic-prism analyze --config fixtures/groups_synthetic/paired_config.json --output runs/my-paired-run
.venv/bin/agentic-prism analyze --config fixtures/multigroup_synthetic/config_dunnett.json --output runs/my-multigroup-run
.venv/bin/agentic-prism analyze --config fixtures/elisa_dilution_5pl/config.json --output runs/my-5pl-run
.venv/bin/agentic-prism analyze --config fixtures/kinetics_single_cycle/config.json --output runs/my-sck-run
.venv/bin/agentic-prism import-plate --manifest fixtures/plate_import_example/plate_manifest.json --output runs/my-plate-import
.venv/bin/agentic-prism analyze --config fixtures/method_validation/config_accuracy_precision.json --output runs/my-ap-run
.venv/bin/agentic-prism analyze --config fixtures/ada_performance/config_sensitivity.json --output runs/my-ada-sensitivity
```

Each analysis needs a new output directory; earlier results are never
overwritten. `render` only restyles figures: it verifies the result hashes and
never refits. Curve-fitting reports embed figures and downloads in one HTML file. ADA and
precision reports display results offline but link to adjacent artifact files;
copy their complete run directory to preserve downloads.

Start new data from a fixture config, and fill in the applicability evidence
according to the input contract:

- Equilibrium: `fixtures/synthetic/config.json`, [input contract](skills/equilibrium-binding/references/input-schema.md)
- Kinetics: `fixtures/kinetics_synthetic/config.json`, [input contract](skills/binding-kinetics/references/input-and-model.md)
- 4PL: `fixtures/dose_synthetic/config_*.json`, [input contract](skills/dose-response/references/input-and-model.md)
- ELISA: `fixtures/elisa_synthetic/config.json`, [input contract](skills/elisa-quantification/references/input-and-model.md)
- Group statistics: `fixtures/groups_synthetic/*_config.json`, `fixtures/multigroup_synthetic/*.json`, [input contract](skills/group-comparison/references/input-and-model.md)
- Plate-reader grids: `fixtures/plate_import_example/plate_manifest.json`, [grid format](skills/elisa-quantification/references/input-and-model.md#plate-reader-grids-070)
- Method validation: `fixtures/method_validation/config_*.json`, [input contract](skills/method-validation/references/contract.md)
- ADA cut points, sensitivity, drug tolerance: `fixtures/ada_synthetic/`, `fixtures/ada_performance/`, [input contract](skills/ada-cut-point/references/contract.md)

Do not copy the `true` applicability flags from simulated data as if they had
been verified for a real experiment.

## Example reports

Every command in [Usage](#usage) writes a new run directory with `report.html`
that opens offline in a browser. Keep the run directory for complete artifacts. Precomputed example reports are not included in this
repository.

## Outputs and exit codes

A run contains an input snapshot, the resolved configuration, input/code/
environment hashes, normalized observations, results, diagnostics,
SVG/PDF/PNG figures in two themes, and the HTML report. The curve-fitting modules
also save predictions, residuals and parameter intervals where supported.
The equilibrium and dose-response modules add per-sample
summaries; kinetics saves per-sensor parameters and joint bootstrap samples;
dose-response comparisons add relative-potency tables. ELISA adds standard-QC, unknown-well (with intervals), sample and dilution-linearity tables;
two-group comparisons add a mean-difference table; multi-group comparisons add
omnibus, group-summary and contrast tables; equivalence parallelism adds
`parallelism_equivalence.csv`. Outputs separate audit
values from reportable estimates. Out-of-range values are never turned into
upper or lower bounds that imply statistical guarantees.

- Exit 0: workflow completed with no hard curve failures. Some results may still be limited or out of range.
- Exit 3: some curves failed; results and report are saved for inspection. The simulated equilibrium example contains a flat curve on purpose.
- Exit 2: configuration, input or execution error; a started run directory contains `failure.json`.

## Scope and limitations

- **Equilibrium:** intervals are computed for KD only. There are no baseline or
  amplitude intervals, curve confidence bands or prediction intervals. A
  baseline estimated from controls and then fixed does not propagate its
  uncertainty, and the report says so.
- **Kinetics:** the public example is an author-processed Octet RED384
  `Results.txt` layout. Reference results come from TitrationAnalysis, not the
  Octet vendor software. The 300 s and 600 s runs are two fit windows of one
  experiment. `.frd` files are not read. SPR data can use the canonical CSV,
  but there is no dedicated Biacore importer. Single-cycle curves declare each
  row's own injection start, end and concentration; double referencing
  `(sample − reference) − (blank − blank reference)` runs only when configured.
- **Dose response:** the 4PL fits raw per-well responses and reports the
  **relative** EC50/IC50, halfway between the fitted (or fixed) plateaus, with
  a profile-F interval. It equals the absolute Y=50 concentration only when
  the plateaus are fixed at 0 and 100, and it is not a KD or Ki. Relative
  potency RP = C50_reference / C50_test comes from a parallel model with shared
  plateaus and Hill slope. The default parallelism F test is not the equivalence test
  USP <1032>/<1034> recommend; `parallelism_method=equivalence` requires the
  Hill-ratio (and optional plateau-difference) intervals to lie inside
  predeclared margins, and `rp_acceptance_limits` adds an RP pass/fail check.
  Margins and limits must come from the lab's historical data or protocol. The public example has no physical dose unit
  and only demonstrates the numerics and workflow. Real assays still need
  explicit dose units, readout meaning, replicate structure and normalization
  history.
- **ELISA:** per-plate 4PL (or predeclared 5PL) standards. Each standard level is back-calculated;
  by default a level passes at 80–120% mean recovery (75–125% at the LLOQ and
  ULOQ) and replicate CV ≤ 20%, and a plate needs at least 75% of levels and
  at least six passing (all limits configurable). These are software screening
  rules, not full ICH M10 validation; LLOQ/ULOQ are legacy keys for the current
  plate screening bounds. Independently prepared QC can additionally gate unknowns. Failing
  standards are reported, not excluded automatically. Unknown intervals (delta
  method) reflect only this plate's calibration and well noise, not plate,
  matrix or dilution error. Dilution linearity is checked within a plate; the
  hook/matrix flag is a pattern, not a diagnosis.
- **Group statistics:** one predeclared Welch or paired t contrast, or a one-way
  ANOVA over ≥3 independent groups with a predeclared contrast family, on
  unit-level observations. Repeated measures across more than two conditions and scoped random-intercept
  mixed models use the separate repeated-measures specialist. Rank tests use the separate nonparametric contract; baseline numeric covariates
  are supported only in the MMRM contract. Welch
  procedures with a group below six units are flagged as possibly liberal.
- **Data provenance:** public fixtures keep the authors' data and source
  checksums, and reproduce specific published settings. They are not general
  defaults for antibody experiments and do not establish independence between
  experiments. Where a separate data license could not be confirmed, none is
  asserted over third-party data.
- **In vivo and longitudinal:** repeated measures use compound-symmetry random
  intercepts, split-plot ANOVA or marginal US/AR(1) MMRM; tumor growth is
  exponential on the log scale with MAR removal; time-to-event has no competing
  risks, frailty or interval censoring. Small-sample misses are listed in the
  validation records.
- **Bioanalysis and immunogenicity:** method validation starts from
  back-calculated concentrations (no calibration-curve fitting, carry-over or
  incurred-sample reanalysis); acceptance criteria are always user-declared.
  ADA cut points need complete balanced negative panels from one reagent lot;
  dynamic cut points are not implemented. Positive-control sensitivity is not
  the sensitivity for patients' antibodies.
