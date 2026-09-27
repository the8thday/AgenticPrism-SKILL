# AgenticPrism — modular scientific analysis skills

**English** | [简体中文](README.zh-CN.md)

Agent Skills that reproduce GraphPad Prism-style analyses common in antibody
drug development, backed by one versioned Python package for fitting, statistics
and reporting. An AI agent reads a Skill, checks whether the analysis applies to
the experiment, writes an explicit JSON config and runs the command-line tool.
Every run saves its inputs, configuration, hashes, results, diagnostics and a
self-contained offline HTML report.

Local development version **0.7.1**.

> **AI agents:** to install this collection for a user, follow
> [INSTALL_FOR_AGENTS.md](INSTALL_FOR_AGENTS.md).

| Analysis | Skill | What it reports |
|---|---|---|
| Equilibrium binding | [equilibrium-binding](skills/equilibrium-binding/SKILL.md) | Single-site KD with profile-F interval; summary across independent experiments |
| BLI / SPR kinetics | [binding-kinetics](skills/binding-kinetics/SKILL.md) | Global 1:1 kon, koff and kinetic KD from independent cycles or single-cycle series; explicit reference or double referencing; block-bootstrap intervals; Octet `Results.txt` import |
| Dose response | [dose-response](skills/dose-response/SKILL.md) | Relative EC50/IC50 from a symmetric 4PL (optional fixed plateaus, relative 1/Y² weighting); summary across experiments; reference-vs-test relative potency with parallelism by F test or by predeclared equivalence margins, optional RP acceptance limits |
| ELISA quantification | [elisa-quantification](skills/elisa-quantification/SKILL.md) | Plate-specific 4PL or predeclared 5PL standard curve; standard back-calculation QC, independent controls and per-plate screening bounds; unknown interpolation within that range with dilution correction, calibration-conditional intervals and dilution-linearity checks |
| Group comparison | [group-comparison](skills/group-comparison/SKILL.md) | Two groups: Welch or paired t test. Three or more independent groups: classic or Welch one-way ANOVA with a predeclared Dunnett, Tukey-Kramer, Games-Howell or Holm-Welch family |
| Router | [agentic-prism](skills/agentic-prism/SKILL.md) | Picks the specialist by experimental purpose |

Plate-reader grids plus a plate map convert to the long input tables with
`agentic-prism import-plate`. Repeated-measures ANOVA, mixed-effects models and
survival analysis are not implemented. See the
[module registry](skills/agentic-prism/references/module-registry.md).

**No numerical equivalence to GraphPad Prism is claimed.** What has actually
been checked, and its limits, is recorded in [validation/README.md](validation/README.md).

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
The six Skills and the package are one version and must stay together.

Platform evidence: install and results were checked on macOS (uv with Python
3.13, and pip with Python 3.14; identical results). Linux and Windows are
expected to work but the CI workflow has not yet run. See [validation/RELEASE_0.7.1.md](validation/RELEASE_0.7.1.md).

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
```

Each analysis needs a new output directory; earlier results are never
overwritten. `render` only restyles figures: it verifies the result hashes and
never refits. Reports are single self-contained HTML files with all figures and
downloads embedded, so they can be copied on their own.

Start new data from a fixture config, and fill in the applicability evidence
according to the input contract:

- Equilibrium: `fixtures/synthetic/config.json`, [input contract](skills/equilibrium-binding/references/input-schema.md)
- Kinetics: `fixtures/kinetics_synthetic/config.json`, [input contract](skills/binding-kinetics/references/input-and-model.md)
- 4PL: `fixtures/dose_synthetic/config_*.json`, [input contract](skills/dose-response/references/input-and-model.md)
- ELISA: `fixtures/elisa_synthetic/config.json`, [input contract](skills/elisa-quantification/references/input-and-model.md)
- Group statistics: `fixtures/groups_synthetic/*_config.json`, `fixtures/multigroup_synthetic/*.json`, [input contract](skills/group-comparison/references/input-and-model.md)
- Plate-reader grids: `fixtures/plate_import_example/plate_manifest.json`, [grid format](skills/elisa-quantification/references/input-and-model.md#plate-reader-grids-070)

Do not copy the `true` applicability flags from simulated data as if they had
been verified for a real experiment.

## Example reports

Every command in [Usage](#usage) writes a new run directory with a single
self-contained `report.html` (figures and downloads embedded) that opens
offline in a browser. Precomputed example reports are not included in this
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
  unit-level observations. Repeated measures across more than two conditions,
  mixed effects, covariates and nonparametric tests are not included. Welch
  procedures with a group below six units are flagged as possibly liberal.
- **Data provenance:** public fixtures keep the authors' data and source
  checksums, and reproduce specific published settings. They are not general
  defaults for antibody experiments and do not establish independence between
  experiments. Where a separate data license could not be confirmed, none is
  asserted over third-party data.
