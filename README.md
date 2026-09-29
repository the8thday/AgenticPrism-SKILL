# AgenticPrism — modular scientific analysis skills

**English** | [简体中文](README.zh-CN.md) | [Changelog](CHANGELOG.md)

Agent Skills for the statistical analyses of antibody and biologics research
and development, from binding characterization to in vivo efficacy,
bioanalytical validation, immunogenicity and CMC quality (stability,
potency, comparability, specifications). An AI agent reads a Skill, checks
that the analysis fits the experiment, writes an explicit JSON config and runs
one versioned Python package. Every run saves its inputs, configuration,
hashes, results, diagnostics and an offline HTML report. Every module also
writes `interpretation_facts.json`: the saved results the agent must base its
narrative on (primary estimates, what is reportable or withheld and why,
required disclosures, limitations). Earlier saved runs are not rewritten.

## Status: public beta 0.13.0

AgenticPrism is for research analysis. It has **not** been validated as a
computerized system for GLP/GMP work or regulatory submissions (for example
21 CFR Part 11), and it is not a substitute for a statistician's review of
decisions that depend on the result. Installation and results have been checked
on macOS only. No numerical comparison with GraphPad Prism has been made, so no
Prism equivalence is claimed. All validation so far uses simulated and public
data; feedback from real project data is welcome through GitHub issues.

| Maturity | Modules |
|---|---|
| Validated against reference implementations and calibrated | Equilibrium KD, 1:1 kinetics, 4PL dose response and relative potency, ELISA, two-group and multi-group comparisons, rank tests, repeated measures and MMRM, time-to-event, tumor growth, ADA cut points, variance components, method validation, stability, potency assay, comparability, specifications, routine statistics (0.11.1) |
| Usable with stated limits | Affinity depth and cell binding (no published worked example), advanced surface kinetics (bootstrap intervals not calibrated at default settings), epitope binning (no published worked example), drug combination (Loewe/ZIP not identical to synergyfinder) |
| Exploratory | HTS hit calling (false-hit rate above nominal under the global null) |

Each module's exact limits are in the
[module registry](skills/agentic-prism/references/module-registry.md) and the
linked release records; every calibration miss is kept in `validation/`.

## 0.13.0

HTS nominal BH screening has not demonstrated FDR control: the registered global-null false-hit rate was 0.100 (bound 0.063784); review traced it to median-polish residuals (4.3% without polishing, 9.3% with). HSA interval coverage also missed; see the retained evidence.

Named Bliss, Loewe, HSA and ZIP combination references with independent-matrix uncertainty; declared HTS plate QC, median-polish B scores and exploratory FDR hits. Dose-response and ELISA facts complete the oldest-six retrofit.

See [release evidence](validation/RELEASE_0.13.0.md) for retained misses and unmet gates. Live-agent review: 4/4 scenarios passed. PK/PD remains deferred.

## 0.12.1

Directed epitope binning with declared controls and thresholds, asymmetric-pair diagnostics, conditional bootstrap clustering stability and reciprocal-block communities.

See [release evidence](validation/RELEASE_0.12.1.md) for retained misses and unmet gates. Live-agent review: 2/3 passed; the by-eye bin-merge scenario failed in 2 of 3 attempts, when the agent answered without reading the Skill. PK/PD remains deferred.

> **AI agents:** to install this collection for a user, follow
> [INSTALL_FOR_AGENTS.md](INSTALL_FOR_AGENTS.md).

## What it covers

| Stage | Question | Skill | Methods |
|---|---|---|---|
| Pharmacology and screening | Combination reference scores | [drug-combination](skills/drug-combination/SKILL.md) | Bliss/Loewe/HSA/ZIP with independent-matrix intervals |
| | HTS quality and exploratory hits | [import-plate](skills/import-plate/SKILL.md) | Z-prime, SSMD, median-polish B scores and FDR |
| Binding characterization | Epitope competition and bins | [epitope-binning](skills/epitope-binning/SKILL.md) | Directed controlled blocking, conditional bootstrap stability and reciprocal-block graphs |
| Binding characterization | Equilibrium affinity (KD) | [equilibrium-binding](skills/equilibrium-binding/SKILL.md) | Hyperbola or opt-in exact depletion, SET n-curves, three-state competition Ki; profile-F and independent-experiment summaries |
| | Cell-surface apparent affinity | [cell-binding](skills/cell-binding/SKILL.md) | Joint total/control fit, declared background, hyperbolic or quadratic receptor depletion; apparent KD only |
| | Association and dissociation rates | [binding-kinetics](skills/binding-kinetics/SKILL.md) | Global 1:1 fit of independent-cycle or single-cycle BLI/SPR data, reference or double referencing, block-bootstrap intervals with reliability gates; opt-in drift, heterogeneous ligand, bivalent analyte, mass transport and off-rate screening; scoped Octet/T200/Carterra imports |
| In vitro function | EC50/IC50 and relative potency | [dose-response](skills/dose-response/SKILL.md) | Symmetric 4PL, relative midpoint with profile interval, parallel-line relative potency with F-test or equivalence-margin parallelism |
| | Concentration from a standard curve | [elisa-quantification](skills/elisa-quantification/SKILL.md) | Per-plate 4PL/5PL, back-calculation QC, independent QC gates, delta-method unknown intervals, dilution linearity; plate-grid import |
| In vivo efficacy | Tumor volumes over time | [tumor-growth](skills/tumor-growth/SKILL.md) | Log-volume random-slope model: growth rates, doubling times, rate differences, model T/C; observed TGI% and T/C% with Fieller limits and dropout diagnostics |
| | Survival, time to humane endpoint | [time-to-event](skills/time-to-event/SKILL.md) | Kaplan–Meier, log-rank (asymptotic or permutation), Cox hazard ratios, proportional-hazards test |
| | Body weight or other scheduled measurements | [repeated-measures](skills/repeated-measures/SKILL.md) | RM and split-plot ANOVA with GG correction, random-intercept models with Satterthwaite, marginal US/AR(1) MMRM with optional Kenward–Roger |
| Bioanalysis and immunogenicity | Method validation (ICH M10-style, ligand-binding assays) | [method-validation](skills/method-validation/SKILL.md) | Accuracy/precision with total error, dilution linearity and hook effect, parallelism with trend check, selectivity, specificity, stability; accuracy profile |
| | ADA cut points, sensitivity, drug tolerance | [ada-cut-point](skills/ada-cut-point/SKILL.md) | Screening/confirmatory/titer cut points, fixed or floating, confidence lower bounds; positive-control sensitivity with a run-to-run prediction limit; drug tolerance |
| | Repeatability and intermediate precision | [variance-components](skills/variance-components/SKILL.md) | REML for nested/crossed factors, unbalanced data, MLS/MOVER intervals |
| CMC and quality | Long-term stability and shelf life | [stability](skills/stability/SKILL.md) | Q1E linear regression, slope-first poolability, mean confidence bounds and declaration-gated extrapolation |
| | Relative potency across runs and validation | [potency-assay](skills/potency-assay/SKILL.md) | Log-RP random-run REML, MLS/MOVER intermediate precision, bias, linearity equivalence and tested range; failing runs retained |
| | Lot comparability and biosimilarity (one attribute) | [comparability](skills/comparability/SKILL.md) | TOST equivalence of lot means against a declared margin, quality range with declared k, or descriptive; tier recorded, never inferred |
| | Tolerance intervals and process capability | [specifications](skills/specifications/SKILL.md) | Exact normal and nonparametric tolerance intervals (minimum n when too few lots); Pp/Ppk and within-subgroup Cp/Cpk with intervals |
| Any stage | Comparing groups | [group-comparison](skills/group-comparison/SKILL.md) | Welch or paired t, one-way ANOVA with Dunnett/Tukey/Games-Howell/Holm families, Mann–Whitney and signed-rank with Hodges–Lehmann, Kruskal–Wallis with Dunn, Friedman |

The [router Skill](skills/agentic-prism/SKILL.md) picks the specialist from the
experiment; the [module registry](skills/agentic-prism/references/module-registry.md)
lists exact capabilities and boundaries.

## 0.12.0 surface kinetic depth

Opt-in drift, heterogeneous ligand, bivalent analyte and mass transport;
dissociation-only koff ranking; verified T200/Carterra XY imports. Complex
mechanisms require predeclaration and identifiability. Default 1:1 stays unchanged.
Kinetics now saves interpretation facts. Calibration misses and unmet public
worked-example gates are retained; live-agent review 5/5 passed.
**Advanced-model bootstrap intervals are not calibrated at the default 200 replicates.** The registered calibration used 50 replicates and covered 87–92%; with fast mass transport, fits that passed every gate covered only 69–77%. Treat advanced-model intervals as provisional; a default-setting calibration is deferred.
See [release evidence](validation/RELEASE_0.12.0.md). PK/PD remains deferred.

## 0.11.1 routine statistics

Adds independent two-way ANOVA, contingency tables/proportions, correlation and
linear regression, and Deming/Passing-Bablok/Bland-Altman method comparison.
The new `correlation-regression` Skill distinguishes association from agreement.
Two/multi-group runs now save interpretation facts without changing prior
scientific artifacts. SS type, contrast family, independent units and Deming
error ratio must be declared. Numerical agreement is scoped; calibration misses
are retained; live-agent review 6/6 passed.
See [release evidence](validation/RELEASE_0.11.1.md).

## 0.11.0 affinity depth

New models are opt-in. Active Pt needs a source and activity basis; fitted Pt is
profiled, with points withheld in the fitted titration regime or lower-open
profile. Cell binding has its own Skill and always reports apparent KD, with
measured nonspecific controls and explicit wash/detection/valency gates.
Kinetics can additionally extract Req from declared eligible plateau windows;
KD_ss/KD_kin is diagnostic and is never averaged.

Evidence is limited: all named new published worked-example gates are **unmet**.
The registered 14,000 simulations retained five failed criteria: known-Pt
coverage 93.1% (Pt/KD=1) and 80.0% (100), fitted-Pt titration diagnostic
withholding 48.5% (10) and 91.5% (100), and cell apparent-KD coverage 93.4%,
against 93.62%. At fitted Pt/KD=100, any point withholding was 95.7%; that does
not replace the preregistered diagnostic-specific miss. SET and Ki coverage were
95.2% and 95.8%. Relative-weighted and quadratic cell fits were not separately
coverage-calibrated. A post-hoc breakdown (not a replacement) shows the 80.0%
row's design titrated only up to Pt: 158 of 1,000 estimates fell above the
titrated range and were all withheld, and the 842 reportable fits covered 95.0%.
Runs now carry a non-blocking warning for that design. SET also requires a
declared valency and readout (bivalent IgG captured as molecules with any free
site is refused), antigen-negative control cells get their own baseline, and
steady state requires a shared kinetic Rmax.
See [the complete release record](validation/RELEASE_0.11.0.md).

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
| CMC stability | R `lm`/`anova`/`emmeans`; Koleva regression example | 3.49e-13 relative | 6/10 rows pass; four coverage misses disclosed; printed shelf life not reproduced |
| Potency across runs | R lme4; independent MLS/MOVER calculation | 6.90e-7 relative | 12/12 rows pass; matching published worked example unavailable |
| Comparability, specifications | R `tolerance` (EXACT factors), `t.test`; NIST/SEMATECH printed tolerance and capability values | 3.3e-9 relative | 13/15 pass; the two n = 10 capability rows miss by Monte Carlo error (100,000-dataset check: 94.9%, 95.4%) |

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
The 21 Skills and the package are one version and must stay together.

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
.venv/bin/agentic-prism analyze --config fixtures/stability/config_pooled.json --output runs/my-shelf-life
.venv/bin/agentic-prism analyze --config fixtures/comparability/config_tost_absolute.json --output runs/my-comparability
.venv/bin/agentic-prism analyze --config fixtures/specification/config_normal_exact.json --output runs/my-tolerance-interval
```

Each analysis needs a new output directory; earlier results are never
overwritten. `render` only restyles figures: it verifies the result hashes and
never refits. Curve-fitting reports embed figures and downloads in one HTML file. ADA and
precision reports display results offline but link to adjacent artifact files;
copy their complete run directory to preserve downloads.

Start new data from a fixture config, and fill in the applicability evidence
according to the input contract:

- Equilibrium: `fixtures/synthetic/config.json`, [input contract](skills/equilibrium-binding/references/input-schema.md)
- Affinity depth: `fixtures/affinity_0110/`, [depletion, SET and competition contract](skills/equilibrium-binding/references/affinity-depth.md)
- Cell binding: `fixtures/affinity_0110/config_cell_*.json`, [input contract](skills/cell-binding/references/input-contract.md)
- Kinetics: `fixtures/kinetics_synthetic/config.json`, [input contract](skills/binding-kinetics/references/input-and-model.md)
- 4PL: `fixtures/dose_synthetic/config_*.json`, [input contract](skills/dose-response/references/input-and-model.md)
- ELISA: `fixtures/elisa_synthetic/config.json`, [input contract](skills/elisa-quantification/references/input-and-model.md)
- Group statistics: `fixtures/groups_synthetic/*_config.json`, `fixtures/multigroup_synthetic/*.json`, [input contract](skills/group-comparison/references/input-and-model.md)
- Plate-reader grids: `fixtures/plate_import_example/plate_manifest.json`, [grid format](skills/elisa-quantification/references/input-and-model.md#plate-reader-grids-070)
- Method validation: `fixtures/method_validation/config_*.json`, [input contract](skills/method-validation/references/contract.md)
- Repeated measures, time to event, tumor growth: `fixtures/repeated_synthetic/`, `fixtures/survival_synthetic/`, `fixtures/tumor_growth_synthetic/`, input contracts for [repeated measures](skills/repeated-measures/references/input-and-model.md), [time to event](skills/time-to-event/references/input-and-model.md) and [tumor growth](skills/tumor-growth/references/input-and-model.md)
- ADA cut points, sensitivity, drug tolerance: `fixtures/ada_synthetic/`, `fixtures/ada_performance/`, [input contract](skills/ada-cut-point/references/contract.md)
- Precision components: `fixtures/variance_reml/`, [input contract](skills/variance-components/references/contract.md)
- CMC: `fixtures/stability/`, `fixtures/potency_assay/`, `fixtures/comparability/`, `fixtures/specification/`, input contracts for [stability](skills/stability/references/input-contract.md), [potency](skills/potency-assay/references/input-contract.md), [comparability](skills/comparability/references/input-contract.md) and [specifications](skills/specifications/references/input-contract.md)

Do not copy the `true` applicability flags from simulated data as if they had
been verified for a real experiment.

## Example reports

Every command in [Usage](#usage) writes a new run directory with `report.html`
that opens offline in a browser. Keep the run directory for complete artifacts. Precomputed example reports are not included in this
repository.

## Outputs and exit codes

A run contains an input snapshot, the resolved configuration, input/code/
environment hashes, normalized observations, results, diagnostics,
SVG/PDF/PNG figures with selectable themes, and the HTML report. The curve-fitting modules
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
  with scoped T200 XY text and Carterra XY workbook imports in 0.12.0. Other layouts remain unsupported. Single-cycle curves declare each
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
- **CMC:** stability uses ICH Q1E linear regression with the α = 0.25
  poolability pre-tests, whose model selection lowers coverage in some designs
  (disclosed); extrapolation needs declared supporting data. Potency
  combination needs replicated determinations per run. Comparability covers one
  attribute at a time and never gives an overall biosimilarity verdict;
  tolerance intervals describe data and do not set specifications. Tier,
  margins, k and specification limits are always user declarations.
