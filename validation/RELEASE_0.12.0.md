# Release evidence — local development 0.12.0

Base dev/0.11.1 at 2e4e991. Local only; PK/PD deferred.

## Phase 0: sources and real export audit (before numerical implementation)

- Myszka et al. 1998, doi:10.1016/S0006-3495(98)77549-6:
  https://pubmed.ncbi.nlm.nih.gov/9675161/ describes simulated transport data.
  No licensed raw experimental series verified. Published mass-transport gate
  UNMET; no synthetic fixture will replace it.
- Karlsson & Falt 1997, doi:10.1016/S0022-1759(96)00195-0:
  publisher abstract accessible, numeric export/supplement and reuse licence
  not verified. Heterogeneous/bivalent candidate gate UNMET for this source.
  https://www.sciencedirect.com/science/article/pii/S0022175996001950
- Found actual T200 XY exported text (1,024,899 bytes), information workbook and
  published TitrationAnalysis parameter table, CC-BY-4.0. Repository commit
  cb8112503658d24c4c6576503243fa506bf99e27. Text has explicit paired _X/_Y
  headers; concentrations and injection times must come from declared metadata,
  never inferred from the curve or mass-concentration header text.
  https://github.com/DukeCHSI/TitrationAnalysis
- Found actual author-downselected Carterra XY workbook and information CSV in
  Nguyen et al.'s bivalent model repository, CC-BY-4.0, commit
  68fefdc9c1ace6e2030017404ac1e2e9c31e0953. It has a title row, ligand/curve
  labels, then X/Y headers. Regenerative workbook is strict OOXML (openpyxl
  does not recognize its namespace); a standard XML reader can read verified
  cells without changing the original workbook. Support will be scoped to this
  observed layout, not every Carterra software version.
  https://github.com/DukeCHSI/Bivalent-Analyte-Modeling
- All source bytes/licences/hashes are in fixtures/kinetics_0120/public/.
  These real samples remove the no-sample blocker for the two scoped layouts.
  No Biacore 8K/Insight-specific or binary .blr/.bme layout is claimed.
- Public bivalent readout: singly and doubly bound analyte contribute one
  analyte mass each, signal X1+X2; surface-site conservation L+X1+2X2=Rmax.
  First-arm association includes the factor two for free analyte arms; second
  association has surface-response inverse units. Rebinding during dissociation
  is a declared choice. No single bivalent KD is defined.

## Implemented scope and deliberate limits

Opt-in models: one_to_one_drift, heterogeneous_ligand, bivalent_analyte,
mass_transport, off_rate_screening. Runtime solve_ivp primitive uses nM/seconds/
response units, with explicit conversion for transport. Advanced global models
currently require regenerated multi-cycle curves and a common response scale
with shared Rmax. Per-curve capacities, single-cycle complex models and automatic
mechanism selection are not supported. Default 1:1 remains the original fitter.
Both imports require explicit selected column pairs, molarity and phase times;
only the two pinned layouts are implemented. No guessed 8K/Insight/binary parser.

Profile endpoint checks are support diagnostics over a declared log span, not
complete profile confidence intervals. Confidence intervals use segment-wise
residual block bootstrap. Sampling points/curves are not independent experiments.
Screening reports apparent koff only; no affinity conversion. Reportable scalar
kinetic KD can be plotted on ka/kd iso-affinity axes. Bivalent has no scalar KD.

## Numerical evidence

- `kinetics_benchmarks_0.12.0.json`: four forward-model comparisons passed.
  deSolve bivalent maximum relative difference 5.604716907002508e-7;
  mass transport 2.1146161534995322e-7. Analytic drift/heterogeneous differences
  below 5e-15. Real source imports preserve observation values and record hashes.
- `kinetics_refits_0.12.0_first_attempt.json` retains failed initial R nlminb
  refits. The default finite-difference step was too small relative to ODE
  tolerance. The independent oracle now uses central 1e-4 log-parameter steps;
  this is a numerical benchmark correction, not calibration retuning.
- Final fixture R refits: heterogeneous max relative difference 5.24096e-5,
  transport 2.92236e-6. Bivalent differs 2.99712e-5 but R convergence code1:
  **UNMET**, not promoted to passing from small parameter differences.
- Public T200 numerical objective checks (first three column pairs, declared
  pre-injection mean subtraction, stride50, injection0/180seconds): drift
  7.48656e-6 and transport7.19177e-7 pass. Heterogeneous2.32704e-2 and
  bivalent2.30597e-2 fail the 1% parameter agreement criterion. Similar SSE
  does not establish identifiable constants. These arbitrary-model checks do
  not assert biological applicability or reproduce author-reported constants.
- The independent refits start from a deterministic perturbation of the Python
  optimum. They validate a local objective/refit, not an independent global search.
- Public author data are usable, but published complex-fit parameter reproduction
  is **UNMET**: the source bivalent analysis permits per-curve capacities and
  different baseline conventions. Current shared-capacity results are not a
  matched worked example. TitrationAnalysis's unchanged 1:1 benchmark was not
  rerun. Anabel was not used where a matching declared mechanism was unverified.

## Calibration registration

`scripts/validate_kinetics_0120_calibration.py` committed as `ab0c931` before
execution: seed120260930, nine rows, 1000 draws each, 50 bootstrap refits each,
three concentrations1/30/1000nM, 600s association/600s dissociation, explicit
short-association stress. Every draw, error and diagnostic is retained.
Coverage is assessed over evaluable and reportable-only runs; zero reportable
runs mean no conditional coverage claim. True1:1 false complex selection is
registered separately; model comparison uses improvement10 plus reliability.
Mass-transport rows vary km=.002/.02/2 s^-1; these record identification and
coverage, not an invented universal transport-limitation threshold.

## Live agents and facts

Five misuse scenarios are **PENDING**, with fixtures/prompts committed. No Claude
CLI or substitute agents ran. Reviewer command (on dev/0.12.0):

```sh
.venv/bin/python scripts/run_agent_scenarios.py /tmp/agentic-review-0120 --release 0.12.0 --attempt reviewer-1
```

New advanced results have saved-result facts with source pointers and reliability
states. Legacy kinetics now extracts facts after all scientific files, before the
manifest; optional steady-state facts are copied without refitting. Failed and
withheld estimates cannot become reportable facts. State-coverage tests cover
reportable, limited and failed states and unchanged source bytes.

## Self-review

- **Valency/avidity:** bivalent model requires bivalent declaration and second
  rebinding choice; monovalent mechanisms refuse bivalent declarations. First-arm
  statistical factor2 and surface units are explicit; no single bivalent KD.
- **Readout:** require bound-analyte mass response. Bivalent signal X1+X2 differs
  from site occupancy X1+2X2; the conservation/zero-rebinding limit is tested.
- **Shared scale:** require declared comparable response scales and shared Rmax.
  This restriction prevents claiming reproduction of the source's per-curve-Rmax
  published analysis. It should be reviewed before extending sharing options.
- **Independence:** regenerated cycles are a design condition, not independent
  biological replicates. Segment-wise bootstrap models trace residuals only.
- **Design range:** minimum counts do not establish identifiability. Profiles,
  correlation, little dissociation, window changes and serial residuals can
  withhold rates despite optimizer convergence. Strong transport or heterogeneous
  rates may remain nonidentifiable. All calibration withholding is retained.
- **Reference/drift:** drift needs a source and valid reference; it cannot rescue
  a known failed reference subtraction.
- **Transport:** conversion from response to surface nM needs a source. Identified
  km does not alone prove transport limitation; no binary biological conclusion
  is inferred from a better AICc.

## Design-range limitation observed during the frozen calibration

The heterogeneous row had 0/1000 reportable fits: every draw triggered
`profile_flat_or_open_kd2`. Its concentration range spans the two true KDs, but
the high-concentration fast association is undersampled by the registered
20-point association grid. Thus the intended adequate-design row does not
establish adequate time resolution. The row, its original label and all misses
are retained; it has not been relabelled after observing results or replaced by
a tuned scenario. This is a deviation from the requested adequate-design
calibration and leaves that validation scope unmet.

The calibration uses 50 bootstrap resamples per simulated dataset, explicitly
registered before execution. It does not calibrate the default200-resample
setting. Low-resample percentile intervals have observed undercoverage; changing
the resample count after observing these outcomes would not replace this record.
All new results/facts must disclose the recorded misses and this scope limit.

Additional observed checks: dissociation-only apparent koff matched base R
optimize over four synthetic and three public T200 traces, maximum relative
difference3.733090263502627e-7. Full pytest after the final facts-state correction:
601 passed,8 existing warnings,270.07seconds. Both changed Skills passed
quick_validate. Legacy check against dev/0.11.1:136 configurations and1055
pre-existing scientific artifacts byte-identical; only kinetics facts added.

### Registration implementation clarification (recorded before the true-1:1 row completes)

Code inspection found that the generic coverage collector also computes nominal
ka2/kd2 coverage in the true-1:1/heterogeneous-fit row. There is no second-site
truth in that generating model, and the labelled first-site constants are also
not uniquely defined there. Those coverage fields are **not interpretable
coverage criteria**. They will remain in the original record, including every
numeric miss, but only the separately registered false-complex-selection metric
is applicable to that row. This clarification changes no simulations, seeds,
fits, selection thresholds or outcomes. Runtime evidence distinguishes these
invalid coverage fields from the valid false-selection criterion.

Report-only checks: screening figure visually inspected; rerender preserved all8 scientific files and did not create an affinity plot. Screening retains per-curve failures in facts and failing_items even when other curves remain rankable.

Observed true-1:1 result:987 evaluable fits,13 numerical failures,0 reportable
complex fits, false-selection rate0 (MCSE0; registered upper bound0.06387453).
No site-specific coverage fields were ultimately produced in this row because
none had accepted parameter intervals. The pre-outcome clarification remains
as an audit of the generic collector, not an omitted numeric failure.

Documentation correction: the validation index title was still labelled 0.11.0
on the 0.11.1 branch, although its new release entry and package version were
0.11.1. The current index title is corrected in 0.12.0; prior release commits and
older index entries are preserved. This was a missed version declaration in the
0.11.1 release mechanics, not a numerical change.

Per-draw calibration storage is explicit: iteration/row, evaluability,
reportability, coverage indicators, AICc improvement and diagnostics/errors are
retained for every simulation. Raw simulated traces and interval endpoints are
not serialized in this registration; its fixed seed/row/iteration reconstructs
the generating draw. Thus this is a complete per-draw outcome record, not a
stored collection of every fitted interval or raw sensorgram.

False-selection scope: the registered true-1:1 row fits the heterogeneous-ligand
model using `aicc_and_identifiability`. Its 0/987 result does not calibrate
selection among every advanced mechanism or the default diagnostic-only
comparison rule. At original registration this omitted drift, bivalent and transport alternatives;
the separately registered scope-completion supplement below supplies those checks.

### Prespecified scope-completion supplement

After reviewing the original registration's scope, a separate
`supplement_kinetics_0120_selection.py` fixes seed120261001 and1000 true1:1 draws
for each drift/bivalent/transport alternative. It uses the original noise,
time/concentration design,50 bootstrap resamples and AICc+identifiability rule.
It is committed before execution. Conditional selection among already-reportable
fits is descriptive only; the false-selection gate uses all evaluable fits.
This is an explicitly later supplement, not a replacement for the original
nine rows or their misses. It will fill the missing mechanism scope above only
to the extent actually observed; it does not validate diagnostic-only comparison.

Denominator clarification from the frozen fast-transport row:1000 fits completed,
705 had accepted bootstrap intervals and222 were reportable. The registration's
metric population label `evaluable` means interval-evaluable for coverage, while
its row-level `evaluable` counts completed fits. Every metric includes its actual
n, and interval availability is copied into runtime evidence. These populations
must not be conflated; the295 missing intervals and778 unreportable fits remain
visible. Coverage fails in both reported metric populations.

Supplement performance revision (no statistical retuning): full bootstrap for
already-rejected null fits was prohibitively slow. A validation-only
`fit_group(..., _selection_only=True)` now skips bootstrap only after an existing
profile/AICc/reliability rejection. This preserves selection exactly because
bootstrap can only append rejection notes; it cannot rescue such a fit. Fits
without a prior rejection still execute all50 resamples. A test compares the
full and shortcut decisions and identical parameters/profiles on rejected data,
and byte-equivalent result objects on a reportable case (passed before use).
Existing full-bootstrap draw records are preserved and resumed; seeds, data,
models, bounds, AICc rule, sample sizes and denominators are unchanged. Skipped
intervals are explicit and are not used for a coverage claim. The original
nine-row coverage calibration and default analysis path remain unchanged.

Supplement execution provenance: original registration `a3dc3f8`; tested
performance revision `747c02f`. The full-bootstrap attempt retained1104 completed
records (1000 drift,104 bivalent) in
`kinetics_selection_full_bootstrap_attempt_0.12.0.draws.jsonl`, with its original
partial summary. Resume reuses those records and completes the same registered
iterations; original coverage process was not stopped.

## Final original calibration outcomes

All nine registered rows completed: 9000 draws; 44 coverage criteria missed,
and the heterogeneous null-selection criterion passed. Stress-row misses are
retained. MCSE below is the observed binomial Monte Carlo standard error; each
metric retains its own denominator. No row or bound was retuned.

| Scenario | Parameter / population | n | Coverage/rate | MCSE | Bound | Result |
|---|---|---:|---:|---:|---:|---|
| one_to_one_drift | ka_coverage / evaluable | 1000 | 0.899000000 | 0.009528851 | 0.936215951 | MISS |
| one_to_one_drift | kd_coverage / evaluable | 1000 | 0.877000000 | 0.010386096 | 0.936215951 | MISS |
| one_to_one_drift | ka_coverage / reportable | 879 | 0.902161547 | 0.010020803 | 0.935297799 | MISS |
| one_to_one_drift | kd_coverage / reportable | 879 | 0.882821388 | 0.010848412 | 0.935297799 | MISS |
| heterogeneous_ligand | ka_coverage / evaluable | 1000 | 0.902000000 | 0.009401915 | 0.936215951 | MISS |
| heterogeneous_ligand | ka2_coverage / evaluable | 1000 | 0.906000000 | 0.009228434 | 0.936215951 | MISS |
| heterogeneous_ligand | kd_coverage / evaluable | 1000 | 0.903000000 | 0.009359006 | 0.936215951 | MISS |
| heterogeneous_ligand | kd2_coverage / evaluable | 1000 | 0.920000000 | 0.008579044 | 0.936215951 | MISS |
| bivalent_analyte | ka_coverage / evaluable | 1000 | 0.889000000 | 0.009933730 | 0.936215951 | MISS |
| bivalent_analyte | ka2_coverage / evaluable | 1000 | 0.881000000 | 0.010239092 | 0.936215951 | MISS |
| bivalent_analyte | kd_coverage / evaluable | 1000 | 0.896000000 | 0.009653186 | 0.936215951 | MISS |
| bivalent_analyte | kd2_coverage / evaluable | 1000 | 0.885000000 | 0.010088360 | 0.936215951 | MISS |
| bivalent_analyte | ka_coverage / reportable | 883 | 0.886749717 | 0.010664483 | 0.935331137 | MISS |
| bivalent_analyte | ka2_coverage / reportable | 883 | 0.877689694 | 0.011026094 | 0.935331137 | MISS |
| bivalent_analyte | kd_coverage / reportable | 883 | 0.894677237 | 0.010330325 | 0.935331137 | MISS |
| bivalent_analyte | kd2_coverage / reportable | 883 | 0.886749717 | 0.010664483 | 0.935331137 | MISS |
| mass_transport | ka_coverage / evaluable | 1000 | 0.888000000 | 0.009972763 | 0.936215951 | MISS |
| mass_transport | kd_coverage / evaluable | 1000 | 0.882000000 | 0.010201765 | 0.936215951 | MISS |
| mass_transport | km_coverage / evaluable | 1000 | 0.873000000 | 0.010529530 | 0.936215951 | MISS |
| mass_transport | ka_coverage / reportable | 890 | 0.887640449 | 0.010585918 | 0.935388938 | MISS |
| mass_transport | kd_coverage / reportable | 890 | 0.882022472 | 0.010812957 | 0.935388938 | MISS |
| mass_transport | km_coverage / reportable | 890 | 0.874157303 | 0.011117671 | 0.935388938 | MISS |
| off_rate_screening | kd_coverage / evaluable | 1000 | 0.870000000 | 0.010634848 | 0.936215951 | MISS |
| off_rate_screening | kd_coverage / reportable | 979 | 0.872318693 | 0.010666204 | 0.936068899 | MISS |
| true_one_site_fitted_heterogeneous | false_complex_selection / evaluable | 987 | 0.000000000 | 0.000000000 | 0.063874528 | PASS |
| transport_slow | ka_coverage / evaluable | 1000 | 0.879000000 | 0.010313050 | 0.936215951 | MISS |
| transport_slow | kd_coverage / evaluable | 1000 | 0.881000000 | 0.010239092 | 0.936215951 | MISS |
| transport_slow | km_coverage / evaluable | 1000 | 0.882000000 | 0.010201765 | 0.936215951 | MISS |
| transport_slow | ka_coverage / reportable | 886 | 0.876975169 | 0.011035021 | 0.935355993 | MISS |
| transport_slow | kd_coverage / reportable | 886 | 0.886004515 | 0.010676894 | 0.935355993 | MISS |
| transport_slow | km_coverage / reportable | 886 | 0.874717833 | 0.011121459 | 0.935355993 | MISS |
| transport_fast | ka_coverage / evaluable | 705 | 0.887943262 | 0.011880016 | 0.933583437 | MISS |
| transport_fast | kd_coverage / evaluable | 705 | 0.880851064 | 0.012201178 | 0.933583437 | MISS |
| transport_fast | km_coverage / evaluable | 705 | 0.879432624 | 0.012263703 | 0.933583437 | MISS |
| transport_fast | ka_coverage / reportable | 222 | 0.765765766 | 0.028424751 | 0.920744986 | MISS |
| transport_fast | kd_coverage / reportable | 222 | 0.725225225 | 0.029960454 | 0.920744986 | MISS |
| transport_fast | km_coverage / reportable | 222 | 0.693693694 | 0.030937510 | 0.920744986 | MISS |
| bivalent_short_association_stress | ka_coverage / evaluable | 955 | 0.877486911 | 0.010609866 | 0.935894934 | MISS |
| bivalent_short_association_stress | ka2_coverage / evaluable | 955 | 0.910994764 | 0.009214336 | 0.935894934 | MISS |
| bivalent_short_association_stress | kd_coverage / evaluable | 955 | 0.878534031 | 0.010570728 | 0.935894934 | MISS |
| bivalent_short_association_stress | kd2_coverage / evaluable | 955 | 0.894240838 | 0.009951409 | 0.935894934 | MISS |
| bivalent_short_association_stress | ka_coverage / reportable | 554 | 0.848375451 | 0.015237858 | 0.931480813 | MISS |
| bivalent_short_association_stress | ka2_coverage / reportable | 554 | 0.880866426 | 0.013763131 | 0.931480813 | MISS |
| bivalent_short_association_stress | kd_coverage / reportable | 554 | 0.877256318 | 0.013941450 | 0.931480813 | MISS |
| bivalent_short_association_stress | kd2_coverage / reportable | 554 | 0.889891697 | 0.013299145 | 0.931480813 | MISS |

Short-association stress: 1000 fits completed, 955 interval-evaluable and 554 reportable. Its 45 missing intervals and all eight coverage misses are retained.

Registry wording distinguishes experimental code availability from completion
of the roadmap five-gate definition. Those scientific gates remain incomplete;
the local development release and mechanical checks do not supersede them.

## Final supplemental null-selection outcomes

Registration a3dc3f8, execution-only shortcut 747c02f; 3000 registered draws,
including all 1104 original full-bootstrap records. Reportable-only selection
is descriptive because reportability includes selection. This validates only
the predeclared AICc-plus-identifiability rule on this generating design.

| Alternative | Evaluable | Selected | Rate | MCSE | Bound | Result |
|---|---:|---:|---:|---:|---:|---|
| one_to_one_drift | 1000 | 1 | 0.001000000 | 0.000999500 | 0.063784049 | PASS |
| bivalent_analyte | 1000 | 1 | 0.001000000 | 0.000999500 | 0.063784049 | PASS |
| mass_transport | 1000 | 0 | 0.000000000 | 0.000000000 | 0.063784049 | PASS |

Final full pytest: 602 passed, 8 existing warnings, 513.80 seconds. Both changed
Skills pass quick_validate; unknown arguments reject before output mutation in
the calibration, supplement, evidence updater, benchmark, legacy and build entries.

## Release checks

Final post-evidence full pytest: 602 passed, 8 warnings, 207.48 seconds.
Clean git-archive build from `c7a733fb681aeadaef2d3de00f38f19b7b475f2a`; 76 packaged Python files match archived source.
Separate external wheel-only environment smoke: 141 configurations passed.
Install-exact and clean-environment checks passed; see their committed JSON records.
Legacy: 136 configurations / 1055 scientific artifacts byte-identical against dev/0.11.1.
No unchanged-method benchmarks/calibrations rerun; no runs/, build/ or historical validation records modified.
No remote operation or PK/PD implementation.
