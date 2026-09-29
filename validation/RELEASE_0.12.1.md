# Release 0.12.1 — directed epitope binning

Local development release from the final dev/0.12.0 commit. PK/PD remains deferred.

## Scope

New independently usable epitope-binning Skill and analysis_type=epitope_binning.
Tandem, premix and sandwich require their matching readout, explicit background,
matched controls and prespecified thresholds. Directed matrices preserve missing,
ambiguous and asymmetric cells. Self-blocking must suppress reference response;
normalizing a self cell to one alone cannot pass the control. Equal experiment
weights and t intervals; complete-profile average-linkage clustering; conditional
panel-feature bootstrap BP; explicitly reciprocal-block connected communities.
No multiscale AU or structural epitope identity claim. No guessed vendor parser.

## Evidence gates

1. Cross-implementation: R pvclust at resampling ratio 1, using exactly the R
   resample indices in Python. Public pvclust lung data are a generic clustering
   oracle, not an epitope worked example. Record numbers after execution.
2. Published epitope worked example: UNMET. See source audit below.
3. Calibration: fixed seed121260930,1000 per row, independent n4/n8/n16 and
   missing-self-control stress. Script committed before execution. Both
   evaluable and reportable denominators retained. Record observed rows below.
4. Live-agent misuse: PENDING; no Claude CLI or substitute agents run.
5. Facts: saved results/config source hashes, directed pair states, asymmetries,
   control failures and bootstrap limitations; no refitting.

Reviewer command (run on the matching release branch; --release selects scenarios, not runtime code):
```sh
.venv/bin/python scripts/run_agent_scenarios.py /tmp/agentic-review-0121 --release 0.12.1 --attempt reviewer-1
```

## Self-review

- Valency/avidity: declared and preserved; asymmetric competition may reflect
  avidity, incomplete saturation or assay order, not a distinct structural site.
- Readout: format-specific detector meaning required; premix detects available
  antigen and tandem/sandwich detect additional binding. Input responses must
  already be an appropriate measured endpoint, not silently extracted traces.
- Shared scale: affirmative declaration, sourced span and matched reference/self
  response per observation; subtraction explicitly declared and no clipping.
- Independence: one ordered pair per declared independent experiment; duplicate
  wells refused. Panel-feature bootstrap is not experiment replication.
- Design range: missing/failed controls withhold affected pairs and clustering.
  Graph components can be non-cliques and incomplete membership is provisional.
  Threshold uncertainty and structural epitope identity are outside scope.
# Release evidence — 0.12.1 epitope binning

## Sources checked before numerical implementation

The roadmap's "Abdiche 2014 mAbs" candidate resolves in the primary literature
to Abdiche et al., PLoS ONE9:e92451 (2014), doi10.1371/journal.pone.0092451.
https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0092451
The article explicitly permits reuse under Creative Commons Attribution.
Its heat maps and node plots show classified pairs, but the checked article does
not provide a downloadable numeric response matrix with matching reference and
self-block controls. The separate Brooks/Miles/Abdiche2014 review is Drug
Discovery Today19:1040–1044, doi10.1016/j.drudis.2014.05.011, not a numeric
worked dataset. Published normalized-response worked-example gate is UNMET.
No synthetic responses will be substituted for the missing measured matrix.

R pvclust2.2.0 is installed in ignored .r-lib. Its bootstrap resamples matrix
rows/features, clusters columns, and compares original branch membership. This
is clustering stability conditional on the observed panel, not a biological
experiment confidence interval. Numerical oracle checks must name that scope.

The checked Table 3 lists antibody identities rather than measured responses and
matched controls. The article's network links blocking seen in either direction;
this implementation deliberately requires reciprocal blocking for an edge and
retains both directed measurements. Its communities therefore do not reproduce
the article's graph rule, and are not claimed as a matched published result.

## Observed oracle checks

All 3 checks passed. With 1000 matched R bootstrap samples, BP differences
were zero for 4 synthetic and 10 public-lung branches. Maximum relative
linkage-height differences were 1.966876662e-15 and 3.095360458e-15.
Directed normalization maximum absolute difference was 4.884981308e-15
(relative 1.249926126e-14). This validates conditional BP, not multiscale AU.
Targeted state/numerical tests: 10 passed. Both changed Skills validate.

Development drafts were preflighted outside the repository while 0.12.0
calibration ran; source availability was checked first. Release branching and
calibration execution remain sequential, and this calibration had not run at preregistration commit680482b.

## Observed calibration

Registration `680482b` committed before execution. Seed 121260930; 4000 draws, 1000 per row. All original outcomes and Monte Carlo SE retained; no retuning.

| Scenario | Metric / population | n | Rate | MCSE | Bound | Result |
|---|---|---:|---:|---:|---:|---|
| independent_n4 | blocking_interval_coverage / evaluable | 1000 | 0.943000000 | 0.007331507 | 0.936215951 | PASS |
| independent_n4 | blocking_interval_coverage / reportable | 1000 | 0.943000000 | 0.007331507 | 0.936215951 | PASS |
| independent_n4 | correct_communities / all | 1000 | 1.000000000 | 0.000000000 | 0.936215951 | PASS |
| independent_n8 | blocking_interval_coverage / evaluable | 1000 | 0.951000000 | 0.006826346 | 0.936215951 | PASS |
| independent_n8 | blocking_interval_coverage / reportable | 1000 | 0.951000000 | 0.006826346 | 0.936215951 | PASS |
| independent_n8 | correct_communities / all | 1000 | 0.999000000 | 0.000999500 | 0.936215951 | PASS |
| independent_n16 | blocking_interval_coverage / evaluable | 998 | 0.955911824 | 0.006498374 | 0.936202146 | PASS |
| independent_n16 | blocking_interval_coverage / reportable | 998 | 0.955911824 | 0.006498374 | 0.936202146 | PASS |
| independent_n16 | correct_communities / all | 1000 | 0.996000000 | 0.001995996 | 0.936215951 | PASS |
| missing_self_control | blocking_interval_coverage / evaluable | 0 | unavailable | unavailable | unavailable | descriptive / unavailable |
| missing_self_control | blocking_interval_coverage / reportable | 0 | unavailable | unavailable | unavailable | descriptive / unavailable |
| missing_self_control | withheld / all | 1000 | 1.000000000 | 0.000000000 | 0.936215951 | PASS |

All10 evaluable calibration criteria passed; two missing-control coverage
entries are unavailable (n=0), not passes. In the n16 row, 998/1000 AB pairs
had reportable intervals; both excluded outcomes remain in the draw record.
Community recovery was 1.000,0.999,0.996 for n4,n8,n16. Missing-control
withholding was1000/1000. Five new script unknown-argument guards passed.

## Facts state coverage

`validate_epitope_0121_facts.py` passed all six fixture states: complete,
asymmetric, missing, failed-control, premix and tandem. Extraction ran with
the fitter patched to raise, copied primary values/disclosures/failures from
saved results, and preserved all30 scientific files byte-for-byte. This
validation-only supplement changes no packaged runtime code. Its unknown-argument
guard also passed before any output mutation.

## Release checks

Full pytest: 612 passed, 8 existing warnings, 207.02 seconds. Changed Skills pass quick_validate.
Legacy: 141 configurations / 1098 scientific artifacts byte-identical against previous final release.
Clean git-archive build `abd39718d7db302fce5f62c0fca22e87d8f90dcf`: 79 packaged Python files exact. Separate external wheel-only smoke: 147 configurations; install-exact: 42 configurations / 253 artifacts. Clean-environment check passed.
No unchanged-method benchmark/calibration rerun, no changes to runs/, build/ or historical validation records; no remote operation or PK/PD implementation.
