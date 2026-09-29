# Release 0.13.0 — combination reference models and HTS plate QC

Local development release from final dev/0.12.1. No PK/PD/NCA/TMDD/NLME added.

## Scope

New drug-combination and import-plate Skills. Combination uses the existing
separable 4PL fitter, named Bliss/Loewe/HSA/ZIP references and equal-weight
independent-matrix Student t uncertainty. Unsupported fits withhold dependent
scores; no independent uncertainty from public repeated rows lacking experiment
identifiers. HTS has explicit roles/layout, Z-prime, plug-in SSMD, median-polish
B scores, row/column/edge diagnostics and exploratory predictive-t/BH hits.
Dose-response and ELISA saved-result facts complete the oldest-six retrofit.

## Evidence gates

1. Cross-implementation: synergyfinder3.20.0 on fixture and packaged ONEIL public
   data; base R medpolish/pt/p.adjust for HTS. cellHTS2 unavailable for installed
   Bioconductor3.23; base R medpolish used as expressly allowed. Record results.
2. Published: ONEIL packaged public example checked with its licence. Its four
   repeated rows per dose pair lack independent-experiment IDs, so no biological
   replication or confidence interval is inferred. Separate matched published HTS
   worked-result gate UNMET unless source audit establishes one.
3. Calibration: seed130260930,1000 per row, fixed/free single-agent plateaus,
   independent matrices, null HTS hits and failed-control stress. Commit before
   execution; all draws and evaluable/reportable coverage retained.
4. Live-agent misuse: PENDING; no Claude CLI or substitute agents run.
5. Facts: saved artifacts only, audit quantities excluded from primary claims;
   new-model and old dose/ELISA withholding states tested.

Reviewer command (run on the matching release branch; --release selects scenarios, not runtime code):
```sh
.venv/bin/python scripts/run_agent_scenarios.py /tmp/agentic-review-0130 --release 0.13.0 --attempt reviewer-1
```

## Self-review

- Valency/avidity: no binding KD is inferred from inhibition; combination scores
  do not identify a molecular interaction mechanism or clinical benefit.
- Readout: a required description states the detector/biological endpoint or explicitly records that it is unknown; saved facts retain it. Explicitly normalized percent inhibition on one scale, no implicit
  viability conversion. HTS direction and units require declarations.
- Shared scale: all matrices share conditions/dose grids; Loewe needs compatible
  invertible effect ranges and both declared and fitted plateau checks.
- Independence: matrix IDs alone do not prove independence; supported declaration
  required for intervals. Wells are screening units, not independent biological
  confirmation. Repeated unit IDs are refused in HTS.
- Design range: complete dose grids include single-agent axes and positive doses;
  4PL midpoint, boundaries, flatness and rank gates withhold dependent models.
  Unreplicated or unsupported matrices cannot yield a synergy claim.
- Plate effects: median polish requires randomized layout and mostly inactive
  wells; cannot remove a biological row/column effect honestly. Failed QC and
  correction failures withhold hits. Predictive t/BH is exploratory and its
  dependence/normality assumptions and observed calibration must be disclosed.
# Release evidence — 0.13.0 pharmacology and HTS QC

## Source and licence audit before numerics

Bioconductor synergyfinder3.20.0 is installed in ignored .r-lib. Its package
licence is Mozilla Public License2.0. Primary package and vignette:
https://bioconductor.org/packages/release/bioc/html/synergyfinder.html
The packaged ONEIL_screening_data has200 rows, two drug-pair blocks and repeated
concentration pairs. Its help cites O'Neil et al.2016, Molecular Cancer
Therapeutics15:1155–1162, PMID26983881. Reuse is scoped to the package-distributed
example under its stated licence; original full-screen rights are not inferred.
The example lacks explicit independent-experiment identities. It can validate
point calculations, but repeated rows alone do not establish biological
replication or support replicate-level uncertainty claims.

The package ZIP implementation fits single-agent4PLs, forms the Bliss reference
from those fitted responses, fits conditional dose-response curves along both
matrix directions with fixed condition baseline, and averages the two fitted
responses before subtracting the reference. ZIP is not implemented as simply
observed response minus fitted Bliss. Edge single-agent scores are zero.

Base R stats::medpolish is the primary B-score residual oracle. Z-prime and SSMD
will be matched to explicit formulae and sources. A public measured plate with
usable licence and matched printed QC results has not yet been verified; until
then its published worked-example gate is UNMET. No synthetic plate substitutes
for this evidence.

Optional cellHTS2 install was attempted with BiocManager1.30.27/R4.6.0/
Bioconductor3.23. Repository response: package cellHTS2 is not available for
Bioconductor3.23. Use base R medpolish as explicitly requested; no substitute
cellHTS2 equivalence claim. Z-prime source: Zhang/Chung/Oldenburg1999,
doi10.1177/108705719900400206. SSMD source: Zhang2007,
PMID17276655 (https://pubmed.ncbi.nlm.nih.gov/17276655/). These establish equations,
not a matched publicly licensed raw-plate worked example.
Additional public numerical example: archived cellHTS2 2.68.0 RELEASE_3_19
commit e73938774eea1c409e1205ad781d9aa97551430d includes KcViab
FT01-G01.txt (384 measured ATP-luminescence wells), described in Boutros2004
Science303:832–835, PMID14764878 and the package vignette
https://bioconductor.posit.co/packages/3.19/bioc/vignettes/cellHTS2/inst/doc/cellhts2.pdf .
Package licence Artistic-2.0; source files/hashes retained. Dataset Description
has an empty separate License field: reuse is limited to the package-distributed
example, without inferring rights to the original full screen. Plateconf has
only B01 negative and B02 positive, so a complete HTS workflow is refused
(insufficient replicated controls). Base R medpolish/B-score numerical check
is usable; this does not establish a matched published QC/hit worked result.

Pre-release scratch oracle check retained: ONEIL block1 agrees for Bliss/HSA,
but its Loewe inverse is incomplete and its ZIP local fits differ beyond the
registered 0.05 percentage-point tolerance. Both are withheld by applicability
gates. The final public benchmark includes both original blocks separately;
block1 failures are retained, and no tolerance was changed to make them pass.
A copied synthetic Loewe declaration in the unselected public config was removed
during preflight: public maximal-effect compatibility is not established.

The second public ONEIL block also fails Loewe completeness (17/25 finite
cells) and ZIP agreement (maximum absolute difference 76.3297188901 percentage
points; relative 0.7824330464). The constrained existing 4PL workflow withholds
these unsupported fits. No public ZIP/Loewe equivalence is claimed, no tolerance
was widened, and the discrepancy remains an unmet numerical scope. Formal
on-branch results will be recorded separately.

Design-range preflight correction: HTS now requires declared 96/384-well
dimensions and refuses an entirely absent row/column. An observed smaller
rectangle cannot redefine plate edges or silently pass completeness. This was
fixed before calibration registration/execution; simulation values are unchanged.

ZIP retains every conditional-slice 4PL parameter, fixed-other-agent dose,
SSE, support flag and failure reason in saved results. Unsupported slice fits
remain inspectable without promoting their scores to reportable summaries.

## Observed numerical comparisons

16/20 checks pass; four public Loewe/ZIP checks fail. Every original
discrepancy is retained, including incomplete inverse ranges. Public point
agreement is established for Bliss/HSA only. Tolerances were not widened.

| Dataset | Model/metric | Max absolute difference | Max relative difference | Result |
|---|---|---:|---:|---|
| combination | Bliss | 2.88657986403e-14 | 7.62020062026e-14 | PASS |
| combination | Loewe | 0.466202909285 | 1.29775488431 | PASS |
| combination | HSA | 4.97379915032e-14 | 3.25752460299e-14 | PASS |
| combination | ZIP | 0.00078142662943 | 0.00839985632418 | PASS |
| public_oneil | Bliss | 4.97379915032e-14 | 2.90037858425e-15 | PASS |
| public_oneil | Loewe | 0.663743847824 | 0.0253336364093 | FAIL |
| public_oneil | HSA | 7.1054273576e-15 | 2.17939507798e-16 | PASS |
| public_oneil | ZIP | 0.299139007025 | 0.0238235190605 | FAIL |
| public_oneil2 | Bliss | 4.26325641456e-14 | 1.93327675572e-15 | PASS |
| public_oneil2 | Loewe | 0.45217231864 | 0.007650319798 | FAIL |
| public_oneil2 | HSA | 1.42108547152e-14 | 1.56005818646e-16 | PASS |
| public_oneil2 | ZIP | 76.3297188901 | 0.782433046361 | FAIL |
| synthetic_hts | residual | 5.3290705182e-15 | 4.77145653359e-15 | PASS |
| synthetic_hts | b_score | 4.61852778244e-14 | 4.31104100754e-15 | PASS |
| synthetic_hts | p_value | 5.55111512313e-16 | 2.02238848741e-15 | PASS |
| synthetic_hts | adjusted_p | 4.4408920985e-16 | 5.44780073568e-15 | PASS |
| synthetic_hts | zprime | 4.4408920985e-16 | 5.35316217114e-16 | PASS |
| synthetic_hts | ssmd | 3.19744231092e-14 | 1.2859707052e-15 | PASS |
| public_KcViab_FT01_G01 | residual | 0 | 0 | PASS |
| public_KcViab_FT01_G01 | b_score | 4.79616346638e-14 | 1.37093371871e-14 | PASS |

Loewe absolute tolerance1.1 percentage points reflects synergyfinder’s
100-response-level grid; synthetic max absolute difference0.466203 is not
exact numerical equivalence. ZIP tolerance0.05 points passes the synthetic
fixture but fails both public blocks. The workflow withholds those unsupported
fits; withholding is not treated as agreement with the oracle.

34 targeted tests passed, covering seven new-model/result states, saved-only
facts in each, analytic reference limits and 12 dose/ELISA withholding states.
All five changed Skills pass quick_validate. Development preflight occurred
in an external scratch copy while earlier checks ran; release branches and
calibration execution remain sequential. Calibration had not run at preregistration commit4134ce5.

Additional status-only facts checks: six passed for dose and ELISA summaries
without a reportable boolean. Saved withholding/independence states and absent
intervals remain explicit; a valid one-experiment point stays descriptive.
No statistical implementation or registered calibration scenario changed.
HTS raw-plate and row/column plots were visually inspected; rerender with the
fitter disabled preserved all5 scientific files.

## Observed calibration

Registration `4134ce5` committed before execution. Seed 130260930; 4000 draws, 1000 per row. All original outcomes and Monte Carlo SE retained; no retuning.

| Scenario | Metric / population | n | Rate | MCSE | Bound | Result |
|---|---|---:|---:|---:|---:|---|
| combination_fixed_plateaus | Bliss_coverage / evaluable | 1000 | 0.949000000 | 0.006956939 | 0.936215951 | PASS |
| combination_fixed_plateaus | Bliss_coverage / reportable | 1000 | 0.949000000 | 0.006956939 | 0.936215951 | PASS |
| combination_fixed_plateaus | Loewe_coverage / evaluable | 1000 | 0.942000000 | 0.007391617 | 0.936215951 | PASS |
| combination_fixed_plateaus | Loewe_coverage / reportable | 1000 | 0.942000000 | 0.007391617 | 0.936215951 | PASS |
| combination_fixed_plateaus | HSA_coverage / evaluable | 1000 | 0.915000000 | 0.008819014 | 0.936215951 | MISS |
| combination_fixed_plateaus | HSA_coverage / reportable | 1000 | 0.915000000 | 0.008819014 | 0.936215951 | MISS |
| combination_fixed_plateaus | ZIP_coverage / evaluable | 1000 | 0.951000000 | 0.006826346 | 0.936215951 | PASS |
| combination_fixed_plateaus | ZIP_coverage / reportable | 1000 | 0.951000000 | 0.006826346 | 0.936215951 | PASS |
| combination_free_plateaus | Bliss_coverage / evaluable | 1000 | 0.956000000 | 0.006485677 | 0.936215951 | PASS |
| combination_free_plateaus | Bliss_coverage / reportable | 1000 | 0.956000000 | 0.006485677 | 0.936215951 | PASS |
| combination_free_plateaus | Loewe_coverage / evaluable | 1000 | 0.953000000 | 0.006692608 | 0.936215951 | PASS |
| combination_free_plateaus | Loewe_coverage / reportable | 1000 | 0.953000000 | 0.006692608 | 0.936215951 | PASS |
| combination_free_plateaus | HSA_coverage / evaluable | 1000 | 0.928000000 | 0.008174105 | 0.936215951 | MISS |
| combination_free_plateaus | HSA_coverage / reportable | 1000 | 0.928000000 | 0.008174105 | 0.936215951 | MISS |
| combination_free_plateaus | ZIP_coverage / evaluable | 1000 | 0.953000000 | 0.006692608 | 0.936215951 | PASS |
| combination_free_plateaus | ZIP_coverage / reportable | 1000 | 0.953000000 | 0.006692608 | 0.936215951 | PASS |
| hts_null_hits | false_hit / evaluable | 1000 | 0.100000000 | 0.009486833 | 0.063784049 | MISS |
| hts_null_hits | false_hit / reportable | 1000 | 0.100000000 | 0.009486833 | 0.063784049 | MISS |
| hts_failed_controls | withheld / evaluable | 1000 | 1.000000000 | 0.000000000 | 0.936215951 | PASS |
| hts_failed_controls | withheld / reportable | 0 | unavailable | unavailable | unavailable | descriptive / unavailable |

## Calibration interpretation and unmet statistical scope

19 pass/fail criteria were evaluable:13 passed and6 missed; the conditional
withholding row has n=0 and is descriptive/unavailable. The six misses are
HSA coverage0.915 (MCSE0.008819014) with fixed plateaus, HSA0.928
(MCSE0.008174105) with free plateaus, and HTS global-null false-hit probability
0.100 (MCSE0.009486833), each repeated for all-evaluable/reportable populations
(n=1000 in each). Coverage bound0.936215951; false-hit bound0.063784049.

Under the global null, FDR equals the chance of at least one false discovery.
Thus the requested validated FDR-control scope is UNMET for these nominal
predictive-t/BH screens. Existing indicators remain explicitly experimental,
with the miss copied into results/facts; they are not validated hits. No cutoff,
fit, input design or simulation row was changed after observing this failure.
Likewise, HSA intervals have not demonstrated nominal95% coverage. The other
registered combination intervals passed only their two declared designs.
All1000 failed-control stress plates were withheld. No unreported rows were dropped.

## Release checks

Full pytest: 652 passed, 8 existing warnings, 208.87 seconds. Changed Skills pass quick_validate.
Legacy: 147 configurations / 1128 scientific artifacts byte-identical against previous final release.
Clean git-archive build `20f0a4cb30231d104f924a8f0d9b81aa622335c7`: 83 packaged Python files exact. Separate external wheel-only smoke: 154 configurations; install-exact: 55 configurations / 370 artifacts. Clean-environment check passed.
No unchanged-method benchmark/calibration rerun, no changes to runs/, build/ or historical validation records; no remote operation or PK/PD implementation.

[Four-release consolidated review](RELEASE_SERIES_0.11.1-0.13.0.md) lists
all76 calibration misses, source/gate limits, self-review and reviewer commands.
`series_artifact_sources_0.13.0.json` verifies all321 tracked runtime files
(including R resources) across the four wheels against their checked branches.

## Reviewer live-agent run (2026-09-29)

The PENDING live-agent gate above was run by the reviewer with Claude Code: **4/4 passed**. See [agent-scenarios/0.13.0/RESULTS.md](agent-scenarios/0.13.0/RESULTS.md). The implementer did not run or score these sessions.
