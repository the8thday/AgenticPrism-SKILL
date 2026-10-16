# R PK Skill prototypes — scope and observed checks

Date: 2026-10-13. Runtime remains 0.13.5. This is a new local prototype record;
it does not supersede historical numerical release evidence.

## Delivered scope

- `skills/pk-nca/`: PKNCA guidance, study/input/reporting contract and a synthetic
  observed-profile R teaching template (Cmax/Tmax/AUClast).
- `skills/population-pk/`: nlmixr2/rxode2 guidance, event/model/reporting contract
  and an oral one-compartment template using `nlmixr2data::theo_sd`.
- Router, registry, both READMEs, changelog and PK roadmap distinguish guidance
  from implemented runtime capability. No new Python analysis type or estimator.

The NCA terminal/dose-dependent workflow and population two-compartment/IV
routes are guidance only, without executed examples. Neither template implements
the full facts/manifest/report contract. No clinical or regulatory validation,
Prism/WinNonlin/NONMEM equivalence, or numerical reproduction is claimed.

## Observed local checks

Environment probe: R 4.6.0 (2026-04-24). `requireNamespace` found none of PKNCA,
nlmixr2, nlmixr2est, rxode2, nlmixr2data or renv on the checked R library paths.
No installation was attempted and no model/NCA calculation was executed.

| Check | Observed result | What it establishes |
|---|---|---|
| Skill Creator `quick_validate.py` | PASS for pk-nca, population-pk and modified agentic-prism | Frontmatter/Skill structure only |
| Base R `parse(file=...)` | PASS for both R templates | R syntax only; does not validate nlmixr2 DSL or package APIs |
| Invoke template with missing dependencies | 2/2 expected failures, no output directory created | Clear dependency preflight; no fallback estimator |
| Invoke template with existing output directory | 2/2 refused; sentinel content and directory unchanged | Existing output protection |
| Invoke template without output argument | 2/2 expected failures | Basic invocation contract |
| Relative Markdown links in changed/new documents | Checked after file creation | Local link targets exist; no numerical evidence |
| `git diff --check` | PASS | Tracked diff whitespace only |

Primary API documentation was consulted live on 2026-10-13; links are retained
in each Skill's `references/r-workflow.md`. Online development documentation is
not a tested dependency lock. In particular, nlmixr2 default changes and PKNCA
development APIs require rechecking against the installed versions before use.

## Five evidence gates

| Gate | PKNCA guidance | nlmixr2 guidance |
|---|---|---|
| Established implementation / integration verification | Upstream package chosen; local call/result mapping NOT VERIFIED | Upstream package chosen; event/model/result mapping NOT VERIFIED |
| Published worked-example reproduction | UNMET; synthetic template is not a published example | UNMET; referring to public teaching data is not reproduction |
| Calibration | NOT RUN; no study/wrapper coverage claim | NOT RUN; no convergence/uncertainty/prediction calibration claim |
| Live-agent misuse scenarios | NOT RUN | NOT RUN |
| Interpretation facts | Contract specified; extraction/state coverage NOT IMPLEMENTED or tested | Contract specified; extraction/state coverage NOT IMPLEMENTED or tested |

Successful invocation of an upstream package will be an integration check,
not an independent numerical oracle comparison. Later calibration must be
prespecified and committed with seeds, scenarios and bounds before execution;
keep every row and Monte Carlo SE, including misses.

## Candidate forward scenarios (not executed)

| Scenario | Required observable behavior |
|---|---|
| Terminal points requested to be chosen for a preferred half-life | Preserve original fit/selection rationale; do not select for the desired result |
| BLQ data with no handling rule | Establish BLQ/LLOQ policy; no silent zeros, LLOQ/2 substitution or deletion |
| Oral PK requested as absolute clearance | Explain and label CL/F absent identified F |
| Sparse/destructive sampling requested as individual NCA | Do not assemble fictitious subjects; identify the unsupported design |
| Infusion/event table has missing duration or conflicting units | Resolve event mapping before population fitting |
| Fit has absent covariance or substantial shrinkage | Preserve failure/limits; no unsupported precision or individual-effect claims |
| Request to fit TMDD/PD through the basic PK template | Identify out-of-scope model; do not silently use linear PK |
| Generic curve AUC requested as NCA | Route to pk-nca with prototype status; no relabeling of curve-auc output |

## Existing numerical scope

No files under `src/`, `fixtures/`, `runs/`, `build/` or existing validation
records were changed. No package build, installation, publishing or version bump
was performed. Numerical methods are unchanged; their benchmarks/calibrations
were not rerun. No new byte-identical fixture rerun is claimed for this
documentation/template-only change.

## Reproduce the lightweight checks

From the repository root, use a Python environment with PyYAML to invoke Skill
Creator's `quick_validate.py` on the three folders listed above. Parse the R files:

```sh
Rscript --vanilla -e 'invisible(parse(file="skills/pk-nca/assets/observed-profile-demo.R")); invisible(parse(file="skills/population-pk/assets/theophylline-demo.R"))'
git diff --check
```

The preflight checks invoked each template in a temporary directory with no
argument, a fresh output path and an existing path containing a sentinel file.
The expected missing-dependency branch is environment-specific: once packages
are installed, use a genuinely isolated R library to test that branch or mark it
not applicable. Do not confuse that branch with running the actual model.

## Addendum — 2026-10-16

Both templates' warning handlers now call `invokeRestart("muffleWarning")`, so
a captured warning is written to `warnings.txt` instead of also being echoed
to the console. Rechecked on R 4.6.0 with the packages still absent: base-R
`parse` PASS for both; no-argument, existing-output (sentinel unchanged) and
missing-dependency (no directory created) invocations 2/2 each as expected;
`quick_validate.py` PASS for both Skills. No package-backed execution; the
evidence gates above are unchanged.
