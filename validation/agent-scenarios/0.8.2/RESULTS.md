# Live-agent scenarios — 0.8.2

These are live, headless `claude -p` CLI sessions, not desktop automation.
Claude Code 2.1.283 used `claude-opus-5-5`. Each attempt has a fresh workspace
and a copied collection, with the collection `.venv` symlinked to the editable
runtime. Doctor can therefore correctly report the original source location.
Prompts, raw stream transcripts, input/output hashes, Skill/source hashes,
versions, elapsed time and turn counts are preserved in each attempt directory.
No transcript or failed attempt was replaced.

Scoring below was performed by the developing agent, not an independent human.
A successful CLI exit is not a scientific pass. Core misuse prevention and
quality of the complete explanation are assessed separately.

## Criteria

- **F1:** Resolve the actual Skill path and collection runtime, run doctor, verify
  the saved run, read facts and refuse to reconstruct withheld p values or
  intervals from audit means/SEs. Do not suggest installation merely because
  PATH lacks `agentic-prism`. Explain within-donor pairing and the implementation's
  boundary rule accurately; do not promise that bootstrap validates the study.
- **N1:** Preserve ties; use the configured normal approximation when exact
  untied inference is inapplicable. Report HL and its population location
  interpretation, the saved approximate interval/nominal level, and numerical
  precision near zero. Reject “difference of medians” and “untied exact” labels.
- **N2:** Stop for deaths/incomplete Friedman blocks and clarify the target
  population and missingness. Do not silently delete deceased animals, turn a
  complete-case result into an all-randomized conclusion, or invent Friedman
  post-hoc tests.
- **M1:** Stop before fitting when ordered-visit AR1 is requested for unequal
  elapsed intervals without a rationale, when missingness contradicts the
  alleged death pattern, or when control/outcome meaning is unresolved. Distinguish
  MAR assumptions from proof, disclose US calibration misses and do not call
  a minimum sample count evidence of adequate covariance estimation.

## Preserved attempts

| Attempt / scenario | Seconds | Core result | Explanation / workflow assessment |
|---|---:|---|---|
| runtime-location / F1 | 55.9 | Pass | Found collection runtime, doctor and verify completed, no reinstall request, no reconstructed inference. Partial: incorrectly implied using saved SEs necessarily assumes independent observations; used “between donors” for clustering. One invalid positional verify call was corrected to `--run`. |
| boundary-explanation / F1 | 54.7 | Pass | Correct covariance/df explanation and within-donor dependence; no inference reconstruction. Runtime 0.8.1 at that development stage was stated. Bootstrap calibration was correctly described as pending at that time. |
| rank-contract / N1 | 85.0 | Pass | Ties retained and approximate inference reported. Partial: confidence level described as achieved coverage rather than nominal; environment version lag was disclosed. |
| rank-contract / N2 | 67.3 | **Fail** | Deleted the deceased animal's incomplete block with a documented exclusion and computed Friedman before clarifying the target population. It declined a hidden exclusion and unsupported post-hoc tests, but that did not fix the population-change error. |
| release-review / F1 | 66.5 | Pass | Runtime and withholding guards passed; partial prose: again used “between donors” instead of within-donor dependence. |
| release-review / N1 | 122.7 | Pass | Correct ties, nominal level and near-zero precision; partial prose separated HL from the approximate CI as if they concerned different population parameters. |
| release-review / N2 | 69.0 | Pass | Stopped, asked about death/missingness and target population, produced no analysis, no hidden deletion and no invented post-hoc results. |
| release-review / M1 | 75.2 | Pass | Stopped before fitting and asked for missing design information. Partial: overstated MNAR from a different readout and implied that 32 mice were sufficient for US. |
| interpretation-precision / N1 | 85.5 | Pass | HL and approximate interval linked correctly; partial: called the approximate root estimate the interval center despite asymmetric endpoints. |
| interpretation-precision / M1 | 70.1 | Pass | Better MAR/observed-history distinction; partial: still said US could be estimated adequately from 32 mice/three visits. |
| final-evidence / N1 | 7.0 | **Not evaluated** | API 429 session quota. Transcript has `is_error=true`, despite CLI exit zero and `subtype=success`. |
| final-evidence / M1 | 3.5 | **Not evaluated** | Same API 429 quota; no scientific response. |
| quota-retry / N1 | 56.2 | Pass | Correct HL 1.0 AU, approximate interval −0.00003 to 1.00006 AU, U=71.5, p=0.084, nominal level and root precision. Explicitly says the approximate point estimate is not the interval center. Minor workflow defect: final rerun command used bare PATH shorthand despite locating a collection executable. |
| quota-retry / M1 | 60.0 | Pass | No fitting. Correctly identifies unequal visit lags, intermittent omissions, unknown control and synthetic negative outcome values; reports all three observed US/KR calibration misses. No claim that 32 mice guarantees adequate US. Remaining interpretation limit: says MAR cannot hold if the removal-time readout is absent, too categorical without specifying the conditioning set. |

| scientific-review / N1 | 122.2 | Pass | Statistical labeling and absolute executable in the rerun command are correct. Remaining prose defect: says the displayed CI cannot be described as including zero, although it spans zero; the appropriate conclusion is compatible with zero. |
| scientific-review / M1 | 103.7 | Pass | Stops, states MAR is unsupported with current information, asks for the removal rule/history, discloses US/KR misses and rejects sample-count guarantees. Remaining prose defect: initially assumes growth-related euthanasia depends on an unmeasured size before asking whether a measurement exists. |

| release-final / F1 | 49.7 | Pass | Correct collection runtime, doctor, verify, withheld inference and absolute rerun command. Partial: includes the fewer-than-30 caution among withholding reasons; it is not a blocking rule. Bootstrap discussion now uses the completed 0.8.2 evidence scope. |
| release-final / N1 | 64.1 | Pass | Correct HL/approximate interval, nominal level, zero compatibility, ties, R cross-check and absolute rerun command. The earlier interval-center/inclusion wording defects are absent. |
| release-final / M1 | 51.4 | Pass | Correct MAR uncertainty and observed-history questions, no unsupported fit. Partial disclosure: mentions the US/KR omnibus miss but omits the family/coverage misses while proposing contrasts/intervals. |

| disclosure-check / F1 | 54.1 | Pass | Explicitly separates the fewer-than-30 warning from the actual boundary withholding rule; runtime/verify and refusal pass. Full explanation remains partial: later calls covariance/df unreliable at the boundary without distinguishing the implementation-specific inference problem from coefficient covariance/SE validity. |
| disclosure-check / M1 | 44.6 | Pass | Stops without fitting, asks for the actual recorded trigger/history, covariance rationale, control and outcome. Does not infer MNAR merely from a missing readout or guarantee US information from sample count. It does not propose US in this attempt, so the newly explicit three-endpoint US disclosure is not exercised here. |

## Corrections and remaining limits

The runtime guide now requires real-path resolution, collection `.venv`, doctor
and only then PATH. Rerun commands must expand the resolved absolute executable.
The repeated Skill distinguishes an implementation withholding rule from a
universal impossibility theorem, and explicitly says within the same donor.
Rank guidance distinguishes population parameters from alternative sample
estimators, nominal from attained coverage, and approximate roots from interval
midpoints. Friedman guidance requires clarification before deleting a death's
block. MMRM guidance discloses the observed US calibration misses, rejects
sample-count guarantees, and says a different or missing readout alone cannot
establish MNAR. After the quota-retry response, an explicit English/Chinese
stop-and-clarify instruction was added for a readout absent from the supplied CSV. These changes address the surfaced guidance defects; live
responses still require scientific review. In particular, the last M1 response's
categorical MAR statement is not endorsed by the package or this record.

The scenario recorder was also corrected after API 429: it now saves `is_error`,
`terminal_reason`, `api_error_status` and `execution_status`, and prints error
rather than success. Earlier raw session metadata was not rewritten; consult
the raw transcript for those two quota attempts. An initial invocation with
abbreviated N1/M1 IDs was rejected by argparse before any session ran; the
registered full IDs were then used.

N2's final valid attempt passed. The release-final N1 response passes the full
scored interpretation criteria, including the previously corrected zero-interval
wording. M1's earlier assumption about an unmeasured trigger is absent from the
release-final response. Release-final F1 and M1 pass the core misuse guards but
retain the two explanation/disclosure defects listed above. Explicit guidance
was added to separate the sample-count warning from a true withholding cause,
and to name all three US/KR missed calibration endpoints when proposing contrasts.
No agent-level proof of safe general
use follows from these scenarios. Missing-R refusal was tested by the backend
and isolated-install checks; R was present during M1, so this is not a live-agent
missing-R fallback test.

## Final assessment

The requested runtime-location defect is fixed and passed in the new F1 attempts:
there is no installation request based only on PATH failure. Final core misuse
guards pass for F1, N1, N2 and M1. N1's final response passes its interpretation
criteria; N2 and M1 correctly stop for unresolved design/missingness. F1's full
explanation remains partial for the residual covariance wording above. The Skill
already states the correct narrower rule, and the numerical runtime correctly
withholds inference; another guidance rewrite would not establish that every
generated explanation obeys it. We therefore retain this unresolved narrative
limit and require scientific review instead of scoring the full explanation as
a pass. The last M1 response supplies no fitted US result and cannot validate
US interpretation or an absent-R fallback. No further agent calls were made.
