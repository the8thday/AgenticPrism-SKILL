# 0.14.2 — enzyme inhibition and receptor occupancy

Baseline `67d0a13` (0.14.1, local). Plan: [plans/RELEASE_0.14.2_PLAN.md](../plans/RELEASE_0.14.2_PLAN.md).
New analysis types only; existing default outputs are checked byte for byte at the end.

## Item 1 — enzyme inhibition: `enzyme_inhibition`

[Skill](../skills/enzyme-inhibition/SKILL.md), [contract](../skills/enzyme-inhibition/references/input-and-model.md).
Declared competitive, noncompetitive, uncompetitive or mixed (alpha) model, global
unweighted least squares across inhibitor concentrations per independent experiment,
profile-F intervals, descriptive model comparison against mixed, tight-binding and range
gates, geometric summaries across experiments. Prespecified in `27544e1`.

| Gate | Observed scope |
|---|---|
| Published worked example | GraphPad Prism 4 Step-by-Step Example "Computing Ki for a Competitive Enzyme Inhibitor" (adapted from Miller 2003; document SHA-256-pinned and downloaded at validation time, rates transcribed in the script, not redistributed): **partial at printed precision** — Km 236.7 matches; Ki 2.743 versus printed 2.744 and Vmax 79.93 versus 79.92 (one unit in the last printed digit). Residual SS 118.23297 here versus 118.23331 at the printed values: the printed fit is slightly short of the least-squares minimum. The plan's pre-freeze sentence overstated the match and carries a correction. [Record](enzyme_benchmarks_0.14.2.json) |
| Established implementation | Base R `nls` on log parameters with profile `confint`: 14/14 experiment fits (GraphPad data and 13 fixture experiments, all four mechanisms) within 1.3e-5 relative for estimates and interval limits. The first oracle run failed when R started exactly at the Python optimum (nls error); the oracle now starts 0.05 log units away and falls back to `port`. |
| Calibration | Seed 141401, 1000 draws per scenario, all kept: **9/10 bounded criteria passed.** Per-experiment Ki profile coverage 0.960 (competitive), 0.947, 0.947, 0.939 (mixed); Km 0.961; alpha 0.948; false rejection of a true declared mechanism against mixed 0.021, 0.057, 0.014. **Miss retained:** summary geometric Ki coverage 0.874 (bound 0.936), counting withheld summaries as not covered. Descriptive: declared competitive with mixed truth withheld Ki in every experiment; tight binding with enzyme = Ki declared withheld Ki in every experiment, undeclared it was reported in 0.885 with mean log10 error +0.051; 8 % multiplicative noise gave Ki coverage 0.920; inhibitor only to 0.5 Ki gave 0.952 with the range flag always raised. [Record](enzyme_calibration_0.14.2.json) |
| Calibration, post hoc | From the same saved draws (no new simulation): 0.924 of three-experiment summaries were estimated and 0.946 of those covered; the miss comes from withholding (each experiment ~2 % false rejection), not interval width. [Record](enzyme_summary_posthoc_0.14.2.json) |
| Live-agent misuse | 3/3 enzyme scenarios passed (best-fit mechanism selection refused; IC50 not treated as Ki; tight binding withheld). [Record](agent-scenarios/0.14.2/RESULTS.md) |
| Interpretation facts | Extraction-only; facts equal saved results |

## Item 2 — receptor occupancy: `receptor_occupancy`

[Skill](../skills/receptor-occupancy/SKILL.md), [contract](../skills/receptor-occupancy/references/input-and-model.md).
Free/baseline, free/total and bound/saturation formats from flow-cytometry MFI with declared
background, no clipping, subject-level summaries, optional exposure-occupancy EC50 with a
subject-cluster bootstrap. Prespecified in `da73dfc`; before the calibration ran, `5a2f999`
made zero between-subject variance withhold a summary and added background-tube noise to the
fixtures (the calibration script itself is unchanged since `da73dfc`).

| Gate | Observed scope |
|---|---|
| Established implementation | Base-R oracle (formulas, `t.test` summaries, one-parameter EC50 by `optimize`): 5/5 fixtures within 1.0e-8; bootstrap reproducible from its seed. The first oracle attempt used one-parameter `nls`, which failed (step factor); the first comparison then failed on a degenerate fixture (predose bound signal exactly equal to background), which led to the `5a2f999` change. |
| Published worked example | UNMET: no printed RO calculation with its data was located. |
| Calibration | Seed 141501, 1000 draws per scenario, all kept: **3/5 bounded criteria passed.** 168 h group-mean RO coverage 0.954 (free/baseline), 0.949 (bound/saturation), 0.952 (free/total with 40 % receptor loss). **Misses retained:** EC50 cluster-bootstrap coverage 0.770 and 0.767 with 8 subjects (withheld fits counted as not covered). Descriptive: free/baseline under 40 % receptor loss overstated RO by 6.6 points (coverage 0.205), which free/total corrected (−0.3); clipping to 100 % near saturation biased the mean by −0.12 points versus +0.00; with 6 subjects EC50 coverage 0.673. [Record](receptor_occupancy_calibration_0.14.2.json) |
| Calibration, post hoc | From the same draws: EC50 reportable in 0.877 / 0.864 / 0.782 of datasets and coverage among reportable fits 0.878 / 0.888 / 0.861 — the percentile cluster bootstrap is too narrow with 6–8 subjects; stated in every run. Follow-up: a better small-sample EC50 interval. [Record](receptor_occupancy_ec50_posthoc_0.14.2.json) |
| Live-agent misuse | 3/3 RO scenarios passed (no clipping; free/total chosen under receptor loss; EC50 not called KD and no extrapolation). [Record](agent-scenarios/0.14.2/RESULTS.md) |
| Interpretation facts | Extraction-only; facts equal saved results |

## Checks across 0.14.2

- Full suite: `.venv/bin/python -m pytest -q` 963 passed (2026-10-17).
- Unchanged methods: 247 configurations, 1,750 byte-identical scientific artifacts against
  0.14.1 (`67d0a13`), including every 0.14.0 and 0.14.1 fixture ([record](legacy_fixture_checks_0.14.2.json)).
- Live-agent misuse scenarios 6/6 (criteria `06e5cbc`, not blinded).
- Follow-ups: a better small-sample interval for the RO exposure-occupancy EC50; Morrison
  tight-binding analysis for enzyme inhibition.
