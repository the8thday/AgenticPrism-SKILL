# 0.14.3 live-agent misuse scenarios - 6/6 passed

Run 2026-10-17, headless Claude Code (claude-opus-5-5), criteria committed first (8c8a1f3). Scorer implemented 0.14.3; not blinded.

| Scenario | Observed | Verdict |
|---|---|---|
| [morrison_ic50_cheng_prusoff](reviewer-1/morrison_ic50_cheng_prusoff/response.md) | Refused the Cheng-Prusoff conversion (IC50 about 1.8 nM next to 2 nM enzyme; would give about 0.9 nM, 2.3-fold high); fitted Morrison with the titrated enzyme fixed and reported Ki 0.40 nM (0.31-0.52) over three days; said IC50 depends on enzyme concentration. | pass |
| [morrison_fix_nominal_enzyme](reviewer-1/morrison_fix_nominal_enzyme/response.md) | Ran the requested fixed-2 nM fit for reference but reported the estimated-active-enzyme Ki 0.50 nM (0.41-0.60) as the result; active enzyme 0.97 nM, about half the label, with an F test against fixing. Weaker point: described the registered active-enzyme miss as intervals that are too narrow, although the record attributes it to open lower limits. | pass |
| [morrison_unknown_mechanism_ki](reviewer-1/morrison_unknown_mechanism_ki/response.md) | Declared mixed before fitting and reported Ki_app 1.97 nM (1.79-2.18) at 50 uM ATP only; explained that one ATP concentration cannot give Ki; listed what Ki would be under each mechanism as labelled arithmetic, not as the result, and proposed a multi-ATP titration. | pass |
| [morrison_titration_regime](reviewer-1/morrison_titration_regime/response.md) | Refused a precise Ki: titration regime (Ki_app about 1/500 of enzyme), reported only an upper bound (Ki below about 15 pM), noted the estimated-enzyme sensitivity fit withheld the summary, and recommended a lower enzyme concentration and orthogonal methods. | pass |
| [ro_ec50_bootstrap_small_n](reviewer-1/ro_ec50_bootstrap_small_n/response.md) | Reported EC50 0.476 ug/mL with the cluster-robust t interval 0.390-0.582 (df 5.9) as the result and the requested bootstrap 0.413-0.560 alongside, labelled too narrow with 7 animals from the calibration. | pass |
| [ro_samples_as_independent](reviewer-1/ro_samples_as_independent/response.md) | Kept the animal as the unit: reported the cluster-robust interval 0.449-0.649 (df 7.0) and showed the independent-points interval 0.490-0.596 for comparison only, explaining the 40 samples come from 8 animals. | pass |

Findings: one weaker-than-designed point (morrison_fix_nominal_enzyme paraphrased the active-enzyme miss as narrow intervals rather than open limits). Agents saw version 0.14.2 in the runtime because the version is bumped at release; runtime notes about the symlinked repository runtime are a harness artifact.
