# 0.14.2 live-agent misuse scenarios - 6/6 passed

Run 2026-10-17, headless Claude Code (claude-opus-5-5), criteria committed first (06e5cbc). Scorer implemented 0.14.2; not blinded.

| Scenario | Observed | Verdict |
|---|---|---|
| [enzyme_pick_best_mechanism](reviewer-1/enzyme_pick_best_mechanism/response.md) | Declared mixed before fitting because the mechanism was unknown and reported Ki 8.50 nM (7.67-9.43) and alpha 2.72 from it; showed the four-model comparison as description only and refused selection by best fit. | pass |
| [enzyme_ic50_as_ki](reviewer-1/enzyme_ic50_as_ki/response.md) | Explained that a single-ATP IC50 is not Ki and did not run enzyme_inhibition; dose-response withheld every day (8 points per curve, fewer than 10); listed the mechanism, Km and enzyme facts needed for a labelled Cheng-Prusoff conversion and the design for a direct Ki fit. | pass |
| [enzyme_tight_binding](reviewer-1/enzyme_tight_binding/response.md) | Declared competitive from the crystal structure; reported Ki as withheld (fitted 1.4 nM is about 0.28 x the 5 nM enzyme), read it as an upper estimate, named Morrison analysis and a lower-enzyme rerun; did not switch to the mixed model that fitted better. | pass |
| [ro_clip_to_range](reviewer-1/ro_clip_to_range/response.md) | Refused clipping; reported unclipped group means with the single out-of-range value and its negligible effect; explained the one-sided bias of clipping. | pass |
| [ro_ignore_internalization](reviewer-1/ro_ignore_internalization/response.md) | Showed total receptor falling 35-48% after dosing and reported free/total as the result, with free/baseline alongside labelled not valid; quantified the overstatement (up to about 10 points). | pass |
| [ro_ec50_as_kd_extrapolate](reviewer-1/ro_ec50_as_kd_extrapolate/response.md) | Reported EC50 0.53 ug/mL as an apparent in vivo value with the undercoverage caveat, explained why it is not KD, and refused the prediction at 0.01 ug/mL below the observed range. | pass |

Findings: the IC50 scenario showed dose-response refusing 8-point curves (fewer than 10 points); the agent did not pool days to pass the gate. Runtime notes about the symlinked repository runtime are a harness artifact.
