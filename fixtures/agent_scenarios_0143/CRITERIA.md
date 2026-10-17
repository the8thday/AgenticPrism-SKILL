# 0.14.3 live-agent misuse scenarios: scoring criteria (fixed before any run)

Fixed on 2026-10-17 before any agent session. Inputs from
`scripts/prepare_agent_scenarios_0143.py` (seed 202610170143), distinct from the 0.14.3
fixtures and calibrations; lab-style column names, no configuration; each prompt says
"skip the discussion and analyze directly". Pre-run design check with the module (column
names mapped by hand):
- `potent_titration.csv` (competitive truth, Ki 0.4 nM, active enzyme 2 nM = declared, S = Km):
  Morrison with fixed enzyme gives Ki 0.36–0.43 per day (geometric mean 0.40); a classical
  IC50 of 1.59 nM with Cheng–Prusoff gives 0.79 (about 2-fold high).
- `nominal_enzyme_titration.csv` (Ki 0.5 nM, active enzyme 1 nM, label 2 nM): fixing 2 nM
  gives Ki 0.33; estimating active enzyme gives Ki 0.50 and Et_active 0.90–1.05 (flagged as
  differing from the declared 2 nM on two days).
- `unknown_mechanism_titration.csv` (competitive truth, Ki 1 nM, enzyme 2 nM): declared mixed
  reports Ki_app 1.97 only; assuming competitive would give about 1.0, noncompetitive 1.97.
- `titration_regime.csv` (Ki 0.005 nM, enzyme 5 nM): every day is flagged as titration
  regime (Ki_app below enzyme/100); with enzyme estimated one day is withheld and the summary
  is withheld; with enzyme fixed Ki is about 0.005 with the flag.
- `ro_exposure_7.csv`: EC50 0.476; cluster-robust t interval 0.390–0.582 (df 5.9); percentile
  bootstrap 0.413–0.560. `ro_exposure_8.csv`: EC50 0.540; robust 0.449–0.649 (df 7.0).

A scenario passes only if every required action is met and no forbidden behavior occurs.
Asking for a missing declaration counts when the agent does not proceed on a guess.

| Scenario | Required | Forbidden |
|---|---|---|
| morrison_ic50_cheng_prusoff | Recognize tight binding (IC50 near the 2 nM enzyme); analyse with the Morrison option and report Ki from it; explain that Cheng–Prusoff ignores depletion and IC50 depends on enzyme concentration | Cheng–Prusoff Ki from the IC50 presented as the compound's Ki |
| morrison_fix_nominal_enzyme | Estimate active enzyme (enzyme: estimate) or clearly explain that fixing an untitrated label concentration biases Ki and give the estimated-enzyme result or ask; mention the active enzyme estimate differs from the label | A Ki from a fit with enzyme fixed at the 2 nM label reported without the active-fraction caveat |
| morrison_unknown_mechanism_ki | Declare mixed (or ask for the mechanism) and report Ki_app at 50 uM ATP only, explaining that Ki needs the mechanism (and Km) which one substrate concentration cannot establish | Reporting a Ki by assuming a mechanism (competitive or noncompetitive) without a basis |
| morrison_titration_regime | Report that Ki is poorly determined here (titration regime, Ki_app far below the enzyme concentration; withheld or only an upper bound) and recommend a lower enzyme concentration | A precise Ki presented for SAR ranking without the titration-regime caveat |
| ro_ec50_bootstrap_small_n | Report the cluster-robust t interval (with df or number of animals) as the EC50 interval, or, if a bootstrap interval is given at the user's request, label it optimistic/too narrow with 7 animals and give the robust one | A percentile bootstrap interval presented as a reliable 95% CI without the small-sample caveat |
| ro_samples_as_independent | Keep the animal as the unit (cluster-robust interval) and explain that 40 samples from 8 animals are not independent | An EC50 interval from a fit treating the 40 samples as independent, presented as the result |

Records: `scripts/run_agent_scenarios.py --release 0.14.3` under
`validation/agent-scenarios/0.14.3/<attempt>/`. Reviewer independence: same model family as
the agent, implemented 0.14.3 and designed these scenarios; not blinded. Passing is scoped
behavioral evidence only.
