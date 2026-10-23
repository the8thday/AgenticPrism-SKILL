"""Observed 0.15.2 evidence, extracted by scripts/update_evidence_0152.py; do not edit by hand."""
EVIDENCE = {'fih_dose': {'benchmark': 'validation/fih_benchmarks_0.15.2.json',
              'calibration': 'not applicable (deterministic)',
              'published_worked_example_status': 'MET for NOAEL-HED-MRSD; UNMET for MABEL and allometry',
              'live_agent_status': '6/6 passed',
              'must_mention': ['Published worked example (NOAEL to HED only): Nair & Jacob 2016 (J Basic '
                               'Clin Pharm 7:27) rat examples reproduced at the printed precision (18 mg/kg '
                               'at 0.15 kg -> 2.5 mg/kg -> 150 mg -> 15 mg with a safety factor of 10; 50 '
                               'mg/kg / 6.2 -> 8.1 mg/kg). The FDA 2005 Table 1 divisors were transcribed '
                               'from secondary reproductions (the FDA PDF could not be downloaded) and '
                               'checked against divisor = 37/Km within the one-decimal rounding (at most '
                               '0.05). No MABEL or antibody-allometry worked example was reproduced (UNMET).',
                               'An independent base-R/deSolve recomputation agreed on 29/29 reported '
                               'numbers. Calibration is not applicable to these deterministic calculations; '
                               'declared-input sensitivity tables are reported instead.',
                               'Live-agent misuse scenarios: 6/6 passed (not blinded); see '
                               'validation/agent-scenarios/0.15.2/RESULTS.md.']}}
