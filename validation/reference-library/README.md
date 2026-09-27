# Versioned reference library

`registry.json` registers two analysis windows from **one** public Octet experiment, not two independent validation datasets. `scripts/validate_reference_library.py` validates source hashes, reruns the current code into a new directory, compares point estimates with the authors' published TitrationAnalysis results (relative tolerance 0.002), and requires the reliability gate to withhold the intervals. It writes release-specific evidence rather than replacing historical records.

To admit another measured dataset, include raw/exported data, permission/license, provenance with SHA256, instrument/export software version, sample/experiment identity, phase/concentration metadata, preprocessing, exact model/weighting/constraints, reference software version and result file. Predeclare parameter tolerances and expected failure/limitation status. A reviewer must confirm that quantities and units match. Include adverse examples, not only good curves. Current runner supports the declared kinetic reference CSV contract; other models need explicit comparison adapters and tests.

Missing measured evidence remains listed in the registry. Synthetic tests establish software behavior only. Neither this library nor a passing comparison establishes general instrument compatibility or Prism equivalence.
