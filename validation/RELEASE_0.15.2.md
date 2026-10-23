# 0.15.2 — First-in-human dose calculations: `fih_dose`

Baseline `3215010` (0.15.1). Plan: [RELEASE_0.15.2_PLAN.md](../plans/RELEASE_0.15.2_PLAN.md).
Evidence module `src/agentic_prism/evidence_0152.py` (`scripts/update_evidence_0152.py`).

## Scope

New design-only analysis type `fih_dose` (Skill `first-in-human-dose`, Python): NOAEL to HED by
the FDA 2005 Table 1 body-surface-area divisors, the body-weight exponent form or mg/kg (refused
otherwise for proteins > 100 kDa), declared species selection (a non-most-sensitive choice is a
failing item) and safety factor (below 10 is a failing item), MRSD in mg/kg and mg; MABEL from
receptor occupancy (1:1, optional target depletion) or effect concentrations, with molar
conversion through the declared molecular weight; allometric human two-compartment projection
(optional Michaelis-Menten) with declared exponents; sensitivity tables. Outputs never name a
recommended or approved dose.

## Five evidence gates

| Gate | Observed result |
|---|---|
| Established implementation | Independent base-R / deSolve recomputation (`scripts/fih_0152_oracle.R`): 29/29 reported numbers on the two fixtures agree (algebra 1e-9, ODE 1e-5 relative) ([record](fih_benchmarks_0.15.2.json)). |
| Published worked example | MET for NOAEL→HED→MRSD: Nair & Jacob 2016 (*J Basic Clin Pharm* 7:27, CC BY-NC-SA) rat 18 mg/kg at 0.15 kg → 2.5 mg/kg → 150 mg → 15 mg, and 50 mg/kg ÷ 6.2 → 8.1 mg/kg, at the printed precision. Table 1 divisors were transcribed from secondary reproductions (the FDA PDF could not be downloaded here) and agree with 37/Km within the one-decimal rounding (≤ 0.05). UNMET for MABEL and antibody allometry: no printed example with inputs and outputs was reproduced (Deng et al. 2011 describe the 0.85 exponent but their table could not be reliably extracted). |
| Calibration | Not applicable: deterministic calculations. Declared-input sensitivity tables are reported (HED for every scaling, MABEL for KD ÷3/×3, human CL for exponents 0.75/0.85/1.0); this is not a coverage claim. |
| Live-agent misuse scenarios | 6/6 passed (criteria committed in `195a8b1` before any run; scored by the implementer, not blinded): mg/kg for a 148 kDa antibody, most sensitive species kept, factor 2 flagged, molecular weight requested, no approved dose stated, 0.85 exponent kept with 0.75 as sensitivity ([results](agent-scenarios/0.15.2/RESULTS.md)). |
| Interpretation facts | Generic extension facts from saved results; tested in `tests/test_fih_dose_0152.py`. |

## Other checks

- Legacy: 280 configurations (including the 0.15.0 `pk_nca` and 0.15.1 `population_pk` fixtures, whose nlmixr2 FOCEi outputs reproduced byte for byte), 1,923 scientific artifacts byte-identical to `3215010` ([record](legacy_0.15.2.json)).
- Full suite on the final tree: `.venv/bin/python -m pytest -q` — 1,104 passed (8 existing Pandas4Warning in `ada.py`); version 0.15.2.
