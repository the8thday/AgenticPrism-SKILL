# Enzyme-inhibition input contract (0.14.2; Morrison 0.14.3)

```json
{"analysis_type": "enzyme_inhibition", "input": "rates.csv", "source": "...",
 "design": {"enzyme": "Recombinant kinase domain", "substrate_unit": "uM", "inhibitor_unit": "nM", "rate_unit": "nmol/min/mg",
            "initial_rate_rationale": "Linear product formation over 10 min, under 10% conversion",
            "independent_unit": "experiment (separate enzyme dilution and day)", "independence_source": "Notebook ELN-123"},
 "enzyme_concentration": {"value": 0.2, "unit": "nM"},
 "mechanism": {"declared": "competitive", "rationale": "ATP-site inhibitor; ATP is the varied substrate"},
 "uncertainty": {"level": 0.95}}
```

`enzyme_concentration` is `{value, unit}` (molar: M, mM, uM, nM, pM; the inhibitor unit
must then be molar too) or `{rationale}` when unknown.

Columns: `experiment`, `substrate`, `inhibitor`, `rate`; replicate rows allowed; optional
`exclude` with `exclusion_reason`. Each experiment needs inhibitor 0, at least 3 inhibitor
concentrations and 4 nonzero substrate concentrations.

| Mechanism | Rate equation |
|---|---|
| competitive | Vmax S / (Km (1 + I/Ki) + S) |
| noncompetitive | Vmax S / ((Km + S)(1 + I/Ki)) |
| uncompetitive | Vmax S / (Km + S (1 + I/Ki)) |
| mixed | Vmax S / (Km (1 + I/Ki) + S (1 + I/(alpha Ki))); Kiu = alpha Ki |

Global unweighted least squares within each experiment (log parameters, multistart);
profile-F intervals; an open profile withholds the parameter. Model comparison: each
model's SS, AICc and F against mixed (descriptive; competitive and uncompetitive are
boundary cases). Withheld Ki: declared mechanism rejected against mixed (p < 0.05);
Ki < 10 × total enzyme; rank-deficient fit. Flags: highest substrate below 2 Km, lowest
substrate above Km, highest inhibitor below Ki. Summaries: geometric mean and t interval
over at least 2 experiments, withheld if any experiment's value is withheld.

## Tight binding: Morrison (0.14.3)

```json
{"analysis_type": "enzyme_inhibition", "input": "titration.csv", "source": "...",
 "design": {"enzyme": "Protease catalytic domain", "substrate_unit": "uM", "inhibitor_unit": "nM", "rate_unit": "RFU/min",
            "initial_rate_rationale": "Linear over 15 min, under 5% conversion",
            "independent_unit": "experiment (separate enzyme dilution and day)", "independence_source": "ELN-456"},
 "enzyme_concentration": {"value": 1, "unit": "nM"},
 "mechanism": {"declared": "competitive", "rationale": "Active-site inhibitor (crystal structure)"},
 "tight_binding": {"model": "morrison", "enzyme": "estimate", "rationale": "Ki expected near the enzyme concentration",
                   "km": {"value": 10, "source": "Substrate titration, report ENZ-12"}}}
```

Same columns. Each experiment: one nonzero substrate concentration, inhibitor 0 and at
least 6 nonzero concentrations (spanning the IC50 = Ki_app + Et/2). The total enzyme is
required (molar units, inhibitor unit molar).

v = v0 (1 − [EI]/Et), [EI] = ((Et + I + Ki_app) − √((Et + I + Ki_app)² − 4 Et I))/2.
Fitted: log v0, log Ki_app, and log Et when `enzyme: estimate` (reported as `Et_active`).
Profile-F intervals; open profile or rank-deficient fit withholds. Ki: competitive
Ki_app/(1 + S/Km), uncompetitive Ki_app/(1 + Km/S) (both need `km`), noncompetitive
Ki_app; mixed: Ki_app only. Flags: highest inhibitor below 2 × IC50, lowest above IC50,
Ki_app above 100 × Et (depletion negligible), below Et/100 (titration regime), estimated
active enzyme outside 0.5–2 × declared. Summaries as above; Ki_app summary only when the
substrate concentration is shared.
