# 4PL reliability and relative ECx/ICx (0.13.4)

Use for prespecified concentrations achieving a fraction of the response change
from the low-dose plateau to the high-dose plateau. Optional config:

```json
"effect_levels": {
  "percentages": [10, 20, 50, 80, 90],
  "rationale": "Protocol requests these relative effect levels"
}
```

The names inherit `assay.endpoint`: EC50 gives EC10/EC20/EC50/EC80/EC90; IC50
gives IC10/IC20/IC50/IC80/IC90. Any distinct finite percentage strictly between
0 and 100 is accepted. Declare targets before examining their intervals.

For p = percentage/100 and h = abs(Hill),
`log10(Cp) = log10(C50) + log10(p/(1-p))/h`.
IC90 means 90% inhibition along the fitted response span, leaving 10% of the
span on a falling readout. An inhibition-percent readout can rise. The observed
readout direction never changes the meaning of the requested percent effect.
This convention differs from directly setting F=90 for a decreasing curve in
Prism's ECanything equation (that is 90% remaining response).

These are relative endpoints, not absolute 10%/90% killing or a fixed Y threshold.
An incomplete maximal response cannot be converted to 90% absolute inhibition
by relabelling a relative IC90. Bell curves have multiple crossings; 5PL and
bell-phase ECx/ICx are outside this new contract and are rejected explicitly.

## Intervals and withholding

Each Cp is reparameterized and profiled, reoptimizing Hill and plateaus using
the original bounds and objective. Scaling C50 interval ends by the estimated
Hill factor would omit Hill uncertainty and is not used. C50 itself reuses its
unchanged interval. Intervals are pointwise across endpoints, not simultaneous.
Each endpoint must lie in the measured positive-dose range and have a closed
profile without nuisance bounds when intervals are requested. A reportable
EC50 does not make EC90 reportable. An invalid parent fit withholds every Cp.
`uncertainty.method=none` retains its explicit no-interval meaning.

`effect_endpoints.csv` and `results.json#/effect_endpoints` preserve status and
diagnostics. Saved numeric candidates on withheld rows are audit values;
interpretation facts and the main HTML table do not promote them. Summaries
average log Cp within an experiment before combining independent experiments;
any withheld curve withholds that sample/endpoint summary. They use the existing
declared-independence and comparable-conditions gates.

## Default 4PL adequacy gate

The 4PL fitter and original intervals are unchanged. New diagnostics use:

- Replicate pure-error lack-of-fit F, alpha .005, for unweighted fits with
  positive pure-error and lack-of-fit degrees of freedom. Nonlinear F is
  approximate. Relative weighting does not use this homoscedastic F test.
- Lower-tail runs calculation on residual signs of ordered dose means, alpha
  .005. Exact sign-count combinatorics do not make fitted residuals independent;
  use the registered calibration scope rather than an exact-test claim.

Either failure withholds the curve, downstream RP and requested endpoints.
Nonestimable checks are saved explicitly. Passing or ineligible checks do not
prove 4PL adequacy or exclude a hook outside the measured doses. Do not switch
to 5PL/bell-shaped automatically or delete the high-dose decline. ELISA keeps
its separate calibration/QC contract unchanged in this release.

References: [GraphPad ECanything](https://www.graphpad.com/guides/prism/latest/curve-fitting/reg_ecanything.htm),
[0.13.4 observed evidence](../../../validation/RELEASE_0.13.4.md).
