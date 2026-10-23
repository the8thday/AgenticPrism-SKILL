# R workflow for custom PKNCA analyses

Since 0.15.0, single-dose NCA runs through `analysis_type: pk_nca` (see
[the contract](contract.md)). Use a saved R script only for designs that analysis
type does not cover (multiple dose, steady state, partial intervals, urine, sparse
sampling), and label it as a custom analysis outside the validated runtime.

Use an isolated project library, preferably with `renv`; do not change the user's
global R library. Check `Rscript --version`, `.libPaths()` and
`requireNamespace("PKNCA", quietly=TRUE)`. Record the installed version and
`citation("PKNCA")`. Missing dependencies must be reported, not replaced by a
Python approximation. If installation is part of the authorized task, use the
official stable distribution and save the resolved dependency versions.

The online PKNCA site can describe a development version. Check function
arguments, interval column names and options against **installed** package help
before using them; do not assert compatibility based only on a web vignette.

[assets/observed-profile-demo.R](../assets/observed-profile-demo.R) is a tiny
synthetic, positive-concentration single-profile example. It calls PKNCA for
observed-range Cmax/Tmax/AUClast only. It demonstrates API composition, not dose
handling, terminal extrapolation or the full reporting contract. It refuses an
existing output path and does not install packages. From the collection root:

```sh
Rscript --vanilla skills/pk-nca/assets/observed-profile-demo.R /tmp/pknca-demo-new
```

For a real study, build `PKNCAconc`, and `PKNCAdose` when dose is available,
with matching grouping keys. Explicitly pass the declared interval table and
study options to `PKNCAdata`, call `pk.nca`, and retain the object and full table.
First inspect the package's exclusion/reason columns before producing summaries.
Do not extend this demo to terminal endpoints merely by switching them on:
complete the terminal and dosing checks in [the contract](contract.md).

Primary references, consulted 2026-10-13:

- [PKNCA introduction](https://humanpred.github.io/pknca/articles/v01-introduction-and-usage.html)
- [Interval selection](https://humanpred.github.io/pknca/articles/Selection-of-Calculation-Intervals.html)
- [PKNCAdata API](https://humanpred.github.io/pknca/reference/PKNCAdata.html)
- [Calculation options](https://humanpred.github.io/pknca/reference/PKNCA.options.html)

On 2026-10-23 the demo was executed with PKNCA 0.12.1 (AUClast 25.056 matched a hand
calculation). A custom script is still not covered by the runtime's evidence.
