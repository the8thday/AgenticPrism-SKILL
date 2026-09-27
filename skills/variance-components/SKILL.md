---
name: variance-components
description: Estimate repeatability and intermediate precision using Gaussian REML random intercepts for declared nested, crossed or combined designs, including unbalanced observations and optional fixed effects. Report variance components, SD/CV and MLS/MOVER sum intervals and optional Satterthwaite intervals; identify unsupported method-validation claims.
---

# Variance components and precision

Locate the runtime using the [ordered procedure](../agentic-prism/references/runtime.md).
Read the [design and input contract](references/contract.md) before configuring.
Use this specialist for precision components; full bioanalytical validation,
accuracy, total error, concentration-dependent precision profiles and SOP
acceptance decisions are not implemented by this module.

## Establish the model

Ask for the response units, lot/run/replicate identities, nesting or crossing,
fixed adjustments, intended population of random levels and independent-error
assumptions. Several runs within each lot can distinguish lot and run variation;
one run per lot cannot. Reused local run IDs need the parent lot in their
random term. Technical replicates remain separate residual observations, not
independent lots. Do not silently average, delete incomplete rows or pool lots.

A random term is a list of factor columns. `[["lot"],["lot","run"]]` means
lot plus runs nested in lot. Adding `["analyst"]` specifies a crossed analyst
effect only if the experiment supports it. Fixed categorical levels and numeric
covariates adjust the mean; they are not random precision components.
Do not choose nesting, exclusions or fixed/random status to obtain a desired CV.

Require literal true applicability declarations and a substantive rationale.
If unspecified, clarify the design before analyzing. CV needs a declared positive
reference and rationale; do not silently substitute the observed mean. A small
or inappropriate denominator can make CV meaningless. No specification limit
is inferred from ICH/CLSI names or from a fixture.

## Run and interpret

Use `analyze --config ... --output NEW`, then `verify --run ...`; read saved
`interpretation_facts.json` and results. Lead with the declared design, component
estimates and intermediate-precision SD. Give repeatability separately, and CV
only with its denominator. Units of variance are squared response units.
Intermediate precision here sums every declared random component plus residual;
state which sources the sum includes, especially lot or analyst.

Satterthwaite intervals are marginal approximations. Report boundary flags and
missing intervals. A zero component estimate is not evidence of no variation;
its Satterthwaite interval can be unavailable. With MLS enabled, report the separate one-sided upper bound; never replace it by [0,0]. A numerically unrepresentable lower
endpoint also makes the whole interval unavailable; it is never a zero lower bound. The Satterthwaite total interval uses the active
component information matrix and can be optimistic at boundaries. No fixed-
effect tests or simultaneous coverage claim is provided. Review residual/model
assumptions from study knowledge; numerical convergence is not a normality test.

From 0.9.3 `sum_interval` defaults to `mls` (Satterthwaite covered only 83–85% in the 0.9.1 calibration; MLS 95.1–95.9% in 0.9.2). The applicability declaration is optional; designs whose mean-square strata are confounded are refused with a message to choose `satterthwaite` explicitly. Configs that need the VCA-matching Satterthwaite interval must say so.

Always disclose relevant calibration misses in facts and
[release evidence](../../validation/RELEASE_0.9.1.md). Do not turn agreement with
R into a general small-sample guarantee or call the entire method validated.
If asked whether results pass, request predeclared SOP limits and explain that
this module estimates precision only; it cannot establish accuracy or total error.

## Observed small-sample limits

In 1000 datasets per design (3 lots, 4 runs per lot, about 3 replicates),
total interval coverage was 85.1% ± 1.1261 percentage points MCSE for balanced
nested, 84.9% ± 1.1322 for unbalanced nested, and 83.4% ± 1.1766 for
nested/crossed with a fixed covariate. All miss the 93.6216% bound.
For true-zero lot variance, 0/387 available component intervals covered
zero (observed binomial MCSE 0); 612 boundary intervals and one numerically
unrepresentable interval were unavailable. Zero plug-in MCSE is not a confidence guarantee.
Do not promise nominal 95% coverage or absence of a variance source.

## 0.9.2 intervals

Balanced orthogonal designs use Graybill–Wang MLS on the actual mean squares.
Unbalanced designs and fixed adjustments use correlated, moment-matched
quadratic-form MOVER, with df/correlation estimated at REML. This is an
application of the published MOVER construction, not an exact unbalanced MLS
result or a VCA implementation. Random-term order is declared in the config.
Positive sums whose mean-square coefficients have mixed signs use signed MOVER.
REML point estimates remain separate from the unconstrained interval center;
at a boundary the estimate need not coincide with that center.

Lead with total SD and its selected interval, then repeatability and components.
Report a one-sided component upper bound as an upper bound, with its confidence
level. It does not replace an unavailable two-sided interval or prove a source
absent. Precision remains distinct from accuracy or method acceptance.

In 1000 simulations per row, total MLS/MOVER coverage was 95.1%, 95.3%, 95.9%,
95.1% and 95.5% for balanced three-lot, unbalanced, nested-crossed with a fixed
covariate, zero-lot and six-lot designs; all pass the 93.6216% bound. The
nested-crossed lot upper bound covers 92.8% (MCSE 0.8174 percentage points):
**miss**. Do not generalize total coverage to every component. Three-lot total
intervals can be extremely wide (mean variance width 169.86 for true total 6);
limited random levels still constrain precision. Details and every MCSE are in
[0.9.2 evidence](../../validation/RELEASE_0.9.2.md) and saved facts.

Sources: [Graybill–Wang 1980](https://doi.org/10.1080/01621459.1980.10477565),
[Burdick–Graybill 1984](https://www.stat.cmu.edu/technometrics/80-89/VOL-26-02/v2602131.pdf),
[Zou et al. 2009](https://publish.uwo.ca/~gzou2/10MOVERs.pdf).
These support the statistical constructions; no USP <1033> compliance claim is made.

Even in a balanced design, “exact mean squares” does NOT mean exact 95%
coverage of the MLS interval. State that MLS coverage is approximate. Match
calibration by design, component truth, fixed effects and boundary status;
a six-lot nonzero-lot simulation does not calibrate a six-lot zero-lot study.
Do not call an upper bound “small” without a meaningful declared scale or SOP.
