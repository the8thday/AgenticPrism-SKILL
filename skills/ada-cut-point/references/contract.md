# ADA input and numerical contract — 0.9.1

Use `analysis_type: "ada_cut_point"`. Copy a matching config from
`fixtures/ada_synthetic/` only as a schema, then set every declaration from the
user's experiment. The authoritative options are `ada.DEFAULTS`; unknown keys
are errors. Boolean gates require JSON `true`.

## Design and CSV

A complete crossed panel: at least 50 independent subjects, each measured in
at least 3 runs, the same number (at least 2) of technical replicates per cell
remaining after explicit exclusions. These are implementation gates, not a
proof of sample-size adequacy. One representative negative population and one
reagent lot. Current decomposition is subject + run + residual; analyst/day/
plate are not independently estimated. More than one source of variation
confounded with run cannot be interpreted separately.

Required CSV columns:

| Column | Meaning |
|---|---|
| observation_id | Unique well ID; strings including NA/001 retained |
| subject_id, run_id, replicate_id | Explicit design IDs, unique triple |
| signal | Positive finite uninhibited signal |
| nc_signal | Positive declared arithmetic aggregate NC for that run; same value on all its rows. Retain NC raw-data provenance upstream; no invented NC |
| inhibited_signal | Positive paired drug-spiked signal for confirmatory; blank otherwise |
| lot, population, treatment_status | One reagent lot, configured population, literal drug_naive |
| exclude, exclusion_reason | Literal true/false and nonempty reason for exclusions |

The public fixture's NC=1 and reagent-lot label are declared artificial metadata
on simulated data, not examples for filling unknown assay metadata.

## Config

`assay`: tier (`screening`, `confirmatory`, `titer`), signal_unit, population,
rationale; literal-true drug_naive, representative_negative_panel,
independent_subjects, balanced_crossed_design and single_reagent_lot.

`cut_point`:

- false_positive_rate and target_rationale are required. The target is not
  silently assigned from a regulatory label. 0.05/0.01/0.001 in the fixtures
  are declared examples, not universal defaults.
- method: parametric (`mean + qnorm(1-FPR)*sample SD`) or nonparametric
  (R type-7 percentile). Computed over equally weighted subject/run cell means.
- normalization: none, ratio (signal/NC), difference (signal-NC).
- transformation: none, log (natural), boxcox with prespecified boxcox_lambda.
  No automatic lambda fitting or positive offset. Thresholds are back-transformed.
- Confirmatory uses `100*(1-mean(inhibited)/mean(uninhibited))`; negative
  inhibition is retained. Only no normalization/no transformation supported.
- selection: auto, fixed, floating, dynamic. Specified choices must match the
  diagnostics or reporting is withheld. Default comparison_alpha=.05;
  nc_correlation_min=.7. Values must be fixed with the protocol rationale.
- A significant normalized run mean or median-Levene variance difference
  selects dynamic and withholds common deployment. Otherwise unnormalized
  data select fixed. Normalized data select floating only with demonstrated
  raw shift and NC tracking. No benefit demonstrated -> withhold floating.
- normality_alpha=.05, skewness_limit=1. Parametric reporting pauses if BOTH
  the pooled Shapiro diagnostic fails (or is unavailable above 5000 cells) AND
  absolute skewness exceeds the limit. No automatic fallback to nonparametric.
- variance_confidence=.95 gives conservative simultaneous mean-square-based
  intervals for the additive crossed random ANOVA components.

`outliers`: literal true prespecified, rationale, analytical none/iqr_flag,
biological none/iqr_flag/iqr_exclude, iqr_multiplier (default 3).
Analytical flags use within-cell raw signal median residuals and run-specific
Tukey fences. They withhold pending investigation; no automatic deletion.
Biological fences act once on subject means on the configured analysis scale;
exclusion removes the entire subject. No refitting fences until values pass.
Explicit CSV exclusions are also audited and must leave the supported design.

`report.plot_style`: prism_like or standard. Rendering verifies saved hashes
and changes only figures/report/render_manifest, never thresholds or facts.

## Unimplemented boundaries

Dynamic future-run estimation, lower cut-point confidence bounds, fitted
Box-Cox lambda, robust-MAD alternatives, automatic analytical exclusion,
unbalanced/missing designs, separately identified analyst/day/plate REML,
confirmatory floating/log-ratio methods, sensitivity and drug tolerance are
not implemented. Extreme nonparametric tails are withheld when fewer than
one independent subject per run is expected beyond the percentile. Parametric
extreme tails are flagged; model-based extrapolation needs scientific review.

## Calibration interpretation

The point calculations are unchanged from 0.9.0. Use the exact conditional-FPR
record `validation/ada_calibration_0.9.1.json`; it supersedes the one-Bernoulli
estimate and the prior floating/IQR failure wording. The titer row misses the
new prespecified bound. Future drift is a separate demonstration. Read MCSE,
reportability denominator and the share exceeding twice target in saved facts.
No lower or upper confidence bound on a cut point is implemented.

## Lower confidence bound: 0.9.2

`cut_point.bound`: `point` (legacy default) or `lower`.
`bound_confidence`: numeric in (.5,1), default .90; state it in the protocol.
`bound_applicable`: literal true and `bound_rationale`: nonempty are required
for lower mode; the rationale must justify the marginal future subject/run
population, transformation and independence. No IQR handling, explicit exclusions
or dynamic deployment in lower mode. Existing point-mode numerics are unchanged.
`independent_pairs`: list of `[subject_id, run_id]` strings for nonparametric
lower mode, >=3 pairs, no repeated subject or run, all present in the panel.
Empty for parametric lower mode. Choose pairs before examining responses.

Parametric: total variance = MS_subject/b + MS_run/a +
(1-1/a-1/b)MS_residual. Its Satterthwaite df uses the independent MS variances.
Grand-mean variance uses clipped ANOVA components divided by a, b and ab.
Effective n = total / grand-mean variance. Lower percentile bound = mean +
sqrt(grand-mean variance) * noncentral-t quantile at 1-confidence, with
noncentrality z_percentile * sqrt(effective n). This moment approximation is
not exact in the crossed model. Raw/component estimates and effective df are saved.
Nonparametric: largest k satisfying P(Binomial(n,percentile)>=k)>=confidence,
then the kth ordered independent response; attained continuous confidence is saved.

The target is FPR AT LEAST the declared rate with confidence, not FPR at most.
Run mean differences do not automatically withhold this marginal random-run
estimand; point-mode selection diagnostics are retained for audit. Floating
still requires NC tracking. Neither method promises confidence for each realized
future run. The bound and its direction appear in interpretation_facts.json.
See validation/ada_bounds_calibration_0.9.2.json for all rows and MCSE, including
the high mean FPR of six-pair nonparametric titer bounds. Earlier statements that
no bound exists describe point mode and releases through 0.9.1 only.
