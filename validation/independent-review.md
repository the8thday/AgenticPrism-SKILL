# Independent forward review — 2026-09-22

Scope: independent behavioral evaluation of the AgenticPrism router, equilibrium-binding specialist and local runtime. No implementation or main test files were edited by the evaluator. Generated experimental artifacts were isolated outside the repository. This is a synthetic workflow review, not a Prism numerical-equivalence benchmark.

## Skill routing and usability

Read the router, registry, specialist, all specialist references and distribution README. Applied these natural-language requests:

| Request | Observed instruction-driven behavior | Finding |
|---|---|---|
| “ELISA plate unknown concentrations” | Identify calibration/interpolation purpose; explain that the ELISA module is planned; request standards, known concentrations, unknown wells, dilution and plate metadata as needed. Do not run the KD model on OD data or invent a completed analysis. | Appropriate boundary |
| “Endpoint table kon/koff” | Explain that endpoint observations do not identify association/dissociation rates; time-course data and a validated kinetic module are required. | Appropriate boundary |
| “Three independent KD experiments” | Use equilibrium-binding after assay and independence checks; fit each curve, collapse technical curves within experiment on log KD, summarize independent experiments with equal weights. Do not substitute a significance test. | Executed successfully |
| “Re-render saved report in standard style” | Resolve the distribution from the specialist folder, invoke its local executable with an absolute saved-run path, verify immutable outputs. | Executed from `/tmp` successfully |

The two-parent traversal from the specialist folder resolves the actual distribution root. Its `.venv/bin/agentic-prism` executable works without relying on the current working directory. A requested `analysis_type: elisa_quantification` exits 2 with an explicit unsupported-schema/module error rather than silently selecting KD.

A Markdown link scan initially found two links to `validation/README.md` unresolved while that file was still being prepared. All other local Markdown links in the skills and root README resolved. Release completion should recheck that validation index.

## Numerical and failure behavior

A separately generated 14-point synthetic curve used true KD = 8 nM, baseline = 12 RFU, amplitude = 150 RFU and fixed small residual perturbations.

- Estimated KD: **7.986861883 nM**; 95% profile-F interval: **7.767342964–8.212616510 nM**.
- An independent reference implementation solved baseline/amplitude with exact linear least squares for each fixed KD, then used scalar minimization and root finding: **7.986861847 nM**, interval **7.767342974–8.212616296 nM**.
- Converting the same concentrations from nM to µM changed KD by **9.16 × 10⁻¹⁰ relative**.
- Flat response, insufficient concentration support and nonpositive log responses produced explicit failed-curve artifacts and renderable reports.
- Negative concentration and absent assay applicability metadata were rejected.
- Three synthetic independent experiments with KD near 4, 8 and 16 nM produced a geometric mean of **7.989176027 nM**, log-t interval **1.440409087–44.31167102 nM**, and `summarized_log_t` status. The geometric mean was independently checked against the individual estimates.

## Bugs discovered and fixes rechecked

1. **Identifier preservation:** permitted `curve_id="NA"` was previously parsed as missing when loading predictions/residuals for rendering. The initial report silently had no model/residual line despite 161 saved predictions. After the correction, both rendered themes contain **161 prediction rows**, and the render command exits 0.
2. **Sensitivity status:** alternative fits intentionally skip CI computation but were previously labeled `limited` solely because no CI was computed. After the correction, both well-behaved sensitivity scenarios report **`fit_only`**, accurately distinguishing them from the primary estimate with interval. Alternate-loss KD = 8.135874499 nM; excluding highest concentration KD = 8.010702161 nM.

The updated renderer successfully re-rendered a saved run created before these fixes. Hashes of every manifest-listed scientific artifact remained identical, and the separate `verify` command passed.

## Visual inspection and limits

The generated Prism-like PNG was opened and visually inspected: observation points, fitted line, readable axes/units, separate zero-control panel and signed residuals were present and legible. This evaluator did not perform browser theme-switch/download interaction or mobile-layout testing. No biological assay applicability, Prism reference project, general coverage guarantee or operating-system portability was established by these checks.

Temporary reproduction scripts: `/tmp/agentic_prism_forwardtest.py`, `/tmp/agentic_prism_profile_reference.py`, `/tmp/agentic_prism_skill_review.py`. Latest skill-workflow artifacts: `/tmp/agentic_prism_skill_krpa6_z_/`. These temporary paths are execution evidence, not packaged release dependencies.
