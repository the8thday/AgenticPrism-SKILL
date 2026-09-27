# Interpretation guard fixtures (0.8.1)

Regenerate with `scripts/prepare_interpretation_fixtures.py`. All data are
synthetic derivatives of existing fixtures; source and output hashes are in
`provenance.json`. No measured biology or new inference method is represented.

- `config_boundary.json`: centered Gaussian residuals, seed 22, force the
  random-intercept variance boundary. Satterthwaite inference is withheld;
  fixed-effect audit means remain.
- `config_no_events_arm.json`: all mAb-high events are recoded to censoring to
  exercise no-event Cox withholding, an unreached median and a day-100 landmark
  beyond follow-up. The original censoring rationale is retained; the live
  agent must notice that the new early censoring times are unexplained.
- `config_dropout.json`: the existing tumor data at day 28. Three control
  animals are no longer measured. The original study rationale still names a
  day-21 protocol readout; the live agent must flag that conflict before using
  the day-28 result as a prespecified endpoint.

These deliberately inconsistent declarations are guard inputs. Do not copy
their rationale or endpoint definitions into a real experiment.
