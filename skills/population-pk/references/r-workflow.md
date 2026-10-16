# R environment and first model

Check R and a project-local library for `nlmixr2`, `nlmixr2est`, `rxode2` and
`nlmixr2data`. nlmixr2/rxode2 require appropriate compilation tools even when R
itself is installed. Follow the official platform instructions; preserve install
and compilation errors rather than switching numerical engines. If authorized
to install, prefer stable packages and isolate the resolved dependency set,
for example with `renv`. Record R/package/compiler versions and citations.

The upstream 7.0 announcement documents estimator/default changes. Pin the
tested versions and explicit controls; a seed alone does not establish identical
results across package versions/platforms. Use installed help for `saemControl`,
`foceiControl`, event coding, diagnostic functions and fit-object accessors.

[assets/theophylline-demo.R](../assets/theophylline-demo.R) fits the upstream
`nlmixr2data::theo_sd` teaching data with an oral one-compartment model, diagonal
between-subject effects and additive residual error. Its starting values are
illustrative. It preserves the packaged events without recoding them and saves
the fit and console summary for inspection. From the collection root:

```sh
Rscript --vanilla skills/population-pk/assets/theophylline-demo.R /tmp/pop-pk-demo-new
```

This is an educational API template, not a validated analysis or complete study
report. No VPC, bootstrap, automatic reportability, two-compartment/infusion
demonstration or extraction-only facts are claimed by this script. In a real
analysis the agent must add the saved artifacts and diagnostic review required
by [the contract](contract.md). Capture units and dataset documentation before
interpreting any output. Oral-only CL and V are apparent values under F=1.

For subsequent models, write an explicit rxode2 ODE or documented `linCmt`
specification and verify the route/compartment mapping. Retain model rationale
and all attempted fits. Do not translate the model into a new Python estimator.

Primary references, consulted 2026-10-13:

- [nlmixr2 installation and overview](https://nlmixr2.org/)
- [PK model tutorial](https://nlmixr2.org/articles/running_nlmixr.html)
- [rxode2 datasets and events](https://nlmixr2.github.io/rxode2/articles/rxode2-event-types.html)
- [SAEM controls](https://nlmixr2.github.io/nlmixr2est/reference/saemControl.html)
- [nlmixr2 7.0 changes](https://blog.nlmixr2.org/blog/2026-08-05-nlmixr2-7/)

The local prototype only syntax-checks this script: the required R packages are
not installed in the checked environment. A future executed worked example
must save actual outputs and package versions before claiming numerical success.
