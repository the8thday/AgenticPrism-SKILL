# Post-hoc check of the 0.15.3 ordinal_or3 profile-coverage miss (not a registered check)

Registered row: `ordinal_or3 profile_coverage` 0.924 (868 replicates; bound [0.93, 0.97]) while
Wald coverage was 0.946. On 2026-10-25, 15 of the 19 replicates where the profile interval missed
and the Wald interval covered were regenerated with the registered seeds and refitted with R
`ordinal::clm` (`confint(type = "profile")`). In all 15 our log odds ratio equalled clm's and both
profile limits agreed within 1.4e-4 on the log scale. These are replicates with a large
overestimate (log OR about 3.0-3.5 against the true 1.10): the profile interval is asymmetric
and shifted upward, so it excludes the truth while the wider symmetric Wald interval includes it.
The miss is therefore a small-sample property of profile intervals in this design (10 subjects per
group, 5 categories), not an implementation error. The registered miss stands; this note does not
replace it.
