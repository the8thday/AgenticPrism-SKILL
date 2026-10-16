# Consultation behavior fixtures

Hand-constructed synthetic data for dialogue/routing checks; not biological evidence or a statistical calibration.

- `paired.csv`: six declared donors, one control and treated value per donor.
- `wells.csv`: raw IDs repeated across conditions and technical-looking rows; meaning of the IDs/wells must come from the scenario. Treatment and plate coincide.
- `limited.csv`: one literal `<LLOQ` qualifier, no numerical quantification limit supplied.

Prompts and evaluation criteria are recorded separately under validation/consultation_20261011.

These CSVs retain the CRLF bytes used during evaluation. The scoped attributes
preserve those bytes and recognize CRLF for whitespace checks; input hashes
were not changed after testing.
