# Scientific and release limitations

**Most important:** The broad-dataset emission-dominance trend does not statistically establish the same effect among the candidate O-bridged conjugated xanthene-like subset. The molecule-bootstrap interval includes 0.5, as does the interval for mean absolute-shift difference around zero. The 7-source DOI bootstrap is unstable for strong statistical claims.

## Source and measurement caveats

- Only DOI, solvent, canonical molecular representation and broad endpoints are tracked. Same DOI does not prove equal pH, concentration, aggregation, temperature, instrument, or species/protomer.
- Different publication/source databases can introduce extraction and measurement differences. Repeated molecule-solvent values were pooled across publications in the RF, but within-publication matching was used for observational comparison.
- Solvent descriptors based on atom counts, logP and TPSA are **not** experimental dielectric constants or Kamlet–Taft parameters. They cannot alone identify hydrogen bonding or charge-transfer mechanisms.
- Source data omit DOI for many records, notably in Epsilon, and the candidate screen is not human-verified. The eight chemical reference self-tests are necessary but far from sufficient.
- SMILES represent a chemical drawing, not guaranteed fluorescent chemical state. Protonation, tautomers, regioisomers, aggregation and excitation/emission assignment remain uncertain.
- Split-group overlaps across random seeds and correlated molecular families make naively pooling p-values incorrect; subgroup DOI count is only seven.

## Release/provenance caveats

- The repo contains executable analysis scripts and *transcribed user-reported summary results*, but **no original FluoDB CSV or full primary output files**. It must not assert that the displayed results were independently recomputed or exactly reproduced from the archived scripts.
- The scripts were consolidated from stepwise research drafts. In particular stage3 serialization and some descriptor ordering were normalized for compatibility; a fresh numerical run may differ from the historical experiment. Run with pinned datasets and dependencies, compare outputs and record the git commit hash before a manuscript-level claim.
- Any future publication should include a detailed license review for each imported dataset, source and extraction conditions, a chemical-identity spot audit, repeated paired-split inference with fixed hypotheses, and preferably new experimental validation.
