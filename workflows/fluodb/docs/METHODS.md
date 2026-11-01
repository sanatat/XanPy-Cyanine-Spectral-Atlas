# Methods (concise)

1. **Curation:** canonical isomeric SMILES, source row + source DOI, solvent as case-sensitive SMILES, labels absorption/emission (200–1500 nm), PLQY (0–1), molar extinction (0–1e8). Do not convert PLQY percentages without evidence. Flag negative Stokes values, do not use them in paired-stokes comparisons. Corrected v0.9 Stage 2 uses canonical solvent SMILES for its `solvent_key` rather than `.lower()`.
2. **Prediction:** RDKit Morgan ECFP4 128-bit fingerprints plus bulk molecular descriptors; solvent structural descriptors computed from SMILES. Random Forest 60 trees, min_samples_leaf=3, max_depth=18, fixed model seed 42. Holdout 20% **grouped by molecule** or **Murcko scaffold**. Target repeats for the same molecule–solvent are median-aggregated. This does not control DOI heterogeneity in targets.
3. **Model comparisons:** paired MAE improvement (structure-only MAE less structure+solvent MAE); clustered bootstrap at molecule or scaffold groups on a fixed split; 3 independently seeded scaffold holdouts for preliminary sensitivity. No cross-lab validation or hyperparameter search.
4. **Solvent association:** within (canonical structure, DOI) compare measured pairs from distinct solvents, median-combine repeats; calculate signed wavelength and wavenumber changes and Stokes shift `1e7(1/absorption_nm - 1/emission_nm)` in cm⁻¹; average over pair contrasts within molecule-DOI before molecule/DOI balanced bootstraps. The different units and signed-vs-absolute summaries must never be conflated.
5. **Family audit:** source tags are descriptive (e.g. `6n6` is broad topology), while later stage screen uses strict fused 6/6/6 arrangement and bridge atom type and distinguishes carbonyl/sp3/conjugated C9. All labels **candidate** until chemical validation.
6. **Focused xanthene-like audit:** Require DOI, paired nonnegative-Stokes absorption/emission, valid single-component solvent SMILES, >=2 distinct solvents per molecule–DOI. 37 groups/7 DOIs for current uploaded results.

## Not included

No new DFT/TD-DFT, no exact QMugs join on this full FluoDB snapshot, no SHAP causal analysis, no empirically measured solvent polarity/H-bond scales, no manual verification of all proposed xanthene-like structures. Ground-state QMugs HOMO–LUMO gaps cannot be treated as optical excitation energies.
