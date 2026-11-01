# Reported results — XanPy-QM Atlas, stages 2–13

_Transcribed from user-shared terminal output on 1 Nov 2026. The raw 6.9MB FluoDB file and 2GB QMugs file were not available in this release workspace; the original numerical runs were not independently rerun. Source rows remain attributed to the FluoDB dataset, and all results are exploratory._

## Scope, quality, predictive checks

- 49,851 FluoDB source rows; 49,602 with at least one valid optical endpoint; 30,839 unique structures. Source SHA-256: `3dde9fc378363fc2344811f4381169bf56bfc754e0374a632390550be0acc954`.
- 31,350 usable absorption values, 25,566 emission values, 19,777 PLQY values, 26,273 molar extinction values. Endpoints overlap; counts are not additive.
- Structure-only vs structure+solvent Random Forest baseline, **scaffold-disjoint** test: emission MAE **39.8383 → 38.3875 nm**, R² **0.6155 → 0.6345**; absorption MAE **42.6283 → 42.2921 nm**.
- Emission MAE gains from three independently seeded **split selections** (overlapping data; same model design): **1.4763, 1.4508, 0.5358 nm**; average **1.1543 nm**. This does not mean three independent experimental datasets.
- Fixed-split clustered bootstrap for scaffold test emission gain **1.4508 nm**, 95% CI **[0.5575, 2.4065]**; for absorption gain **0.3361 nm**, 95% CI **[−0.3922, 0.9979]**.

## Full-dataset within-DOI contrasts (not necessarily causal)

- 50,530 **dependent** solvent-pair contrasts, representing 3,124 molecules, 3,195 molecule–DOI units, and 753 distinct DOIs.
- Median |emission shift| **10 nm**, median |absorption shift| **4 nm**, median absolute Stokes shift change **443.43 cm⁻¹**.
- Molecule-balanced fraction for which |Δemission| > |Δabsorption|: **0.63969**; 95% CI **[0.62666, 0.65201]**. In cm⁻¹ the proportion was **0.60448**; 95% CI **[0.59119, 0.61749]**.
- For seven selected solvent pairs, the toluene → acetonitrile comparison spanned **891** molecules, with mean signed emission shift **+20.9084 nm**. Subgroups diverged: BODIPY (413 molecules) **−5.114 nm**, PAHs (100) **+26.155 nm**, `6n6` architecture label (69) **+69.754 nm**. These labels are source tags, not all verified structural classes.

## Strict structural screen and focused subset

- A topology and bridge-chemistry filter returned 1,026 **candidate O-bridged conjugated** records from 812 unique structures, 439 valid paired spectra and 383 DOI-traceable paired spectra. Classification is not verified dye identity.
- Within this focused subset: **37** molecules/DOI units across **7** distinct DOIs; **124** dependent within-study solvent pairs. No repeated-solvent disagreement over 20 nm was identified by the script's specified check.
- **Molecule-cluster bootstrap:** fraction with |Δemission| > |Δabsorption| in nm = **0.59099**, 95% CI **[0.45671, 0.72252]**. In wavenumbers = **0.55946**, 95% CI **[0.42878, 0.68921]**. Both include **0.5**.
- The molecule-balanced mean of |Δemission| − |Δabsorption| was **+1.0982 nm**, 95% CI **[−0.08502, 2.25412]**, including zero. DOI sensitivity analysis gave a different result but only **7 DOI clusters**, too few for strong inference.
- **Interpretation:** the full-dataset pattern does **not** establish stronger solvent-dependent emission shifts within the candidate xanthene-like subclass. The dataset is insufficient for a quantitative O/Si/Ge bridge comparison.

## Defensible conclusion

> In a curated published fluorophore dataset, molecular structure and solvent-derived descriptors jointly improved exploratory emission predictions under molecule- and scaffold-held-out tests. Study-aware analyses showed a larger emission than absorption response to solvent variation across the broad dataset, with substantial scaffold/source-tag heterogeneity. A stricter candidate O-bridged conjugated-xanthene subset had too few independently sourced, multisolvent observations to establish the same trend, motivating focused experimental validation and better chemical-state metadata.

This describes an **observational secondary-data study**, not a mechanistic/quantum validation or a new molecule design.
