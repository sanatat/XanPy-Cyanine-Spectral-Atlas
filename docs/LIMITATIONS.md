# Limits and responsible claims

1. All raw observations originated in a third-party article, **not** a XanPy laboratory experiment. GitHub figures are reanalyses of the published dataset.
2. The published-paper values and local maximum extraction may differ by 0–1 nm owing to rounding, grid position, and export/processing differences. These are not experimental replicates.
3. All absorption traces are sampled on a ~1 nm wavelength grid. Sub-grid interpolation is an **algorithmic estimate**, not proof of sub-nm instrument accuracy.
4. `apparent_halfmax_bandwidth_nm` is a half-height width of the original exported, potentially overlapped band envelope. It is neither a deconvolved intrinsic FWHM nor corrected for independently measured baseline/instrument profile.
5. Comparing normalized intensities cannot measure fluorescence quantum yield, extinction coefficient, or molecular brightness.
6. Single scans cannot establish statistical significance or an experimental confidence interval; smoothing-window variation is *not* an uncertainty interval.
7. Solvent-matched C1′ pairs do not control all possible molecular and preparation variables. Interpret effects as **substituent-associated** differences rather than fully isolated causal effects.
8. `workflows/fluodb/` includes original project code and summary of user-reported results, but does not contain the original ~50k-row data snapshot. It is an **archived supporting workflow**, not a guaranteed byte-for-byte reproduction of those runs.
9. The idealized Hückel results are mathematical/theoretical sign checks only. Neither the SI-only same-paper regression nor the two published C1′ pairs count as statistical external confirmation of the general 2026 rule in FluoDB.
10. Source ZIP archives are not bundled because the Zenodo dataset's rights field did not show an unambiguous license at package assembly. Researchers should retrieve them from the dataset record and consult its applicable terms.
