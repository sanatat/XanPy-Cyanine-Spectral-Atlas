# Research portfolio summary (for a prospective supervisor)

**Project:** *Beyond Absorption Maxima: Reproducible Spectral-Shape Analysis of Substituent Effects in Cyanine Dyes*

> We independently reprocessed eight archival absorption and emission spectra from four cyanine dyes reported by Bedocchi, Okoročenkova and Klán (*J. Org. Chem.*, 2026). Their original numerical ZIP archives were verified against the Zenodo checksums. A Cy5 pair (OMe to Cl at C1′) exhibited nearly coincident absorption maxima (an interpolated shift of −0.38 nm), while its apparent absorption band-envelope half-height width changed from 85.28 to 51.83 nm. The analogous comparison for Cy7 (OMe to NHAc) gave an absorption blue shift of −4.61 nm and emission blue shift of −17.11 nm. These findings highlight spectral-envelope changes that are not captured by reporting absorption maxima alone. We also reproduced an idealized 16-case Hückel gap-sign pattern motivated by Dunlop, Strada and Šebej (2026), but do not claim a formal independent validation of their general substituent rule. The analysis includes reproducible parsing, provenance checks, numerical robustness tests, openly documented exclusions and quantitative visualizations.

## What this project demonstrates

- Direct handling of original exported spectroscopy data and preservation of source provenance.
- Reproducible Python analysis of peak shifts, apparent band widths and shape metrics.
- Careful treatment of 1 nm absorption sampling, interpolation, single-scan limitations and noncausal structural comparisons.
- A conceptual bridge from substituent chemistry and simple Hückel models to measured spectra without overstating experimental validation.

## Scope boundaries

**Not** newly performed synthesis, fluorescence spectroscopy, DFT/TD-DFT, determination of molecular excited-state lifetimes, or proof of a general Hammett correlation. The numerical spectra and reported molecular structures are from the cited experimental authors.

## References

- Bedocchi, O.; Okoročenkova, J.; Klán, P. *J. Org. Chem.* (2026). DOI: 10.1021/acs.joc.6c00573.
- Raw dataset: https://doi.org/10.5281/zenodo.20119205.
- Dunlop, D.; Strada, R.; Šebej, P. *Adv. Optical Mater.* (2026). DOI: 10.1002/adom.202503219.
