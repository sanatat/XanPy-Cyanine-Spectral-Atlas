# Methods and reproducibility

## Original input and provenance

The primary analysis uses four **unmodified original compound ZIP archives** from Zenodo record 20119205: Cy5-7, Cy5-8, Cy7-3 and Cy7-8. Each archive must match the record's MD5 fingerprint (listed in `data/source_manifest.json`). Only absorption and emission text traces designated by the archived file's `Figure S...` label are parsed. Other files (e.g. spectroscopy supporting files) are not used.

The code selects two-column numeric wavelength/intensity arrays, removes nonnumeric rows and checks strict wavelength ordering. The published visible absorption/emission maxima restrict the search for the target band to ±30 nm and ±45 nm, respectively. This means the pipeline is **not** an unsupervised multi-band peak discovery method.

## Main peak location and checks

Wavelengths are sorted; raw-grid maxima and a locally fitted three-point quadratic vertex are evaluated using the raw trace and Savitzky–Golay windows of 7, 11, 21, 31, 51 points (third-order polynomial). The reported comparison uses the 11-point smoothed quadratic estimate. Variation over windows measures **algorithmic sensitivity only**, not measurement error or independent replication.

## Apparent half-height band-envelope width

Using the 11-point-smoothed **full** trace, locate the peak-adjacent left and right crossings at 50% of the exported peak height. Crossings are linearly interpolated to estimate the full width. It is a *descriptive half-height width of the exported envelope*; it is not necessarily a baseline-corrected, instrument-response-corrected or transition-deconvolved full width at half maximum. This distinction especially matters for cyanine absorption shoulders. Other fractions 40–80% are provided in the JSON output. Missing crossings produce null results rather than fabricated widths.

## Shape metrics and normalization

To compare band shape, profiles are peak-aligned on a relative wavelength grid and divided by their maximum over the shared analysis interval. Shape differences are quantified using relative mean absolute difference, root mean squared difference and intersection/union of the two normalized envelopes. These metrics are mathematical comparisons, not p-values or fluorescence efficiencies.

## Stokes shift

A peak-derived Stokes shift can be calculated as `1e7 * (1/absorption_nm - 1/emission_nm)` in `cm^-1`. This relates two maxima, not full excited-state reorganization energy.

## Model interpretation

The included Hückel workflow computes 16 qualitative gap sign checks for idealized polymethine chains. It is **not DFT** and does not perform a calibrated forward model of experimental spectra. Cy5/Cy7 pair comparisons provide qualitative context for the 2026 Dunlop–Strada–Šebej picture, not a four-length independent validation.

## Uncertainty limitations

One exported absorption/emission trace per compound/condition is not enough to determine experimental precision, reproducibility or a confidence interval. No hypothesis test or statistical significance of 0.4–1 nm peak shifts is claimed. Sample purity, isomer distributions, protonation, aggregation and residual instrumental baseline are not independently quantified by the released numerical analysis.
