# Cyanine Hückel Sign-Rule Reproduction (2026 Perspective)

A self-contained numerical reproduction of the **4n+1 / 4n−1** position-dependent substituent sign pattern discussed in Dunlop, Strada & Šebej, *Advanced Optical Materials* (2026), DOI [10.1002/adom.202503219](https://doi.org/10.1002/adom.202503219).

## Verified outputs

1. A nearest-neighbor Hückel Hamiltonian is diagonalized for 3-, 5-, 7- and 9-site open polymethine chains. The even-electron model has N+1 electrons for odd N, so (N+1)/2 orbitals are occupied.
2. The model applies a **small on-site energy perturbation** to the meso site and adjacent α site. A positive perturbation models EDG destabilization; negative models EWG stabilization. This approximation is not a chemically calibrated Hammett σ or real substituent Hamiltonian.
3. HOMO and LUMO site densities give the first-order frontier-gap response: Δgap ≈ [|c_LUMO,i|² − |c_HOMO,i|²] × Δε_i.
4. Direct diagonalization of the perturbed Hamiltonian tests the same sign. A *larger* gap is interpreted qualitatively as a blue-shift in the simple picture, not a numerical photophysical excitation wavelength prediction.
5. The sixteen qualitative model predictions match the sign pattern in the 2026 Perspective. The four explicitly printed approximate literature examples (Cy5/7, meso OH/NO2/OEt/CN) also match in sign. **These four examples come from the same Perspective; they do not constitute independent experimental validation**.

## Execute

```bash
python3 -m pip install numpy
python3 scripts/huckel_test.py
python3 -m pytest -q tests/
```

For Linux Mint with externally-managed-environment (PEP 668), activate an existing venv:

```bash
~/Downloads/.xanpy-venv/bin/python scripts/huckel_test.py
```

Outputs: `results/huckel_16_sign_tests.csv`, `results/paper_four_examples_sign_check.csv`.

## Experimental Hammett/ξ validation: NOT completed

The article cites **Supporting Information tables S2–S4, S6**, with published ξ/σ_p values and absorption maxima. The supporting DOCX filename advertised on the publisher is `adom70951-sup-0001-SuppMat.docx`. The download endpoint returned HTTP 403 in this environment. **We did not invent SI data or present a regression based on fictitious measurements.**

To validate real SI measurements, enter verified values into `data/hammett_SI_template.csv` after reading those tables, preserving source citation, DOI and solvent. The ξ index must use the paper's alternating +1/-1 positional Hammett weighting (Figure 5), not unverified assumptions about a SMILES or raw FluoDB tag. Then run:

```bash
python3 scripts/hammett_validate.py data/hammett_SI_template.csv --out results/hammett_real_SI_results.json
```

The parser separately fits `absorption_nm ~ xi` for each `(chain_length, solvent)` and checks whether 4n−1 yields positive and 4n+1 yields negative slope. At least 4 distinct compounds and 3 distinct ξ values are required per fit. Splitting by solvent is conservative; cross-solvent pooled fits can be confounded. Stronger inference needs molecule-level clustering, errors/uncertainties and rigorous curation.

## Limitations

- These are **idealized isolated polymethine π-chains**, not fully substituted real cyanines with heteroaromatic terminal groups, solvent, charge redistribution, geometry relaxation, or excited-state response.
- Δε=±0.15 is a model perturbation in units of |β|, **not** a DFT result, measured orbital energy or calibration to optical nm.
- Comparison of predicted qualitative blue/red signs with four examples directly quoted from the 2026 Perspective is a consistency check, not independent experimental validation.
- Not a new discovery: the 2026 paper already established the rule. The value is reproducible numerical implementation, mechanistic checks and a route to testing robustness with curated experimental data.

## Provenance

Paper: David Dunlop, Rebecca Strada, Peter Šebej, *Substituent Effects on Cyanine Dyes Vary with Position and Chain Length*, Advanced Optical Materials 14 (2026) e03219, https://doi.org/10.1002/adom.202503219.

## Additional perturbation-strength sensitivity check

```bash
python3 scripts/robustness_and_figure.py
```

The exact-diagonalization predictions preserve all 16 direction labels at six perturbation strengths (0.01, 0.05, 0.1, 0.15, 0.25, 0.4 in units |β|), **96/96 matching the reference sign pattern**. This is mathematical model robustness to a single parameter, not validation against 96 new compounds. The script also creates `figures/huckel_gap_sign_inversion.png` if matplotlib is available.
