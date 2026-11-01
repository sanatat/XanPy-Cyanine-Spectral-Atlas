# Project history: FluoDB → Hückel → archived spectra

These stages are scientifically related but are not all parts of the current four-spectrum main analysis. Distinguish **past user-reported terminal statistics** from reproducible outputs included in this repository.

| Track | Earlier finding | Scientific status | Source in this repo |
|---|---|---|---|
| FluoDB curation (stages 2–5) | 49,602 cleaned measurement rows; single scaffold-split emission MAE 39.84 nm (structure) versus 38.39 nm (structure+solvent) | Historical user-run result; full input snapshot not included | `workflows/fluodb` |
| Solvent contrasts (stages 6–9) | Molecule-balanced emission-change-magnitude greater than absorption in approx. 64.0% of contrasts; cluster-bootstrap 95% CI approx. 62.7–65.2% | Observational, not solvent causality | `workflows/fluodb` |
| Candidate O-xanthene subgroup (stages 10–13) | 812 conjugated xanthene-like candidate structures; only 37 DOI-matched multisolvated molecule groups from 7 articles; ~59.1% comparison fraction with CI crossing 50% | No conclusive xanthene-specific superiority | `workflows/fluodb` |
| Hückel theory | 16/16 signs and 96/96 perturbation-strength sign cases match the *theoretical model pattern* | Reproduction of idealized model, not novel experimental confirmation | `workflows/huckel` |
| Structural FluoDB cyanine screen (stages 14–17) | 189 accepted structures, 176 assigned ξ with strict rules in all-tag scan; zero qualifying same-DOI single-substitution comparisons detected | Cannot test general rule on this screened FluoDB snapshot; classifier may miss motifs | `workflows/cyanine_identity` |
| Archival spectrum analysis (stages 18–20) | Verified four published compound ZIPs, recovered eight spectra, analysed shape changes and robustness | Directly reproducible from original public records | This repository's main workflow |

### Recommended usage

For a GitHub-facing standalone demonstration or postgraduate application, emphasize what can be rerun here: *four original raw-archive hashes, eight numerical spectra, numerical peak-shape analysis, figures and tests*. Historical FluoDB claims should be labeled as user-run unless the original complete CSV and logs are furnished and verified.
