# XanPy-QM Atlas — solvent effects and structural generalization in fluorophore photophysics

**Research-stage reproducible analysis | v0.9 | 1 November 2026**

An exploratory, CPU-friendly Python workflow to curate published fluorophore measurements, compare **structure-only** vs **structure + solvent** Random Forest models under molecule- and Murcko-scaffold-disjoint holdouts, and quantify solvent-associated absorption/emission/Stokes changes with study-aware sensitivity analyses. A focused structural screen identifies **candidate** O-/Si-/Ge-bridged xanthene-like cores; it is **not** chemical identity verification.

> **Result interpretation:** Across a user-run FluoDB snapshot, the addition of solvent-SMILES descriptors improved held-out emission prediction in three scaffold splits (MAE gains of 1.48, 1.45, and 0.54 nm). A full-database, molecule-balanced solvent comparison found emission changed more than absorption in 64.0% of comparisons (95% cluster-bootstrap CI 62.7–65.2%). **The stricter candidate O-conjugated-xanthene subset did not give a conclusive result**: 37 molecules from 7 DOIs; molecule-bootstrap fraction 59.1% (95% CI 45.7–72.3%). Do **not** present the full-database effect as confirmed for xanthenes.

## Data and attribution

- **FluoDB / FLAME**: Zhu et al., *A modular artificial intelligence framework to facilitate fluorophore design*, Nature Communications (2025), DOI [10.1038/s41467-025-58881-5](https://doi.org/10.1038/s41467-025-58881-5). Original CSV: https://github.com/ChemloverYuchen/FLAME/blob/main/FLAME/dataprocess/FluoDB.csv
- **QMugs**: Isert et al., *QMugs: quantum mechanical properties of drug-like molecules*, Scientific Data (2022), DOI [10.1038/s41597-022-01390-7](https://doi.org/10.1038/s41597-022-01390-7). No QMugs-derived predictions are claimed in these Stage 2–13 results.
- **Motivating xanthene study**: Dunlop et al., *Central Ring Puckering Enhances the Stokes Shift of Xanthene Dyes* (2024), DOI [10.1002/chem.202400024](https://doi.org/10.1002/chem.202400024). We do not attribute its mechanistic findings to FluoDB-derived observational results.

The authors' dataset is **not redistributed** in this repository. The file `results/user_reported/summary.json` transcribes selected **user-provided terminal results**; the full original CSV and full run logs are **not present in this release**, so those statistics were *not independently recomputed here*. The provenance hash recorded in that file is for the specific user-run CSV snapshot, not necessarily the current GitHub version.

## Reproduce (Linux Mint / Ubuntu)

```bash
cd ~/Downloads
python3 -m venv .xanpy-venv
./.xanpy-venv/bin/python -m pip install -r /path/to/XanPy-QM-Atlas-v0.9/requirements.txt
curl -fL --retry 3 -o FluoDB.csv https://raw.githubusercontent.com/ChemloverYuchen/FLAME/main/FLAME/dataprocess/FluoDB.csv
# Optional: pin source CSV with original SHA-256 3dde9fc378363fc2344811f4381169bf56bfc754e0374a632390550be0acc954
export XANPY_ROOT="$HOME/Downloads/XanPy_QM_Atlas"
mkdir -p "$XANPY_ROOT"
cd /path/to/XanPy-QM-Atlas-v0.9
./scripts/run_all.sh
```

`run_all.sh` uses the Python interpreter from `~/Downloads/.xanpy-venv` and executes steps 2–13 sequentially (the 3-seed RF re-fits may take time on low-memory machines). To run focused steps manually:

```bash
~/Downloads/.xanpy-venv/bin/python scripts/stage2.py
~/Downloads/.xanpy-venv/bin/python scripts/stage3.py emission_nm
~/Downloads/.xanpy-venv/bin/python scripts/stage3.py absorption_nm
~/Downloads/.xanpy-venv/bin/python scripts/stage4.py
~/Downloads/.xanpy-venv/bin/python scripts/stage5.py emission_nm
~/Downloads/.xanpy-venv/bin/python scripts/stage6.py
~/Downloads/.xanpy-venv/bin/python scripts/stage7.py
~/Downloads/.xanpy-venv/bin/python scripts/stage8.py
~/Downloads/.xanpy-venv/bin/python scripts/stage9.py
~/Downloads/.xanpy-venv/bin/python scripts/stage10.py
~/Downloads/.xanpy-venv/bin/python scripts/stage11.py
~/Downloads/.xanpy-venv/bin/python scripts/stage12.py
~/Downloads/.xanpy-venv/bin/python scripts/stage13.py
```

To use alternate locations, set `XANPY_ROOT` and `FLUODB_CSV`. Use `python -m pytest -q tests` for smoke checks. **Script variants consolidated from iterative research prototypes; minor output/descriptor-order differences from the user's original scripts are possible**. Metrics should be regenerated before claiming exact numerical reproduction. For original previously run outputs, keep the user's Stage 2–13 scripts and audit files.

## Findings and cautions

See [RESULTS.md](docs/RESULTS.md), [METHODS.md](docs/METHODS.md) and [LIMITATIONS.md](docs/LIMITATIONS.md).

- **Not a validated quantum model.** No new DFT, TD-DFT, HOMO/LUMO, solvent polarity parameters, or optical energy predictions were computed or validated by this Stage 2–13 workflow. Structural solvent-SMILES descriptors are proxies only.
- **Observed association is not causal.** Same DOI ≠ same pH, protomer, temperature, concentration or measurement protocol.
- **Multiple solvent contrasts for one study are dependent.** Inference is sensitivity-tested at molecule and DOI levels; those bootstraps are not external validation.
- **Candidate dye classes not manually verified.** The strict reference tests validate eight hand-selected patterns, not the identity of 812 candidate structures.
- **Small subgroup.** In the O-conjugated candidate subset, the molecule-bootstrap 95% CI includes the null fraction of 0.5. No robust superiority is established for this subset; Si/Ge comparisons have insufficient DOI-matched multisolvent data.
- **Data source tags can be nonchemical:** `6n6` is a broad architecture label, not a specific dye family.

## Folder guide

- `scripts/stage2.py`–`stage13.py`: sequential original-analysis-style scripts, with targeted compatibility edits.
- `results/user_reported/summary.json`: compact source-transcribed statistics (not synthetic per-row data).
- `docs/`: reproducibility, results, and limitations.
- `data/`: empty on purpose; place downloaded data locally and do not commit unreviewed copyrighted datasets.

## Research status

**Exploratory independent analysis of a public compiled dataset.** Suitable as a carefully qualified portfolio/GitHub research demonstration; not a published or externally validated mechanism, not evidence of a novel dye being synthesized, and not proof of DFT accuracy.
