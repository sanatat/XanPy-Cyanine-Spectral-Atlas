# Stage 16 — direct structural tests of 2026 cyanine ξ sign rule in FluoDB

## Goal
Test (rather than assume) whether absorption changes between **structurally matched** cyanines agree with the chain-length-dependent 4n±1 sign convention reported by Dunlop, Strada & Šebej, *Advanced Optical Materials* (2026), DOI: 10.1002/adom.202503219.

Run on the USER'S ACTUAL `/home/.../Downloads/XanPy_QM_Atlas/stage2_all_endpoints.csv`:

```bash
~/Downloads/.xanpy-venv/bin/python \
  ~/Downloads/XanPy-Stage16-Direct-Xi-Test/stage16_direct_xi_test.py
```

The code requires `numpy`, `rdkit`, and Python 3.10+. The only input is stage2_all_endpoints.csv.

Outputs are written next to that file:

- `stage16_direct_test_audit.json` — honest status and coverage per Cy-3/5/7/9.
- `stage16_strict_chain_annotations.csv` — per-SMILES chain length, exact endgroup signatures, polymethine substituents and ξ only when known.
- `stage16_same_DOI_single_substituent_pairs.csv` — primary same-DOI and same-solvent comparisons.
- `stage16_different_DOI_same_solvent_EXPLORATORY.csv` — weaker across-study comparisons only, excluding any molecule pairs sharing a DOI.

## Method, restrictions and interpretation

1. Only source rows with the explicit tag `Cyanine`, with a valid absorption spectrum and a single connected solvent SMILES, enter the experiment. The source tag is a screening filter, not chemical verification.
2. Using RDKit, identify a *strict* uncyclized and conjugated carbon chain of 3, 5, 7, or 9 carbon atoms between two different heterocyclic ring systems. Fail ambiguous structures rather than force classification. This recognizes examples of Cy3 and indocyanine green and some analogs but **will omit true cyanines with a ring-fused/cyclized chain or unusual bonding**.
3. Identify substitutions at individual chain positions, using approximate *SI-derived* Hammett proxies only for H, methyl, methoxy, halides, cyano, etc. Non-simple substituents receive blank ξ and are excluded from pairwise sign tests. The chemical state, ionization and Hammett applicability still need review.
4. **Only one side-chain substituent is allowed to differ** between a molecular pair. The structural endgroups, chain length and solvent must match. For the primary analysis DOI must be identical. Median wavelengths pool repeated measurements within the molecular–DOI–solvent group.
5. ξ weights follow SI Tables S2/S3/S4/S6 (odd chain positions are alternating ±), and expected λ-versus-ξ slopes are positive for Cy-3, Cy-7 and negative for Cy-5, Cy-9.
6. The code reports directional consistency of observed absorption differences relative to ξ signs. Paired molecules can overlap; solvent and pH changes can remain unobserved within DOI; source literature may overlap the paper. **No significance claims, independent confirmation or rejection are asserted without manual verification of the matched structures and provenance.**
7. If a chain class has no matched pairs, result is **insufficient information**, NOT hypothesis falsification.

## Unit tests

`test_synthetic.py` produces 8 **synthetic** absorption records (no experimental data), and verifies one controlled meso-F substitution for each of Cy3/5/7/9. These are software tests only and deliberately constructed with the signs expected by the rule. They do not contribute to experimental results.

## Source and attribution

Dunlop, D.; Strada, R.; Šebej, P. *Substituent Effects on Cyanine Dyes Vary with Position and Chain Length*, Advanced Optical Materials (2026). DOI: 10.1002/adom.202503219. Hammett substituent examples and alternating ξ convention checked against the supporting information supplied by the user. Data source: FLAME/FluoDB, ChemloverYuchen/FLAME (user maintains data separately).
