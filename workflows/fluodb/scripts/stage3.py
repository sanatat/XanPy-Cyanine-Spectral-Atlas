#!/usr/bin/env python3
"""XanPy QM Atlas stage 3: leakage-aware, low-memory fluorescence baselines.

Usage: python xanpy_atlas_stage3.py [emission_nm|absorption_nm]
"""
import csv
import json
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, rdFingerprintGenerator, rdMolDescriptors
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupShuffleSplit

RDLogger.DisableLog('rdApp.error')

TARGET = sys.argv[1] if len(sys.argv) > 1 else 'emission_nm'
if TARGET not in ('emission_nm', 'absorption_nm'):
    raise SystemExit('Choose emission_nm or absorption_nm')

D = Path(os.environ.get('XANPY_ROOT', str(Path.home() / 'Downloads' / 'XanPy_QM_Atlas')))
D.mkdir(parents=True, exist_ok=True)
SOURCE = D / 'stage2_all_endpoints.csv'
if not SOURCE.is_file():
    raise SystemExit(f'Missing {SOURCE}')

FP = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=128)
MOL_NAMES = ['mw', 'logp', 'tpsa', 'hba', 'hbd', 'rotb', 'rings',
             'arom_rings', 'frac_sp3', 'heavy_atoms', 'formal_charge',
             'n_atoms', 'o_atoms', 'halogen_atoms', 's_atoms', 'si_atoms', 'ge_atoms']
SOL_NAMES = ['sol_mw', 'sol_logp', 'sol_tpsa', 'sol_hba', 'sol_hbd',
             'sol_arom_rings', 'sol_heavy_atoms', 'sol_formal_charge', 'sol_parse_failed']


def numeric(x):
    try:
        v = float(x)
        return v if np.isfinite(v) else None
    except (ValueError, TypeError):
        return None


def molecule_features(m):
    counts = Counter(a.GetSymbol() for a in m.GetAtoms())
    halogens = sum(counts[e] for e in ('F', 'Cl', 'Br', 'I'))
    vals = [Descriptors.MolWt(m), Crippen.MolLogP(m),
            rdMolDescriptors.CalcTPSA(m), rdMolDescriptors.CalcNumHBA(m),
            rdMolDescriptors.CalcNumHBD(m),
            rdMolDescriptors.CalcNumRotatableBonds(m),
            rdMolDescriptors.CalcNumRings(m),
            rdMolDescriptors.CalcNumAromaticRings(m),
            rdMolDescriptors.CalcFractionCSP3(m), m.GetNumHeavyAtoms(),
            Chem.GetFormalCharge(m), counts['N'], counts['O'], halogens,
            counts['S'], counts['Si'], counts['Ge']]
    bits = FP.GetFingerprintAsNumPy(m).astype(np.float32)
    return np.concatenate((np.asarray(vals, dtype=np.float32), bits))


def solvent_features(raw):
    # Source uses case-sensitive solvent SMILES; never lower-case these.
    mol = Chem.MolFromSmiles(raw.strip()) if raw.strip() else None
    if mol is None or len(Chem.GetMolFrags(mol)) != 1:
        return np.asarray([0.0] * 8 + [1.0], dtype=np.float32), False
    vals = [Descriptors.MolWt(mol), Crippen.MolLogP(mol),
            rdMolDescriptors.CalcTPSA(mol),
            rdMolDescriptors.CalcNumHBA(mol),
            rdMolDescriptors.CalcNumHBD(mol),
            rdMolDescriptors.CalcNumAromaticRings(mol),
            mol.GetNumHeavyAtoms(), Chem.GetFormalCharge(mol), 0.0]
    return np.asarray(vals, dtype=np.float32), True


# Median-aggregate repetitions only within an exact molecule + solvent pair.
measurements = defaultdict(list)
audit = Counter()
with SOURCE.open(newline='', encoding='utf-8') as file:
    for row in csv.DictReader(file):
        audit['stage2_rows_seen'] += 1
        y = numeric(row.get(TARGET))
        if y is None or not 200 <= y <= 1500:
            continue
        smiles = (row.get('canonical_smiles') or '').strip()
        solvent = (row.get('solvent') or '').strip()
        if not smiles:
            continue
        measurements[(smiles, solvent)].append(y)
        audit['valid_target_rows'] += 1

mol_cache, sol_cache = {}, {}
xs, xsol, yy, structure_ids, scaffolds, solvent_names = [], [], [], [], [], []
for (smiles, solvent), values in measurements.items():
    if smiles not in mol_cache:
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            audit['bad_molecule_groups'] += 1
            mol_cache[smiles] = None
        else:
            try:
                from rdkit.Chem.Scaffolds import MurckoScaffold
                scaffold = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m))
                mol_cache[smiles] = (molecule_features(m), scaffold)
            except Exception:
                mol_cache[smiles] = None
                audit['molecule_descriptor_errors'] += 1
    if mol_cache[smiles] is None:
        continue
    if solvent not in sol_cache:
        try:
            sol_cache[solvent] = solvent_features(solvent)
        except Exception:
            sol_cache[solvent] = (np.asarray([0.0]*8+[1.0], dtype=np.float32), False)
    mf, scaffold = mol_cache[smiles]
    sf, valid_sol = sol_cache[solvent]
    xs.append(mf)
    xsol.append(sf)
    yy.append(median(values))
    structure_ids.append(smiles)
    scaffolds.append(scaffold if scaffold else '__ACYCLIC__')
    solvent_names.append(solvent)
    if not valid_sol:
        audit['unparsed_solvent_pairs'] += 1

if len(yy) < 30:
    raise SystemExit('Not enough labeled molecule-solvent pairs for evaluation')

X = np.asarray(xs, dtype=np.float32)
XS = np.concatenate([X, np.asarray(xsol, dtype=np.float32)], axis=1)
y = np.asarray(yy, dtype=np.float64)
smiles_arr = np.asarray(structure_ids, dtype=object)
scaffold_arr = np.asarray(scaffolds, dtype=object)

metrics, predictions, splits = [], [], {}
for method, grouping in [('molecule', smiles_arr), ('scaffold', scaffold_arr)]:
    if len(set(grouping)) < 3:
        audit[f'{method}_split_skipped_too_few_groups'] += 1
        continue
    tr, te = next(GroupShuffleSplit(
        n_splits=1, test_size=0.20, random_state=42
    ).split(X, y, groups=grouping))
    assert not (set(grouping[tr]) & set(grouping[te])), 'Group leakage'
    assert not (set(smiles_arr[tr]) & set(smiles_arr[te])), 'Molecule leakage'
    splits[method] = {'train_pairs': len(tr), 'test_pairs': len(te),
                      'train_molecules': len(set(smiles_arr[tr])),
                      'test_molecules': len(set(smiles_arr[te])),
                      'train_scaffolds': len(set(scaffold_arr[tr])),
                      'test_scaffolds': len(set(scaffold_arr[te]))}
    for label, features in [('median_baseline', None),
                            ('structure', X), ('structure_solvent', XS)]:
        if features is None:
            yp = np.full(len(te), median(y[tr]), dtype=float)
        else:
            model = RandomForestRegressor(
                n_estimators=60, min_samples_leaf=3, max_depth=18,
                max_features=0.75, n_jobs=2, random_state=42
            )
            model.fit(features[tr], y[tr])
            yp = model.predict(features[te])
        mae = float(mean_absolute_error(y[te], yp))
        rmse = float(np.sqrt(mean_squared_error(y[te], yp)))
        r2 = float(r2_score(y[te], yp))
        metrics.append({'target': TARGET, 'split': method,
                        'model': label, 'n_train': len(tr), 'n_test': len(te),
                        'MAE_nm': round(mae, 4), 'RMSE_nm': round(rmse, 4),
                        'R2': round(r2, 4)})
        print(f'{TARGET} / {method} / {label}: MAE={mae:.2f} nm, RMSE={rmse:.2f} nm, R2={r2:.3f}', flush=True)
        for j, idx in enumerate(te):
            predictions.append({'target': TARGET, 'split': method,
                                'model': label, 'smiles': structure_ids[idx],
                                'solvent': solvent_names[idx],
                                'scaffold': scaffolds[idx],
                                'actual_nm': round(float(y[idx]), 5),
                                'predicted_nm': round(float(yp[j]), 5)})


def save_rows(path, rows):
    if not rows:
        return
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

save_rows(D / f'stage3_{TARGET}_metrics.csv', metrics)
save_rows(D / f'stage3_{TARGET}_predictions.csv', predictions)
report = {'target': TARGET, 'dataset': 'FluoDB GitHub CSV snapshot',
          'source': str(SOURCE), 'data_audit': dict(audit),
          'unique_molecules': len(set(structure_ids)),
          'unique_scaffolds': len(set(scaffolds)),
          'molecule_solvent_pairs': len(y),
          'solvents': len(sol_cache),
          'splits': splits, 'metrics': metrics,
          'feature_sets': {'structure': 'RDKit descriptors + ECFP4 (128 bits)',
                           'structure_solvent': 'structure + parsed solvent-SMILES descriptors'},
          'important_limitations': [
              'No DFT computed or imputed.',
              'Solvent SMILES descriptors are NOT dielectric or polarity measurements.',
              'Source DOI/protonation/experimental conditions may still confound labels.',
              'One grouped 80/20 split per method; uncertainty intervals not measured.',
              'Scaffold split tests extrapolation to unseen Murcko scaffolds.'
          ]}
(D / f'stage3_{TARGET}_audit.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
print('FINISHED:', D / f'stage3_{TARGET}_audit.json')
