import csv, json, os
from collections import defaultdict
from pathlib import Path
from statistics import median
import numpy as np

ROOT = Path(os.environ.get('XANPY_ROOT', str(Path.home() / 'Downloads' / 'XanPy_QM_Atlas')))
SRC = ROOT / 'stage6_within_paper_solvent_pairs.csv'
if not SRC.is_file():
    raise SystemExit(f'Missing input: {SRC}')

by_group = defaultdict(list)
rows = 0
with SRC.open(newline='', encoding='utf-8') as fh:
    reader = csv.DictReader(fh)
    needed = {'canonical_smiles', 'doi', 'delta_abs_nm', 'delta_em_nm',
              'delta_abs_wavenumber_cm1', 'delta_em_wavenumber_cm1'}
    if not needed.issubset(reader.fieldnames or []):
        raise SystemExit(f'Missing columns: {needed - set(reader.fieldnames or [])}')
    for r in reader:
        try:
            a = abs(float(r['delta_abs_nm']))
            e = abs(float(r['delta_em_nm']))
            aw = abs(float(r['delta_abs_wavenumber_cm1']))
            ew = abs(float(r['delta_em_wavenumber_cm1']))
        except (TypeError, ValueError):
            continue
        if not all(np.isfinite(v) for v in (a,e,aw,ew)):
            continue
        by_group[(r['canonical_smiles'], r['doi'])].append((a, e, aw, ew))
        rows += 1

if not by_group:
    raise SystemExit('No valid solvent-pair records')

groups = []
for (smi, doi), entries in by_group.items():
    x = np.asarray(entries, dtype=float)
    groups.append({
        'canonical_smiles': smi,
        'doi': doi,
        'n_pairs': len(entries),
        'fraction_emission_larger_nm': float(np.mean(x[:,1] > x[:,0])),
        'fraction_emission_larger_cm1': float(np.mean(x[:,3] > x[:,2])),
        'mean_abs_delta_em_nm': float(x[:,1].mean()),
        'mean_abs_delta_abs_nm': float(x[:,0].mean()),
        'mean_difference_nm': float(np.mean(x[:,1] - x[:,0])),
        'mean_difference_cm1': float(np.mean(x[:,3] - x[:,2])),
        'median_abs_delta_em_nm': float(np.median(x[:,1])),
        'median_abs_delta_abs_nm': float(np.median(x[:,0]))
    })

by_mol = defaultdict(list)
by_doi = defaultdict(list)
for g in groups:
    by_mol[g['canonical_smiles']].append(g)
    by_doi[g['doi']].append(g)

keys = ['fraction_emission_larger_nm', 'fraction_emission_larger_cm1',
        'mean_difference_nm', 'mean_difference_cm1']

def summarize_units(d):
    return np.array([[np.mean([g[k] for g in val]) for k in keys]
                     for val in d.values()], dtype=float)

rng = np.random.default_rng(20261101)
def boot(arr, reps=2000):
    n = len(arr)
    out = np.empty((reps, arr.shape[1]))
    for i in range(reps):
        ix = rng.integers(0, n, n)
        out[i] = arr[ix].mean(axis=0)
    return {
      k: {
        'estimate': round(float(arr[:,j].mean()), 5),
        'ci95': [round(float(t), 5) for t in np.quantile(out[:,j], [0.025, 0.975])]
      } for j,k in enumerate(keys)
    }

m_balanced = boot(summarize_units(by_mol))
doi_balanced = boot(summarize_units(by_doi))

total_pairs = sum(g['n_pairs'] for g in groups)
share_top10 = sum(g['n_pairs'] for g in sorted(groups, key=lambda g:g['n_pairs'], reverse=True)[:10]) / total_pairs
trimmed = [g for g in groups if g['n_pairs'] <= 20]

report = {
    'valid_pairs': rows,
    'molecule_doi_groups': len(groups),
    'unique_molecules': len(by_mol),
    'unique_dois': len(by_doi),
    'pairs_from_top10_molecule_doi_groups_fraction': round(share_top10, 5),
    'groups_over20_pairs': sum(g['n_pairs'] > 20 for g in groups),
    'pair_weighted_fraction_emission_larger_nm': round(sum(
        g['fraction_emission_larger_nm']*g['n_pairs'] for g in groups) / total_pairs, 5),
    'group_balanced_fraction_emission_larger_nm': round(float(np.mean(
        [g['fraction_emission_larger_nm'] for g in groups])), 5),
    'trimmed_groups_le20_pairs_fraction_emission_larger_nm': round(float(np.mean(
        [g['fraction_emission_larger_nm'] for g in trimmed])), 5) if trimmed else None,
    'molecule_balanced_bootstrap': m_balanced,
    'doi_balanced_bootstrap': doi_balanced,
    'cautions': [
        'Pairwise combinations from one molecule/paper are dependent.',
        'Molecule and DOI bootstrap are alternative sensitivity analyses, not independent replications.',
        'Same DOI does not guarantee identical pH, protonation or experimental conditions.',
        'Solvent shifts are observational associations, not causal proof.',
        'Electronic excited-state relaxation cannot be established from this comparison alone.'
    ]
}

with (ROOT / 'stage7_group_balanced.csv').open('w', newline='', encoding='utf-8') as fh:
    writer = csv.DictWriter(fh, fieldnames=list(groups[0]))
    writer.writeheader()
    writer.writerows(sorted(groups, key=lambda g:g['n_pairs'], reverse=True))
(ROOT / 'stage7_robustness_audit.json').write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding='utf-8')
print(json.dumps(report, indent=2, ensure_ascii=False))
print('FINISHED')
