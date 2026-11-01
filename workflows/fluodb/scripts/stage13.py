#!/usr/bin/env python3
"""XanPy stage 13: DOI-traceable, molecule/study-balanced solvent shifts
for *candidate* O-bridged conjugated xanthene-like structures.
Requires stage12_candidate_records.csv. Does not establish dye identity.
"""
import csv
import json
import math
import os
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
from statistics import median

import numpy as np
from rdkit import Chem, RDLogger

RDLogger.DisableLog('rdApp.error')
ROOT = Path(os.environ.get('XANPY_ROOT', str(Path.home() / 'Downloads/XanPy_QM_Atlas')))
SRC = ROOT / 'stage12_candidate_records.csv'
TARGET = 'O_bridge__conjugated_9C'
if not SRC.is_file():
    raise SystemExit(f'Missing: {SRC}')


def wave(text):
    try:
        x = float(text)
        return x if math.isfinite(x) and 200 <= x <= 1500 else None
    except (ValueError, TypeError):
        return None


def canonical_solvent(value):
    # SMILES is case-sensitive: never lower() solvent strings.
    if not value:
        return ''
    mol = Chem.MolFromSmiles(value.strip())
    if not mol or len(Chem.GetMolFrags(mol)) != 1:
        return ''
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def stokes(a, e):
    return 1e7 * (1.0 / a - 1.0 / e)


def write_csv(name, records, fields):
    with (ROOT / name).open('w', newline='', encoding='utf-8') as out:
        writer = csv.DictWriter(out, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


counts = Counter()
# (full SMILES, DOI) -> standardized solvent -> rows
studies = defaultdict(lambda: defaultdict(list))
original_tags = defaultdict(set)
original_sources = defaultdict(set)

with SRC.open(newline='', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        counts['all_stage12_candidate_rows'] += 1
        if row.get('chemical_state') != TARGET:
            continue
        counts['o_conjugated_candidate_rows'] += 1
        smi = (row.get('smiles') or '').strip()
        paper = (row.get('doi') or '').strip().lower()
        if not smi:
            counts['missing_structure'] += 1
            continue
        if not paper:
            counts['missing_doi'] += 1
            continue
        a = wave(row.get('absorption_nm'))
        e = wave(row.get('emission_nm'))
        if a is None or e is None:
            counts['missing_paired_spectrum'] += 1
            continue
        if e < a:
            counts['negative_stokes'] += 1
            continue
        solvent = canonical_solvent(row.get('solvent') or '')
        if not solvent:
            counts['invalid_solvent'] += 1
            continue
        key = (smi, paper)
        studies[key][solvent].append((a, e, row.get('source_row') or ''))
        original_tags[key].add((row.get('source_tag') or '').strip())
        original_sources[key].add((row.get('source_database') or '').strip())
        counts['valid_input_measurements'] += 1

pairs = []
study_rows = []
review_rows = []

def flatten_rows(group):
    return sum(len(x) for x in group.values())

for (smi, doi), solvents in sorted(studies.items()):
    if len(solvents) < 2:
        continue
    spectral = {}
    duplicate_ranges = []
    for solvent, measurements in solvents.items():
        a = median(v[0] for v in measurements)
        e = median(v[1] for v in measurements)
        spectral[solvent] = {'abs': a, 'em': e, 'stokes': stokes(a, e),
                             'n': len(measurements)}
        if len(measurements) > 1:
            duplicate_ranges.append(max(v[1] for v in measurements) -
                                    min(v[1] for v in measurements))
    group_pairs = []
    for solvent_a, solvent_b in combinations(sorted(spectral), 2):
        a = spectral[solvent_a]
        b = spectral[solvent_b]
        ea, eb = b['em'] - a['em'], b['abs'] - a['abs']
        em_cm = 1e7 / b['em'] - 1e7 / a['em']
        abs_cm = 1e7 / b['abs'] - 1e7 / a['abs']
        group_pairs.append({
            'smiles': smi, 'doi': doi,
            'solvent_A': solvent_a, 'solvent_B': solvent_b,
            'n_measurements_A': a['n'], 'n_measurements_B': b['n'],
            'delta_em_nm': round(ea, 5), 'delta_abs_nm': round(eb, 5),
            'delta_stokes_cm1': round(b['stokes'] - a['stokes'], 5),
            'em_dominant_nm': int(abs(ea) > abs(eb)),
            'em_dominant_cm1': int(abs(em_cm) > abs(abs_cm)),
            'abs_em_shift_nm': round(abs(ea), 5),
            'abs_abs_shift_nm': round(abs(eb), 5),
            'abs_em_shift_cm1': round(abs(em_cm), 5),
            'abs_abs_shift_cm1': round(abs(abs_cm), 5)
        })
    pairs.extend(group_pairs)
    group_stat = {
        'smiles': smi,
        'doi': doi,
        'n_solvents': len(spectral),
        'n_pairs': len(group_pairs),
        'n_measurements': flatten_rows(solvents),
        'emission_larger_fraction_nm': round(float(np.mean([p['em_dominant_nm'] for p in group_pairs])), 6),
        'emission_larger_fraction_cm1': round(float(np.mean([p['em_dominant_cm1'] for p in group_pairs])), 6),
        'mean_abs_em_shift_nm': round(float(np.mean([p['abs_em_shift_nm'] for p in group_pairs])), 5),
        'mean_abs_abs_shift_nm': round(float(np.mean([p['abs_abs_shift_nm'] for p in group_pairs])), 5),
        'mean_abs_em_shift_cm1': round(float(np.mean([p['abs_em_shift_cm1'] for p in group_pairs])), 5),
        'mean_abs_abs_shift_cm1': round(float(np.mean([p['abs_abs_shift_cm1'] for p in group_pairs])), 5),
        'mean_difference_abs_shifts_nm': round(float(np.mean([
            p['abs_em_shift_nm'] - p['abs_abs_shift_nm'] for p in group_pairs
        ])), 5),
        'source_tags': '|'.join(sorted(original_tags[(smi, doi)])),
        'source_databases': '|'.join(sorted(original_sources[(smi, doi)])),
        'solvent_smiles': '|'.join(sorted(spectral)),
        'max_repeated_solvent_em_disagreement_nm': round(max(duplicate_ranges), 3)
            if duplicate_ranges else 0.0
    }
    study_rows.append(group_stat)
    review_rows.append(group_stat)

study_rows.sort(key=lambda r: (-r['n_solvents'], -r['mean_abs_em_shift_nm']))
review_rows.sort(key=lambda r: (-r['max_repeated_solvent_em_disagreement_nm'], -r['mean_abs_em_shift_nm']))

pair_fields = ['smiles', 'doi', 'solvent_A', 'solvent_B',
               'n_measurements_A', 'n_measurements_B',
               'delta_em_nm', 'delta_abs_nm', 'delta_stokes_cm1',
               'em_dominant_nm', 'em_dominant_cm1',
               'abs_em_shift_nm', 'abs_abs_shift_nm',
               'abs_em_shift_cm1', 'abs_abs_shift_cm1']
study_fields = list(study_rows[0]) if study_rows else [
    'smiles', 'doi', 'n_solvents', 'n_pairs', 'n_measurements',
    'emission_larger_fraction_nm', 'emission_larger_fraction_cm1',
    'mean_abs_em_shift_nm', 'mean_abs_abs_shift_nm',
    'mean_abs_em_shift_cm1', 'mean_abs_abs_shift_cm1',
    'mean_difference_abs_shifts_nm', 'source_tags', 'source_databases',
    'solvent_smiles', 'max_repeated_solvent_em_disagreement_nm']
write_csv('stage13_o_conjugated_solvent_pairs.csv', pairs, pair_fields)
write_csv('stage13_o_conjugated_studies.csv', study_rows, study_fields)
write_csv('stage13_manual_review_priority.csv', review_rows, study_fields)


def boot_cluster(key, group_key, n=2000):
    # Average the group statistics *within* each clustering unit first.
    # Molecule-level or DOI-level bootstrap is a sensitivity analysis,
    # not two independent sets of observations.
    group = defaultdict(list)
    for r in study_rows:
        group[r[group_key]].append(float(r[key]))
    x = np.array([np.mean(g) for g in group.values()], dtype=float)
    if len(x) < 2:
        return None
    rng = np.random.default_rng(20261101)
    estimates = np.empty(n)
    for b in range(n):
        estimates[b] = x[rng.integers(0, len(x), size=len(x))].mean()
    return {'estimate': round(float(x.mean()), 5),
            'ci95': [round(float(z), 5) for z in np.quantile(estimates, [0.025, 0.975])],
            'n_clusters': int(len(x))}

stat = {}
for measure in ['emission_larger_fraction_nm', 'emission_larger_fraction_cm1',
                'mean_difference_abs_shifts_nm']:
    stat[measure] = {
        'molecule_cluster': boot_cluster(measure, 'smiles'),
        'doi_cluster': boot_cluster(measure, 'doi')
    }

audit = {
    'source_file': 'stage12_candidate_records.csv',
    'class_exact': TARGET,
    'counts': dict(counts),
    'doi_traceable_molecule_paper_groups_with_at_least_2_solvents': len(study_rows),
    'unique_molecules_with_multisolvent_study': len({r['smiles'] for r in study_rows}),
    'unique_dois_with_multisolvent_study': len({r['doi'] for r in study_rows}),
    'within_study_solvent_pairs_not_independent': len(pairs),
    'studies_with_repeated_solvent_disagreement_over_20_nm': sum(
        r['max_repeated_solvent_em_disagreement_nm'] > 20 for r in study_rows
    ),
    'cluster_bootstrap': stat,
    'limitations': [
        'O-conjugated classification is topological screening, not verified xanthene identity.',
        'No pH, temperature, concentration or protomer control in the source.',
        'Solvent pairs from one molecule/DOI are statistically dependent.',
        'Molecule and DOI cluster bootstrap are sensitivity checks, not external validation.',
        'Do not generalize this subset to Si/Ge: those classes have no DOI-traceable multi-solvent groups.',
        'No new DFT, no HOMO/LUMO values and no mechanistic proof.'
    ]
}
(ROOT / 'stage13_o_conjugated_audit.json').write_text(json.dumps(audit, indent=2), encoding='utf-8')
print(json.dumps(audit, indent=2))
print('FINISHED')
