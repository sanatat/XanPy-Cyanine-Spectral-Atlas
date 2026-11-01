import csv, json, math, re, os
from collections import defaultdict
from pathlib import Path
from statistics import median
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors, Crippen

RDLogger.DisableLog('rdApp.error')
ROOT = Path(os.environ.get('XANPY_ROOT', str(Path.home() / 'Downloads' / 'XanPy_QM_Atlas')))
SRC = ROOT / 'stage2_all_endpoints.csv'
if not SRC.exists():
    raise SystemExit(f'Missing: {SRC}')

def wave(s):
    try:
        x = float(s)
        return x if math.isfinite(x) and 200 <= x <= 1500 else None
    except (ValueError, TypeError):
        return None

def doi_key(s):
    m = re.search(r'10\.\d{4,9}/[^\s,;]+', s or '', re.I)
    return m.group(0).rstrip(').,;').lower() if m else ''

solvent_cache = {}
def canonical_solvent(s):
    if s not in solvent_cache:
        m = Chem.MolFromSmiles(s) if s else None
        solvent_cache[s] = (
            Chem.MolToSmiles(m, isomericSmiles=True)
            if m and len(Chem.GetMolFrags(m)) == 1 else None
        )
    return solvent_cache[s]

# Targeted, fixed solvent contrasts. Second minus first.
raw_contrasts = [
    ('methanol', 'CO', 'acetonitrile', 'CC#N'),
    ('dichloromethane', 'ClCCl', 'methanol', 'CO'),
    ('chloroform', 'ClC(Cl)Cl', 'methanol', 'CO'),
    ('ethanol', 'CCO', 'methanol', 'CO'),
    ('dimethyl_sulfoxide', 'CS(C)=O', 'methanol', 'CO'),
    ('water', 'O', 'methanol', 'CO'),
    ('toluene', 'Cc1ccccc1', 'acetonitrile', 'CC#N'),
]
contrasts = [
    (n1, canonical_solvent(s1), n2, canonical_solvent(s2))
    for n1,s1,n2,s2 in raw_contrasts
]
if any(not s1 or not s2 for _,s1,_,s2 in contrasts):
    raise SystemExit('A reference solvent SMILES could not be parsed')

measurements = defaultdict(list)
counters = defaultdict(int)
with SRC.open(newline='', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        counters['read_rows'] += 1
        a = wave(row.get('absorption_nm'))
        e = wave(row.get('emission_nm'))
        if a is None or e is None or e < a:
            counters['incomplete_or_negative_stokes'] += 1
            continue
        doi = doi_key(row.get('reference_doi'))
        mol = (row.get('canonical_smiles') or '').strip()
        sol = canonical_solvent((row.get('solvent') or '').strip())
        if not (doi and mol and sol):
            counters['invalid_study_or_solvent'] += 1
            continue
        measurements[(mol,doi,sol)].append((a,e))
        counters['eligible_rows'] += 1

by_study = defaultdict(dict)
for (mol,doi,sol),vals in measurements.items():
    # Median each endpoint, do not count repeated rows as independent studies.
    a = median(v[0] for v in vals)
    e = median(v[1] for v in vals)
    stokes = 1e7 * (1/a - 1/e)
    by_study[(mol,doi)][sol] = (a,e,stokes,len(vals))

rows = []
for (mol,doi), spectra in by_study.items():
    for a_name,a_sm,b_name,b_sm in contrasts:
        if a_sm not in spectra or b_sm not in spectra:
            continue
        a = spectra[a_sm]
        b = spectra[b_sm]
        rows.append({
            'canonical_smiles':mol, 'doi':doi,
            'solvent_A':a_name, 'solvent_B':b_name,
            'n_records_A':a[3], 'n_records_B':b[3],
            'delta_abs_nm':round(b[0]-a[0],4),
            'delta_em_nm':round(b[1]-a[1],4),
            'delta_stokes_cm1':round(b[2]-a[2],3),
            'abs_em_greater':int(abs(b[1]-a[1])>abs(b[0]-a[0]))
        })

columns = ['canonical_smiles','doi','solvent_A','solvent_B','n_records_A',
           'n_records_B','delta_abs_nm','delta_em_nm','delta_stokes_cm1',
           'abs_em_greater']
with (ROOT/'stage8_selected_solvent_pairs.csv').open('w',newline='') as f:
    wr=csv.DictWriter(f,fieldnames=columns)
    wr.writeheader(); wr.writerows(rows)

rng = np.random.default_rng(20261101)
def equal_molecule_boot(vals, repetitions=1000):
    # vals = records (molecule_id, signed_delta); cluster by molecule.
    per_molecule = defaultdict(list)
    for mol, x in vals:
        per_molecule[mol].append(x)
    z = np.asarray([np.mean(v) for v in per_molecule.values()],dtype=float)
    if len(z) == 0: return None
    means = np.asarray([z[rng.integers(0,len(z),len(z))].mean()
                        for _ in range(repetitions)])
    return {'estimate':round(float(z.mean()),4),
            'ci95':[round(float(x),4) for x in np.quantile(means,[.025,.975])],
            'unique_molecules':len(z)}

report=[]
for a_name,a_sm,b_name,b_sm in contrasts:
    subset=[r for r in rows if r['solvent_A']==a_name and r['solvent_B']==b_name]
    if not subset: continue
    stats={'solvent_A':a_name,'solvent_B':b_name,
           'n_molecule_doi_groups':len(subset),
           'n_unique_molecules':len(set(r['canonical_smiles'] for r in subset))}
    for feature in ('delta_em_nm','delta_abs_nm','delta_stokes_cm1'):
        stats[feature] = equal_molecule_boot(
            [(r['canonical_smiles'],r[feature]) for r in subset]
        )
    stats['fraction_emission_change_larger'] = equal_molecule_boot(
        [(r['canonical_smiles'],r['abs_em_greater']) for r in subset]
    )
    report.append(stats)

summary={'input_file':'stage2_all_endpoints.csv',
         'counts':dict(counters),
         'study_solvent_groups':len(measurements),
         'selected_contrast_groups':len(rows),
         'selected_pairs':report,
         'interpretation':[
           'Delta is solvent B minus solvent A.',
           'Every molecule-DOI-solvent is collapsed before comparison.',
           'Bootstrap equalizes molecule weights; DOI dependence remains.',
           'Same DOI does not control temperature, pH, dye speciation, or instrumentation.',
           'This is descriptive association, not causal or DFT evidence.'
         ]}
(ROOT/'stage8_selected_solvent_audit.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
print('FINISHED')
