import csv, json, re, sys, os
from pathlib import Path
from collections import Counter, defaultdict
import numpy as np

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(os.environ.get('XANPY_ROOT', str(Path.home() / 'Downloads/XanPy_QM_Atlas')))
SRC = ROOT / 'stage2_all_endpoints.csv'
PAIRS = ROOT / 'stage8_selected_solvent_pairs.csv'
if not SRC.exists() or not PAIRS.exists():
    raise SystemExit('Missing stage2_all_endpoints.csv or stage8_selected_solvent_pairs.csv')

DOI_PATTERN = re.compile(r'10\.\d{4,9}/[^\s,;]+', re.I)
def normalize_doi(text):
    m = DOI_PATTERN.search(text or '')
    return m.group(0).rstrip(').,;').lower() if m else ''

def number(value):
    try:
        v = float(value)
        return v if np.isfinite(v) else None
    except (TypeError, ValueError):
        return None

# Preserve original source tags; tags are NOT necessarily verified chemical families.
labels = defaultdict(Counter)
with SRC.open(newline='', encoding='utf-8') as fh:
    for row in csv.DictReader(fh):
        key = ((row.get('canonical_smiles') or '').strip(),
               normalize_doi(row.get('reference_doi')))
        if not all(key):
            continue
        tag = (row.get('source_tag_name') or '').strip() or 'unclassified'
        labels[key][tag] += 1

contrasts = defaultdict(list)
audit = Counter()
with PAIRS.open(newline='', encoding='utf-8') as fh:
    for row in csv.DictReader(fh):
        audit['input_pairs'] += 1
        smi = (row.get('canonical_smiles') or '').strip()
        doi = normalize_doi(row.get('doi'))
        if not smi or not doi:
            audit['missing_molecule_or_doi'] += 1
            continue
        tagcounts = labels.get((smi, doi), Counter())
        if not tagcounts:
            family = 'unmatched_source'
            audit['unmatched_source'] += 1
        elif len(tagcounts) == 1:
            family = next(iter(tagcounts))
            audit['unambiguous_source'] += 1
        else:
            family = 'ambiguous_source_tag'
            audit['ambiguous_source_tag'] += 1
        em = number(row.get('delta_em_nm'))
        ab = number(row.get('delta_abs_nm'))
        st = number(row.get('delta_stokes_cm1'))
        if None in (em, ab, st):
            audit['invalid_shift'] += 1
            continue
        contrast = (row.get('solvent_A') or '') + ' -> ' + (row.get('solvent_B') or '')
        contrasts[(contrast, family)].append({
            'smiles':smi,'doi':doi,'em':em,'abs':ab,'stokes':st,
            'dominant':int(abs(em)>abs(ab))
        })
        audit['valid_pairs'] += 1

rng = np.random.default_rng(26011)

def summarize(rows):
    # Equal-weight independent molecules; collapse multiple DOIs for a molecule.
    molecules = defaultdict(list)
    papers = defaultdict(list)
    for r in rows:
        molecules[r['smiles']].append(r)
        papers[r['doi']].append(r)
    means = lambda grp: np.array([
        [np.mean([r[k] for r in v]) for k in ('em','abs','stokes','dominant')]
        for v in grp.values()], dtype=float)
    mol = means(molecules)
    doi = means(papers)
    def boot(arr, col, reps=1000):
        v = arr[:,col]
        if len(v) < 2: return [None, None]
        out = np.array([np.mean(v[rng.integers(0,len(v),len(v))])
                        for _ in range(reps)])
        return [round(float(z), 3) for z in np.quantile(out,[0.025,0.975])]
    return {
        'pairs':len(rows), 'molecules':len(mol), 'dois':len(doi),
        'delta_em_mean_nm':round(float(mol[:,0].mean()), 3),
        'delta_em_CI95_molecule':boot(mol,0),
        'delta_abs_mean_nm':round(float(mol[:,1].mean()),3),
        'delta_stokes_mean_cm1':round(float(mol[:,2].mean()),3),
        'em_dominant_fraction':round(float(mol[:,3].mean()),4),
        'em_dominant_CI95_molecule':boot(mol,3),
        'delta_em_mean_nm_DOI_balanced':round(float(doi[:,0].mean()),3),
        'delta_em_CI95_DOI':boot(doi,0),
    }

stats = []
for (contrast, family), rows in sorted(contrasts.items()):
    if len({r['smiles'] for r in rows}) < 25:
        continue
    stats.append({'contrast':contrast,'source_family_label':family,**summarize(rows)})

out = ROOT / 'stage9_family_solvent_comparison.csv'
fieldnames = ['contrast','source_family_label','pairs','molecules','dois',
              'delta_em_mean_nm','delta_em_CI95_molecule','delta_abs_mean_nm',
              'delta_stokes_mean_cm1','em_dominant_fraction',
              'em_dominant_CI95_molecule','delta_em_mean_nm_DOI_balanced',
              'delta_em_CI95_DOI']
with out.open('w',newline='', encoding='utf-8') as fh:
    w = csv.DictWriter(fh,fieldnames=fieldnames)
    w.writeheader()
    for stat in stats:
        copy = stat.copy()
        for k,v in copy.items():
            if isinstance(v,list):copy[k] = json.dumps(v)
        w.writerow(copy)

report = {
    'audit':dict(audit),
    'source_molecule_doi_keys':len(labels),
    'all_contrast_family_groups':len(contrasts),
    'reported_strata_minimum_25_molecules':len(stats),
    'top_by_sample_size':sorted(stats,key=lambda r:-r['molecules'])[:18],
    'important_notes':[
        'source_family_label is inherited from FluoDB and not a verified xanthene/pyronin substructure label',
        'Molecule and DOI intervals are separate resampling sensitivity analyses',
        'Signed shifts depend on A-to-B direction; CI excludes zero is exploratory only',
        'Some family labels are source taxonomies rather than chemical scaffold names',
        'This does not prove a causal solvent mechanism'
    ]
}
(ROOT/'stage9_family_audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print('=== STAGE9 ===')
print('Matched pairs:', audit['valid_pairs'])
print('Unmatched:', audit['unmatched_source'], 'Ambiguous:',audit['ambiguous_source_tag'])
print('Family/contrast groups with >=25 molecules:',len(stats))
for s in sorted(stats,key=lambda r:-r['molecules'])[:12]:
    print(f"{s['contrast']:36} | {s['source_family_label']:16} | N={s['molecules']:4d} | emission={s['delta_em_mean_nm']:+7.2f} nm | DOI-balanced={s['delta_em_mean_nm_DOI_balanced']:+7.2f} nm")
print('Saved:',out)
print('FINISHED')
