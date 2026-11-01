#!/usr/bin/env python3
"""All solvent pairs within structure and DOI; not independent observations."""
import os,csv,json,re,math
from pathlib import Path
from collections import Counter,defaultdict
from statistics import median
from itertools import combinations
from rdkit import Chem,RDLogger
RDLogger.DisableLog('rdApp.error')
ROOT=Path(os.environ.get('XANPY_ROOT',str(Path.home()/'Downloads/XanPy_QM_Atlas')))
SRC=ROOT/'stage2_all_endpoints.csv'
if not SRC.exists():raise SystemExit('Missing '+str(SRC))
def wave(x):
    try:
        v=float(x);return v if math.isfinite(v) and 200<=v<=1500 else None
    except (ValueError,TypeError):return None
def doi(s):
    m=re.search(r'10\.\d{4,9}/[^\s,;]+',s or '',re.I)
    return m.group().rstrip(').,;').lower() if m else ''
cache={}
def sol(s):
    if s not in cache:
        m=Chem.MolFromSmiles(s) if s else None
        cache[s]=Chem.MolToSmiles(m,isomericSmiles=True) if m and len(Chem.GetMolFrags(m))==1 else ''
    return cache[s]
counts=Counter();studies=defaultdict(lambda:defaultdict(list))
with SRC.open(newline='',encoding='utf-8') as f:
    for r in csv.DictReader(f):
        counts['input_rows']+=1
        a=wave(r.get('absorption_nm'));e=wave(r.get('emission_nm'))
        if a is None or e is None:counts['missing_paired_spectrum']+=1;continue
        if e<a:counts['negative_stokes_excluded']+=1;continue
        paper=doi(r.get('reference_doi'))
        if not paper:counts['missing_doi']+=1;continue
        smi=(r.get('canonical_smiles') or '').strip();solvent=sol((r.get('solvent') or '').strip())
        if not smi or not solvent:counts['invalid_solvent_or_structure']+=1;continue
        studies[(smi,paper)][solvent].append((a,e));counts['eligible_rows']+=1
rows=[]
for (smi,paper),solvents in studies.items():
    if len(solvents)<2:continue
    v={s:(median(a for a,e in records),median(e for a,e in records),median(1e7*(1/a-1/e) for a,e in records),len(records)) for s,records in solvents.items()}
    for s1,s2 in combinations(sorted(v),2):
        a,b=v[s1],v[s2]
        rows.append({'canonical_smiles':smi,'doi':paper,'solvent_A_smiles':s1,'solvent_B_smiles':s2,
          'n_measurements_A':a[3],'n_measurements_B':b[3],
          'abs_A_nm':round(a[0],4),'abs_B_nm':round(b[0],4),
          'em_A_nm':round(a[1],4),'em_B_nm':round(b[1],4),
          'delta_abs_nm':round(b[0]-a[0],4),'delta_em_nm':round(b[1]-a[1],4),
          'delta_abs_wavenumber_cm1':round(1e7/b[0]-1e7/a[0],2),
          'delta_em_wavenumber_cm1':round(1e7/b[1]-1e7/a[1],2),
          'delta_stokes_cm1':round(b[2]-a[2],2),
          'abs_delta_stokes_cm1':round(abs(b[2]-a[2]),2)})
fields=['canonical_smiles','doi','solvent_A_smiles','solvent_B_smiles','n_measurements_A','n_measurements_B','abs_A_nm','abs_B_nm','em_A_nm','em_B_nm','delta_abs_nm','delta_em_nm','delta_abs_wavenumber_cm1','delta_em_wavenumber_cm1','delta_stokes_cm1','abs_delta_stokes_cm1']
with (ROOT/'stage6_within_paper_solvent_pairs.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
report={'quality_control':dict(counts),'same_paper_solvent_pairs':len(rows),'unique_molecules':len({r['canonical_smiles'] for r in rows}),
'unique_molecule_doi_groups':len({(r['canonical_smiles'],r['doi']) for r in rows}),
'median_absolute_emission_shift_nm':median(abs(r['delta_em_nm']) for r in rows) if rows else None,
'median_absolute_absorption_shift_nm':median(abs(r['delta_abs_nm']) for r in rows) if rows else None,
'median_absolute_stokes_change_cm1':median(r['abs_delta_stokes_cm1'] for r in rows) if rows else None,
'fraction_emission_shift_larger_than_absorption':round(sum(abs(r['delta_em_nm'])>abs(r['delta_abs_nm']) for r in rows)/len(rows),4) if rows else None,
'note':'Descriptive within-DOI comparisons; pairs are dependent and chemical state not controlled'}
(ROOT/'stage6_audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
