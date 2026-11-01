import csv, json, math, os, re
from pathlib import Path
from collections import Counter, defaultdict
from itertools import combinations
from statistics import median, mean
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold
RDLogger.DisableLog('rdApp.error')
R = Path(os.environ.get('XANPY_ROOT', str(Path.home()/'Downloads/XanPy_QM_Atlas')))
source = R/'stage2_all_endpoints.csv'
if not source.is_file(): raise SystemExit('Missing '+str(source))

def doi(s):
    m = re.search(r'10\.\d{4,9}/[^\s,;]+', s or '', re.I)
    return m.group().rstrip(').,;').lower() if m else ''

def num(x):
    try:
        v=float(x)
        return v if math.isfinite(v) else None
    except (TypeError,ValueError): return None

def motif(m):
    rings = [set(x) for x in m.GetRingInfo().AtomRings() if len(x)==6]
    found=set()
    for center in rings:
        flanks = [ring for ring in rings if ring != center
                  and len(center & ring)==2 and all(
                      m.GetAtomWithIdx(i).GetAtomicNum()==6 and
                      m.GetAtomWithIdx(i).GetIsAromatic() for i in ring)]
        for left,right in combinations(flanks,2):
            if (left & right)-center: continue
            bridge=center-left-right
            if len(bridge)!=2: continue
            symbols=sorted(m.GetAtomWithIdx(i).GetSymbol() for i in bridge)
            if symbols==['C','O']: found.add('O_xanthene_like')
            elif symbols==['C','Si']: found.add('Si_xanthene_like')
            elif symbols==['C','Ge']: found.add('Ge_xanthene_like')
            elif symbols==['C','C']:
                if any(not m.GetAtomWithIdx(i).GetIsAromatic() for i in bridge):
                    found.add('C_bridge_candidate')
                else: found.add('all_C_aromatic_tricycle')
    for name in ('O_xanthene_like','Si_xanthene_like','Ge_xanthene_like',
                 'C_bridge_candidate','all_C_aromatic_tricycle'):
        if name in found:return name
    return 'other'

# Structural unit tests (reference dyes and negative controls).
checks={
 '9H_xanthene':('C1C2=CC=CC=C2OC3=CC=CC=C31','O_xanthene_like'),
 'Rhodamine_B_cation':('CCN(CC)C1=CC2=C(C=C1)C(=C3C=CC(=[N+](CC)CC)C=C3O2)C4=CC=CC=C4C(=O)O','O_xanthene_like'),
 'benzene':('c1ccccc1','other'),
 'Si_analog':('c1ccc2c(c1)[SiH2]c1ccccc1C2','Si_xanthene_like'),
 'Ge_analog':('c1ccc2c(c1)[GeH2]c1ccccc1C2','Ge_xanthene_like')}
for label,(smiles,expected) in checks.items():
    mol=Chem.MolFromSmiles(smiles)
    if mol is None or motif(mol)!=expected:
        raise SystemExit('Structural self-test failed: '+label)

mols={}; solvents={}; n=Counter(); by_motif=defaultdict(lambda: {
    'molecules':set(), 'dois':set(), 'rows':0, 'source_tags':Counter()})
records=[]; spectra=defaultdict(list)
with source.open(newline='',encoding='utf-8') as f:
    for row in csv.DictReader(f):
        n['input_rows']+=1
        smi=(row.get('canonical_smiles') or '').strip()
        if smi not in mols:
            m=Chem.MolFromSmiles(smi) if smi else None
            mols[smi]=(motif(m),Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m))) if m else ('invalid','')
        kind,scaffold=mols[smi]
        state=by_motif[kind]
        state['molecules'].add(smi);state['rows']+=1
        tag=(row.get('source_tag_name') or 'unclassified').strip()
        state['source_tags'][tag]+=1
        d=doi(row.get('reference_doi'))
        if d:state['dois'].add(d)
        if kind not in ('O_xanthene_like','Si_xanthene_like','Ge_xanthene_like',
                        'C_bridge_candidate','all_C_aromatic_tricycle'):
            continue
        n['candidate_rows']+=1
        raw=(row.get('solvent') or '').strip()
        if raw not in solvents:
            sm=Chem.MolFromSmiles(raw) if raw else None
            solvents[raw]=Chem.MolToSmiles(sm,isomericSmiles=True) if sm and len(Chem.GetMolFrags(sm))==1 else ''
        sol=solvents[raw]
        records.append({
          'source_row':row.get('source_row',''), 'kind':kind,'smiles':smi,
          'scaffold':scaffold,'tag':tag,'doi':d,'solvent':sol,
          'abs_nm':row.get('absorption_nm',''),
          'em_nm':row.get('emission_nm',''),
          'plqy':row.get('plqy_fraction','')})
        a,e=num(row.get('absorption_nm')),num(row.get('emission_nm'))
        if a is not None and e is not None and 200<=a<=1500 and 200<=e<=1500 and e>=a and d and sol:
            spectra[(kind,smi,d,sol)].append((a,e))

def save(name,rows,columns):
    with (R/name).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows(rows)
rec_cols=['source_row','kind','smiles','scaffold','tag','doi','solvent','abs_nm','em_nm','plqy']
save('stage11_all_tricycle_candidates.csv',records,rec_cols)
save('stage11_hetero_xanthene_candidates.csv',[
    x for x in records if x['kind'] in ('O_xanthene_like','Si_xanthene_like','Ge_xanthene_like')],rec_cols)

by_study=defaultdict(dict)
for (kind,smi,d,sol),values in spectra.items():
    by_study[(kind,smi,d)][sol]=(median(a for a,e in values),median(e for a,e in values))
pairs=[]; group_stats=defaultdict(list)
for (kind,smi,d),sols in by_study.items():
    if len(sols)<2:continue
    indicators=[]
    for (sa,(aa,ea)),(sb,(ab,eb)) in combinations(sorted(sols.items()),2):
        da=ab-aa;de=eb-ea
        wa=1e7/ab-1e7/aa;we=1e7/eb-1e7/ea
        st=(1e7/ab-1e7/eb)-(1e7/aa-1e7/ea)
        indicators.append((int(abs(de)>abs(da)),int(abs(we)>abs(wa))))
        pairs.append({'kind':kind,'smiles':smi,'doi':d,
          'solvent_A':sa,'solvent_B':sb,'delta_abs_nm':round(da,4),
          'delta_em_nm':round(de,4),'delta_stokes_cm1':round(st,3)})
    group_stats[kind].append({
       'smiles':smi,'doi':d,'n_solvents':len(sols),'n_pairs':len(indicators),
       'fraction_em_gt_abs_nm':mean(x[0] for x in indicators),
       'fraction_em_gt_abs_cm1':mean(x[1] for x in indicators)})
save('stage11_all_candidate_solvent_pairs.csv',pairs,
     ['kind','smiles','doi','solvent_A','solvent_B','delta_abs_nm','delta_em_nm','delta_stokes_cm1'])
stat_rows=[{'kind':k,**rec} for k,rs in group_stats.items() for rec in rs]
save('stage11_molecule_doi_summary.csv',stat_rows,
     ['kind','smiles','doi','n_solvents','n_pairs','fraction_em_gt_abs_nm','fraction_em_gt_abs_cm1'])

summary={}
for kind,state in sorted(by_motif.items()):
    groups=group_stats[kind]
    summary[kind]={
       'records':state['rows'], 'unique_molecules':len(state['molecules']),
       'distinct_dois':len(state['dois']),
       'multi_solvent_molecule_doi_groups':len(groups),
       'multi_solvent_molecules':len({g['smiles'] for g in groups}),
       'multi_solvent_pairs':sum(g['n_pairs'] for g in groups),
       'study_equal_weight_fraction_emission_dominant_nm':round(mean(
            g['fraction_em_gt_abs_nm'] for g in groups),4) if groups else None,
       'study_equal_weight_fraction_emission_dominant_cm1':round(mean(
            g['fraction_em_gt_abs_cm1'] for g in groups),4) if groups else None,
       'top_tags':state['source_tags'].most_common(5)}
result={
  'source':'FULL stage2_all_endpoints.csv (not seven selected contrasts)',
  'input_rows':n['input_rows'],'unique_molecules_scanned':len(mols),
  'candidate_rows':n['candidate_rows'],'candidate_solvent_pairs':len(pairs),
  'by_motif':summary,
  'limitations':[
   'Structural motif screening only; not verified dye class or fluorescent chemical state',
   'All-carbon tricycles are not automatically carbon xanthene analogues',
   'Same DOI does not guarantee controlled pH, temperature or speciation',
   'Each molecule-DOI group receives equal weight in descriptive summary',
   'No DFT descriptors or quantum-property validation performed']}
(R/'stage11_full_database_audit.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
print('FINISHED')
