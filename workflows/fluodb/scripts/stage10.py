import csv, json, re, os
from pathlib import Path
from collections import Counter, defaultdict
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog('rdApp.error')
ROOT=Path(os.environ.get('XANPY_ROOT', str(Path.home()/'Downloads'/'XanPy_QM_Atlas')))
SOURCE=ROOT/'stage2_all_endpoints.csv'
PAIRS=ROOT/'stage8_selected_solvent_pairs.csv'
DOI=re.compile(r'10\.\d{4,9}/[^\s,;]+', re.I)

def doi_key(v):
    m=DOI.search(v or '')
    return m.group().rstrip(').,;').lower() if m else ''

def core_group(mol):
    """Structural candidate motif, not a validated dye class."""
    rings=[set(x) for x in mol.GetRingInfo().AtomRings() if len(x)==6]
    for center in rings:
        aromatic_flanks=[r for r in rings if r != center and len(r & center)==2 and all(
            mol.GetAtomWithIdx(i).GetAtomicNum()==6 and mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
        for i,left in enumerate(aromatic_flanks):
            for right in aromatic_flanks[i+1:]:
                if (left & right)-center: continue
                bridge=center-(left|right)
                if len(bridge)!=2: continue
                atoms=[mol.GetAtomWithIdx(j) for j in bridge]
                symbols=[a.GetSymbol() for a in atoms]
                if 'O' in symbols and 'C' in symbols:
                    return 'xanthene_like_O_candidate'
                if 'Si' in symbols and 'C' in symbols:
                    return 'silicon_xanthene_like_candidate'
                if 'Ge' in symbols and 'C' in symbols:
                    return 'germanium_xanthene_like_candidate'
                if symbols.count('C')==2:
                    for a in atoms:
                        exo=[n for n in a.GetNeighbors() if n.GetIdx() not in center]
                        if len(exo)>=2 and sum(n.GetAtomicNum()==6 for n in exo)>=2:
                            return 'geminal_C_bridged_tricyclic_candidate'
                    return 'carbon_bridged_tricyclic_candidate'
    return 'other_structure'

# Stage 2 source tags are descriptive only; keep ambiguity explicit.
labels=defaultdict(Counter)
with SOURCE.open(newline='',encoding='utf-8') as f:
    for r in csv.DictReader(f):
        smi=r.get('canonical_smiles','').strip()
        doi=doi_key(r.get('reference_doi'))
        if smi and doi:
            labels[(smi,doi)][(r.get('source_tag_name') or 'unclassified').strip()] += 1

mols={}
def features(smi):
    if smi not in mols:
        mol=Chem.MolFromSmiles(smi)
        if mol is None:
            mols[smi]=('invalid','',0)
        else:
            try:
                scaffold=Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol),isomericSmiles=True)
                mols[smi]=(core_group(mol),scaffold,len(mol.GetRingInfo().AtomRings()))
            except Exception:
                mols[smi]=('failed','',0)
    return mols[smi]

rows=[]
with PAIRS.open(newline='',encoding='utf-8') as f:
    for r in csv.DictReader(f):
        smi=r['canonical_smiles'];doi=doi_key(r['doi']);t=labels.get((smi,doi),Counter())
        label=next(iter(t)) if len(t)==1 else ('ambiguous' if t else 'unmatched')
        category,scaffold,n_rings=features(smi)
        try:
            em=float(r['delta_em_nm']);ab=float(r['delta_abs_nm']);st=float(r['delta_stokes_cm1'])
        except (TypeError,ValueError): continue
        if not all(np.isfinite(x) for x in (em,ab,st)): continue
        rows.append({'contrast':r['solvent_A']+' -> '+r['solvent_B'],
                     'tag':label,'core':category,'scaffold':scaffold,
                     'smiles':smi,'doi':doi,'em':em,'ab':ab,'st':st})

def stats(group):
    bymol=defaultdict(list);bydoi=defaultdict(list)
    for x in group:
        bymol[x['smiles']].append(x['em'])
        bydoi[x['doi']].append(x['em'])
    molv=np.array([np.mean(v) for v in bymol.values()])
    doiv=np.array([np.mean(v) for v in bydoi.values()])
    return dict(n_pairs=len(group),n_molecules=len(bymol),n_dois=len(bydoi),
                mean_em_nm=round(float(molv.mean()),3),
                doi_balanced_mean_em_nm=round(float(doiv.mean()),3),
                median_abs_em_nm=round(float(np.median(np.abs(molv))),3))

def dump(path,data):
    if not data: return
    with (ROOT/path).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]))
        w.writeheader();w.writerows(data)

# Report motif prevalence within each source tag and solvent contrast.
bytag=defaultdict(list);bycore=defaultdict(list);byscaffold=defaultdict(list)
for x in rows:
    bytag[(x['contrast'],x['tag'])].append(x)
    bycore[(x['contrast'],x['tag'],x['core'])].append(x)
    byscaffold[(x['contrast'],x['scaffold'])].append(x)

bytag_results=[dict(contrast=a,source_tag=b,**stats(v)) for (a,b),v in bytag.items()]
core_results=[dict(contrast=a,source_tag=b,structure_motif=c,**stats(v))
              for (a,b,c),v in bycore.items() if len({x['smiles'] for x in v})>=5]
scaffold_results=[dict(contrast=a,scaffold=c,**stats(v)) for (a,c),v in byscaffold.items()
                  if len({x['smiles'] for x in v})>=5 and len({x['doi'] for x in v})>=3]
core_results.sort(key=lambda r:(r['contrast'],r['source_tag'],-r['n_molecules']))
scaffold_results.sort(key=lambda r:-r['n_molecules'])
dump('stage10_tag_summary.csv',bytag_results)
dump('stage10_structural_motifs.csv',core_results)
dump('stage10_frequent_scaffolds.csv',scaffold_results)

# A compact inspection table of the most striking broad '6n6' group.
interesting=[x for x in rows if x['contrast']=='toluene -> acetonitrile' and x['tag']=='6n6']
interesting.sort(key=lambda x:-abs(x['em']))
examples=[{k:x[k] for k in ('smiles','doi','tag','core','scaffold','em','ab','st')}
          for x in interesting[:20]]
dump('stage10_6n6_structure_examples.csv',examples)

# Export all xanthene-like candidate pairs for follow-up chemical verification.
candidates=[{k:x[k] for k in ('contrast','tag','core','scaffold','smiles','doi','em','ab','st')}
            for x in rows if ('xanthene_like' in x['core'] or 'bridged_tricyclic' in x['core'])]
dump('stage10_tricyclic_candidates.csv',candidates)

report={'total_pairs':len(rows),'unique_molecules':len(mols),
        'source_tag_6n6_interpretation':'Broad 6+n+6 scaffold templates, NOT validated xanthene dyes',
        'motif_candidates':dict(Counter(x['core'] for x in rows)),
        'xanthene_or_bridged_tricyclic_candidate_pairs':len(candidates),
        'top_toluene_acetonitrile_tags':sorted([r for r in bytag_results
             if r['contrast']=='toluene -> acetonitrile'],key=lambda x:-x['n_pairs'])[:15],
        'n_frequent_scaffold_rows':len(scaffold_results),
        'limitations':['Structural motifs are screening heuristics requiring inspection',
          'Murcko scaffolds can split the same photophysical dye family',
          'No claim about mechanism or solvent causality',
          'Source tags overlap broad structural categories']}
(ROOT/'stage10_audit.json').write_text(json.dumps(report,indent=2,ensure_ascii=False))
print(json.dumps(report,indent=2,ensure_ascii=False))
print('FINISHED')
