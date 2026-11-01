import csv, json, math, hashlib, os
from pathlib import Path
from collections import Counter, defaultdict
from statistics import median
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold
RDLogger.DisableLog('rdApp.error')

D = Path.home()/'Downloads'
SRC = Path(os.environ.get('FLUODB_CSV', str(D/'FluoDB.csv')))
OUT = Path(os.environ.get('XANPY_ROOT', str(D/'XanPy_QM_Atlas')))
OUT.mkdir(parents=True,exist_ok=True)

if not SRC.is_file():
    raise SystemExit('FluoDB.csv not found in ~/Downloads')

def numeric(s):
    try:
        x=float(str(s).strip())
        return x if math.isfinite(x) else None
    except (ValueError, TypeError): return None

def allowed(s, low, high):
    x=numeric(s)
    return x if x is not None and low<=x<=high else None

def yes(x): return x is not None

def blank(x): return '' if x is None else x

fields = ['source_row','mol_id','canonical_smiles','murcko_scaffold','solvent','solvent_key',
          'reference_doi','source_database','source_tag','source_tag_name','original_split',
          'absorption_nm','emission_nm','plqy_fraction','molar_extinction_m_1_cm_1',
          'stokes_cm_1','qc_flags']
counts=Counter(); unique=set(); solvents=defaultdict(lambda:defaultdict(list)); fam=Counter(); solvent_n=Counter()
pair_vals=defaultdict(list)

with SRC.open(newline='', encoding='utf-8-sig', errors='replace') as fi, (OUT/'stage2_all_endpoints.csv').open('w',newline='',encoding='utf-8') as fo:
    reader=csv.DictReader(fi)
    expected={'smiles','solvent','absorption/nm','emission/nm','plqy','e/m-1cm-1'}
    if not expected.issubset(set(reader.fieldnames or [])):
        raise SystemExit('Missing columns: '+str(expected-set(reader.fieldnames or [])))
    wr=csv.DictWriter(fo,fieldnames=fields); wr.writeheader()
    for idx,r in enumerate(reader,2):
        counts['source_rows']+=1
        mol=Chem.MolFromSmiles((r.get('smiles') or '').strip())
        if mol is None:
            counts['invalid_smiles']+=1; continue
        if len(Chem.GetMolFrags(mol))!=1:
            counts['multicomponent']+=1; continue
        canonical=Chem.MolToSmiles(Chem.RemoveHs(mol),isomericSmiles=True)
        a=allowed(r.get('absorption/nm'),200,1500)
        e=allowed(r.get('emission/nm'),200,1500)
        p=allowed(r.get('plqy'),0,1)
        ext=allowed(r.get('e/m-1cm-1'),0,1e8)
        raw_plqy=numeric(r.get('plqy'))
        if raw_plqy is not None and not 0<=raw_plqy<=1: counts['plqy_outside_0_1']+=1
        if all(v is None for v in (a,e,p,ext)):
            counts['no_usable_endpoint']+=1; continue
        if a is None and e is None and (p is not None or ext is not None):
            counts['recovered_plqy_or_epsilon_only']+=1
        flags=[]
        solvent=(r.get('solvent') or '').strip(); solvent_mol=Chem.MolFromSmiles(solvent) if solvent else None
        solvent_key=(Chem.MolToSmiles(solvent_mol,isomericSmiles=True) if solvent_mol and len(Chem.GetMolFrags(solvent_mol))==1 else '')
        if not solvent_key: flags.append('missing_solvent')
        if a is not None and e is not None and e<a:
            flags.append('negative_stokes'); counts['negative_stokes']+=1
        # Flag invalid or ambiguous PLQY values without converting percentages automatically
        if raw_plqy is not None and not 0<=raw_plqy<=1: flags.append('plqy_out_of_range')
        if a is not None: counts['has_absorption']+=1
        if e is not None: counts['has_emission']+=1
        if p is not None: counts['has_plqy']+=1
        if ext is not None: counts['has_extinction']+=1
        scaffold=Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))
        compound_id=hashlib.sha1(canonical.encode()).hexdigest()[:16]
        stokes=1e7*(1/a-1/e) if a is not None and e is not None else None
        tag=(r.get('tag_name') or '').strip() or 'unclassified'
        wr.writerow({'source_row':idx,'mol_id':compound_id,'canonical_smiles':canonical,
          'murcko_scaffold':scaffold,'solvent':solvent,'solvent_key':solvent_key,
          'reference_doi':r.get('reference(doi)') or '',
          'source_database':r.get('source') or '',
          'source_tag':r.get('tag') or '', 'source_tag_name':r.get('tag_name') or '',
          'original_split':r.get('split') or '',
          'absorption_nm':blank(a),'emission_nm':blank(e),'plqy_fraction':blank(p),
          'molar_extinction_m_1_cm_1':blank(ext),
          'stokes_cm_1':blank(round(stokes,3) if stokes is not None else None), 'qc_flags':';'.join(flags)})
        counts['retained_rows']+=1
        unique.add(canonical); fam[tag]+=1
        if solvent_key: solvent_n[solvent_key]+=1
        if e is not None and solvent_key:
            solvents[canonical][solvent_key].append(e)
            pair_vals[(canonical,solvent_key)].append(e)

with (OUT/'stage2_family_counts.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(['source_tag_name','retained_rows']);w.writerows(fam.most_common())
with (OUT/'stage2_solvent_counts.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(['solvent','retained_rows']);w.writerows(solvent_n.most_common())
with (OUT/'stage2_same_molecule_solvents.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(['canonical_smiles','n_solvents','min_median_emission_nm','max_median_emission_nm','emission_range_nm'])
    for smi,s in solvents.items():
        if len(s)<2: continue
        v=[median(x) for x in s.values()]
        w.writerow([smi,len(s),min(v),max(v),round(max(v)-min(v),2)])

n_multi=sum(len(s)>1 for s in solvents.values())
audit={'input_source_rows':counts['source_rows'],'counts':dict(counts),
       'unique_molecules_with_any_endpoint':len(unique),
       'unique_molecule_solvent_with_emission':len(pair_vals),
       'molecules_with_multiple_emission_solvents':n_multi,
       'top_families':fam.most_common(12),'top_solvents':solvent_n.most_common(12),
       'rules':{'optical_nm':'200 to 1500 inclusive','plqy':'0 to 1 (no automatic percentage conversion)',
         'extinction':'0 to 1e8, units from source header',
         'molecule_id':'isomeric canonical SMILES; charge and stereochemistry retained',
         'duplicate_policy':'retain source rows/DOIs; do not silently merge',
         'stokes':'computed only if absorption and emission both on same source row'},
       'source_file':'FluoDB.csv','source_sha256':hashlib.sha256(SRC.read_bytes()).hexdigest(),
       'note':'Results describe the exact GitHub CSV downloaded; not the entire publication release.'}
(OUT/'stage2_audit.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(audit,indent=2,ensure_ascii=False))
print('FINISHED:',OUT)
