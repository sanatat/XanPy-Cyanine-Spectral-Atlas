#!/usr/bin/env python3
"""XanPy-QM Atlas Stage 12: topology / oxidation-state / provenance audit.
Purely structural screening: does not verify fluorescence mechanisms or DFT.
"""
import csv, json, re, math, os
from pathlib import Path
from collections import Counter, defaultdict
from itertools import combinations
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.error')

ROOT = Path(os.environ.get('XANPY_ROOT', str(Path.home()/'Downloads/XanPy_QM_Atlas')))
SRC = ROOT/'stage2_all_endpoints.csv'
if not SRC.exists(): raise SystemExit(f'Missing: {SRC}')
DOI = re.compile(r'10\.\d{4,9}/[^\s,;]+', re.I)

def doi(s):
    m = DOI.search(s or '')
    return m.group().rstrip(').,;').lower() if m else ''

def valid_wave(s):
    try:
        x=float(s)
        return x if math.isfinite(x) and 200<=x<=1500 else None
    except (ValueError,TypeError):return None

def core(m):
    """Find two 6-membered aromatic carbon rings fused to a central 6-ring;
    bridge X and C must each join the two opposite aromatic flanks.
    Classify the bridge C's local bonding, without claiming actual dye identity.
    """
    if m is None:return 'invalid_smiles'
    rings=[set(r) for r in m.GetRingInfo().AtomRings() if len(r)==6]
    found=[]
    for center in rings:
        flanks=[r for r in rings if r!=center and len(center&r)==2 and all(
            m.GetAtomWithIdx(i).GetAtomicNum()==6 and
            m.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
        for a,b in combinations(flanks,2):
            if (a&b)-center:continue
            bridge=center-a-b
            if len(bridge)!=2:continue
            elems={m.GetAtomWithIdx(i).GetSymbol():i for i in bridge}
            if len(elems)!=2 or 'C' not in elems:continue
            hetero=next((x for x in ('O','Si','Ge') if x in elems),None)
            if hetero is None:continue
            x,c=elems[hetero],elems['C']
            shared_a=a&center;shared_b=b&center
            def joins_both(idx):
                neighbors={n.GetIdx() for n in m.GetAtomWithIdx(idx).GetNeighbors()}
                return bool(neighbors&shared_a) and bool(neighbors&shared_b)
            if not (joins_both(x) and joins_both(c)):continue
            if m.GetBondBetweenAtoms(x,c):continue
            ca=m.GetAtomWithIdx(c)
            has_carbonyl=any(
                n.GetAtomicNum()==8 and n.GetIdx() not in center and
                m.GetBondBetweenAtoms(c,n.GetIdx()).GetBondType()==Chem.BondType.DOUBLE
                for n in ca.GetNeighbors())
            state=('carbonyl_9C' if has_carbonyl else
                   'sp3_9C' if ca.GetHybridization()==Chem.HybridizationType.SP3 else
                   'conjugated_9C' if ca.GetHybridization()==Chem.HybridizationType.SP2 else
                   'other_9C')
            found.append(hetero+'_bridge__'+state)
    if not found:return 'outside_strict_fused_core'
    return '|'.join(sorted(set(found))) if len(set(found))>1 else found[0]

REFERENCE={
    'xanthene':('C1C2=CC=CC=C2OC3=CC=CC=C31','O_bridge__sp3_9C'),
    'xanthone':('O=C1c2ccccc2Oc2ccccc21','O_bridge__carbonyl_9C'),
    'rhodamine_B_example':('CCN(CC)C1=CC2=C(C=C1)C(=C3C=CC(=[N+](CC)CC)C=C3O2)C4=CC=CC=C4C(=O)O','O_bridge__conjugated_9C'),
    'Si_analogue':('c1ccc2c(c1)[SiH2]c1ccccc1C2','Si_bridge__sp3_9C'),
    'Ge_analogue':('c1ccc2c(c1)[GeH2]c1ccccc1C2','Ge_bridge__sp3_9C'),
    'acridine':('c1ccc2nc3ccccc3cc2c1','outside_strict_fused_core'),
    'dibenzofuran':('c1ccc2oc3ccccc3c2c1','outside_strict_fused_core'),
    'benzene':('c1ccccc1','outside_strict_fused_core'),
}
for name,(smi,expected) in REFERENCE.items():
    actual=core(Chem.MolFromSmiles(smi))
    if actual!=expected:raise SystemExit(f'Self-test FAILED {name}: {actual}')

molecule_class={}
counts=defaultdict(Counter)
structures=defaultdict(set)
studies=defaultdict(lambda:defaultdict(set))
tags=defaultdict(Counter)
source_names=defaultdict(Counter)
examples={}
rows=[]
with SRC.open(newline='',encoding='utf-8') as fi:
    for r in csv.DictReader(fi):
        smi=(r.get('canonical_smiles') or '').strip()
        if smi not in molecule_class:
            molecule_class[smi]=core(Chem.MolFromSmiles(smi) if smi else None)
        label=molecule_class[smi]
        if label in ('outside_strict_fused_core','invalid_smiles'):continue
        c=counts[label];c['rows']+=1;structures[label].add(smi)
        paper=doi(r.get('reference_doi'))
        if paper:c['doi_rows']+=1
        else:c['missing_doi_rows']+=1
        a=valid_wave(r.get('absorption_nm'));e=valid_wave(r.get('emission_nm'))
        if a is not None:c['has_absorption']+=1
        if e is not None:c['has_emission']+=1
        if a is not None and e is not None and e>=a:
            c['valid_paired_spectra']+=1
            if paper:
                c['paired_spectra_with_doi']+=1
                solvent=(r.get('solvent') or '').strip()
                if solvent:
                    sm=Chem.MolFromSmiles(solvent)
                    if sm and len(Chem.GetMolFrags(sm))==1:
                        studies[label][(smi,paper)].add(Chem.MolToSmiles(sm,isomericSmiles=True))
        tag=(r.get('source_tag_name') or 'unclassified').strip() or 'unclassified'
        src=(r.get('source_database') or 'unclassified').strip() or 'unclassified'
        tags[label][tag]+=1;source_names[label][src]+=1
        record={'chemical_state':label,'smiles':smi,'source_row':r.get('source_row',''),
          'source_tag':tag,'source_database':src,'doi':paper,
          'solvent':r.get('solvent',''),'absorption_nm':r.get('absorption_nm',''),
          'emission_nm':r.get('emission_nm',''),'plqy':r.get('plqy_fraction','')}
        rows.append(record)
        ex_key=(label,tag)
        if paper and (ex_key not in examples or not examples[ex_key]['doi']):examples[ex_key]=record
        elif ex_key not in examples:examples[ex_key]=record

headers=['chemical_state','smiles','source_row','source_tag','source_database','doi',
         'solvent','absorption_nm','emission_nm','plqy']
def write_csv(path, records):
    with (ROOT/path).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=headers);w.writeheader();w.writerows(records)
write_csv('stage12_candidate_records.csv',rows)
write_csv('stage12_review_examples.csv',list(examples.values()))

by_core={}
for label,c in sorted(counts.items()):
    group=studies[label]
    multi={k:sols for k,sols in group.items() if len(sols)>1}
    by_core[label]={**dict(c),'unique_molecules':len(structures[label]),
      'unique_dois_with_paired_spectra':len({key[1] for key in group}),
      'multi_solvent_molecule_doi_groups':len(multi),
      'multi_solvent_unique_molecules':len({key[0] for key in multi}),
      'top_tags':tags[label].most_common(6),
      'top_sources':source_names[label].most_common(6)}

report={'input_rows':sum(1 for _ in open(SRC,encoding='utf-8'))-1,
  'unique_structures_scanned':len(molecule_class),
  'strict_reference_tests_passed':len(REFERENCE),
  'candidate_records':len(rows),'categories':by_core,
  'cautions':['Source tag Acridines is not a verified chemical identity',
  'Carbonyl-containing xanthone-like and conjugated xanthene dyes are separated',
  'Only valid paired spectra with DOI and >=2 distinct solvents count as multisolvent studies',
  'This is topology/oxidation-state screening, not verified fluorophore identity or DFT'],
  'files':['stage12_candidate_records.csv','stage12_review_examples.csv']}
(ROOT/'stage12_audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
print('FINISHED')
