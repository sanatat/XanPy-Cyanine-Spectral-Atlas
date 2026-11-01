#!/usr/bin/env python3
"""XanPy Stage 16: conservative structural test of the cyanine xi-sign rule.

This is NOT a validation from spectral coincidence. Only structurally matched
single-substituent molecular pairs are compared, with same-paper results
strictly separated from across-paper exploratory results.
"""
import argparse, collections, csv, json, math, re
from pathlib import Path
import numpy as np
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.error')

DOI_RE = re.compile(r'10\.\d{4,9}/[^\s,;]+', re.I)
SIGMA = {
    'H': 0.0, 'C': -0.17, 'CO': -0.27, 'F': 0.06,
    'Cl': 0.23, 'Br': 0.23, 'I': 0.18, 'C#N': 0.66,
    'N': -0.66, 'C(F)(F)F': 0.65, '[N+](=O)[O-]': 0.78,
}
LENGTHS = (3,5,7,9)

def doi(s):
    m=DOI_RE.search(s or '')
    return m.group().rstrip(').,;').lower() if m else ''

def valid_float(s):
    try:
        v=float(s)
        return v if math.isfinite(v) else None
    except (TypeError,ValueError):
        return None

def ring_systems(m):
    rings=m.GetRingInfo().AtomRings()
    parent=list(range(m.GetNumAtoms()))
    def find(x):
        while parent[x]!=x:
            parent[x]=parent[parent[x]];x=parent[x]
        return x
    def union(a,b):
        a=find(a);b=find(b)
        if a!=b:parent[b]=a
    for ring in rings:
        for a in ring[1:]:union(ring[0],a)
    # Only atoms that participate in rings get a ring system
    member={a for ring in rings for a in ring}
    return {a:find(a) for a in member}

def connected_component(m,start,forbidden):
    if start in forbidden:return set()
    q=[start]; out={start}
    for i in q:
        for nb in m.GetAtomWithIdx(i).GetNeighbors():
            j=nb.GetIdx()
            if j not in forbidden and j not in out:
                out.add(j);q.append(j)
    return out

def chain_candidates(m):
    sys=ring_systems(m)
    if not sys:return []
    sys_atoms=collections.defaultdict(set)
    for at,s in sys.items():sys_atoms[s].add(at)
    # A conservative aromatic or heterocyclic auxochrome endpoint.
    valid_system={s for s,ids in sys_atoms.items() if any(m.GetAtomWithIdx(i).GetSymbol() in ('N','O','S','Se') for i in ids)}
    found=[]
    # Start on ring carbons joined to a non-ring carbon, then walk only
    # NON-RING sp2 carbon atoms; no interposed heteroatom or ring atom.
    for start,s in sys.items():
        if s not in valid_system or m.GetAtomWithIdx(start).GetAtomicNum()!=6:continue
        for first in m.GetAtomWithIdx(start).GetNeighbors():
            f=first.GetIdx()
            if f in sys or first.GetAtomicNum()!=6:continue
            stack=[(f,[f],{f,start})]
            while stack:
                at,path,visited=stack.pop()
                if len(path)>9:continue
                for nb in m.GetAtomWithIdx(at).GetNeighbors():
                    j=nb.GetIdx()
                    if j in visited:continue
                    if j in sys:
                        if sys[j]==s or sys[j] not in valid_system or nb.GetAtomicNum()!=6:continue
                        if len(path) not in LENGTHS:continue
                        chain=[m.GetAtomWithIdx(x) for x in path]
                        if not all(c.GetHybridization() == Chem.HybridizationType.SP2 for c in chain):continue
                        nodes=[start]+path+[j]
                        bonds=[m.GetBondBetweenAtoms(x,y) for x,y in zip(nodes,nodes[1:])]
                        if any(b.GetBondType()==Chem.BondType.TRIPLE for b in bonds):continue
                        db=sum(b.GetBondType()==Chem.BondType.DOUBLE for b in bonds)
                        if db < (len(path)+1)//2:continue
                        orient=(start,tuple(path),j)
                        rev=(j,tuple(reversed(path)),start)
                        cand=min(orient,rev)
                        found.append((cand,db))
                    elif nb.GetAtomicNum()==6 and len(path)<9:
                        stack.append((j,path+[j],visited|{j}))
    distinct={x[0]:x for x in found}
    return list(distinct.values())

def substitutions(m,start,path,end):
    forbidden=set(path)|{start,end}
    result=[]
    for i,p in enumerate(path):
        atom=m.GetAtomWithIdx(p)
        other=[nb.GetIdx() for nb in atom.GetNeighbors() if nb.GetIdx() not in forbidden]
        if not other:
            if atom.GetTotalNumHs() == 0:return None
            result.append('H');continue
        if len(other)!=1:return None
        comp=connected_component(m,other[0], forbidden)
        if len(comp)>22:return None
        # Substituent must NOT be reconnected to other path atoms or ring endpoint;
        # this would be a cyclized polymethine and needs a different template.
        if any(m.GetAtomWithIdx(k).GetIdx() in (start,end) for k in comp):return None
        frag=Chem.MolFragmentToSmiles(m,atomsToUse=sorted(comp),isomericSmiles=True)
        result.append(frag)
    return result

def characterize(smi):
    m=Chem.MolFromSmiles(smi)
    if m is None:return None,'invalid_smiles'
    if len(Chem.GetMolFrags(m))!=1:return None,'multi_fragment'
    cand=chain_candidates(m)
    if not cand:return None,'no_strict_bridge'
    # Only allow a unique structurally distinguishable polymethine chain.
    unique={tuple(c[0][1]) for c in cand}
    # Reversed orientations are normalized; alternatives of same atoms rejected.
    if len(unique)!=1:return None,'ambiguous_chains'
    ((start,path,end),db)=cand[0]
    R=substitutions(m,start,path,end)
    if R is None:return None,'cyclized_or_complex_substitution'
    # Break the chain and retain EXACT terminal heterocycle fragments.
    removed=set(path)
    a=connected_component(m,start,removed)
    b=connected_component(m,end,removed)
    if a & b:return None,'terminal_groups_connected_elsewhere'
    end1=Chem.MolFragmentToSmiles(m,atomsToUse=sorted(a),isomericSmiles=False)
    end2=Chem.MolFragmentToSmiles(m,atomsToUse=sorted(b),isomericSmiles=False)
    # Re-orient all substituent sites to a canonical terminal order.
    if end1>end2:
        end1,end2=end2,end1
        R=R[::-1]
    elif end1==end2:
        R=min(R,R[::-1])
    n=len(path)
    xi=None
    if all(g in SIGMA for g in R):
        factor=1 if n%4==1 else -1
        xi=round(sum(factor*((-1)**i)*SIGMA[g] for i,g in enumerate(R)),4)
    return {'n':n,'endgroup_signature':end1+' :: '+end2,
            'R':R,'xi':xi,'n_double_bonds':db},'ok'

def aggregate_csv(path):
    with path.open(newline='',encoding='utf-8-sig') as f:
        for r in csv.DictReader(f):yield r

def dump_csv(path,rows,fields):
    with path.open('w',newline='',encoding='utf-8') as f:
        wr=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');wr.writeheader();wr.writerows(rows)

def main():
    ap=argparse.ArgumentParser()
    home=Path.home()/'Downloads'/'XanPy_QM_Atlas'
    ap.add_argument('--input',default=str(home/'stage2_all_endpoints.csv'))
    ap.add_argument('--output',default=str(home))
    args=ap.parse_args(); inp=Path(args.input);out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    if not inp.exists():raise SystemExit('Input file missing: '+str(inp))
    stats=collections.Counter();raw=collections.defaultdict(list)
    source_smiles=set()
    for r in aggregate_csv(inp):
        if (r.get('source_tag_name') or '').strip().casefold()!='cyanine':continue
        stats['cyanine_rows']+=1
        s=(r.get('canonical_smiles') or '').strip()
        if not s:continue
        source_smiles.add(s)
        v=valid_float(r.get('absorption_nm'))
        if v is None or not 200<=v<=1500:continue
        sol=(r.get('solvent') or '').strip()
        paper=doi(r.get('reference_doi'))
        if not sol:continue
        # Keep only existing canonical solvent codes; case-sensitive SMILES.
        sm=Chem.MolFromSmiles(sol)
        if sm is None or len(Chem.GetMolFrags(sm)) != 1:continue
        sol=Chem.MolToSmiles(sm,isomericSmiles=True)
        raw[(s,sol,paper)].append(v)
    ann={};reason=collections.Counter()
    for smi in sorted(source_smiles):
        result,why=characterize(smi)
        reason[why]+=1
        if result is not None:ann[smi]=result
    annotations=[]
    for smi in sorted(source_smiles):
        d=ann.get(smi)
        annotations.append({'canonical_smiles':smi,'status':'strict_candidate' if d else 'unclassified',
                'chain_length':d['n'] if d else '',
                'terminal_signature':d['endgroup_signature'] if d else '',
                'chain_substituents':' | '.join(d['R']) if d else '',
                'xi_SI_simple_groups_only':d['xi'] if d and d['xi'] is not None else '',
                'n_double_bonds':d['n_double_bonds'] if d else ''})
    dump_csv(out/'stage16_strict_chain_annotations.csv',annotations,list(annotations[0]))
    stats['cyanine_molecules_with_absorption'] = len({key[0] for key in raw})
    meds={k:float(np.median(v)) for k,v in raw.items() if k[0] in ann and ann[k[0]]['xi'] is not None}
    # Collapsed to one median peak per molecule and solvent, for cross-publication exploratory comparisons.
    across=collections.defaultdict(list)
    for (s,sol,p),v in meds.items():across[(s,sol)].append((p,v))
    molecule_solvent={(s,sol):float(np.median([v for _,v in vals])) for (s,sol),vals in across.items()}
    def matched(one,two):
        s1,sol,p1,y1=one;s2,_,p2,y2=two
        if s1==s2:return None
        a,b=ann[s1],ann[s2]
        if a['n']!=b['n'] or a['endgroup_signature']!=b['endgroup_signature']:return None
        differences=[i for i,(x,y) in enumerate(zip(a['R'],b['R'])) if x!=y]
        if len(differences)!=1:return None
        j=differences[0]; dxi=b['xi']-a['xi']
        if abs(dxi)<0.045:return None
        dy=y2-y1
        if abs(dy)<1e-8:return None
        predicted=math.copysign(1, (1 if a['n']%4==3 else -1)*dxi)
        return {'n':a['n'],'position':j+1,'meso':int(j==a['n']//2),
                'smiles_A':s1,'smiles_B':s2,'R_A':a['R'][j],'R_B':b['R'][j],
                'xi_A':a['xi'],'xi_B':b['xi'],'delta_xi':round(dxi,3),
                'solvent':sol,'doi_A':p1,'doi_B':p2,
                'abs_A_nm':round(y1,3),'abs_B_nm':round(y2,3),'observed_delta_nm':round(dy,3),
                'predicted_delta_direction':'red' if predicted>0 else 'blue',
                'observed_delta_direction':'red' if dy>0 else 'blue',
                'direction_agrees':int(predicted*dy>0),
                'endgroup_signature':a['endgroup_signature']}
    def pairs(mode):
        sets=collections.defaultdict(list)
        if mode=='same_Doi':
            for (s,sol,p),y in meds.items():
                if not p:continue
                d=ann[s]
                sets[(d['n'],d['endgroup_signature'],sol,p)].append((s,sol,p,y))
        else:
            for (s,sol),y in molecule_solvent.items():
                d=ann[s]
                papers='|'.join(sorted({p for p,_ in across[(s,sol)] if p}))
                sets[(d['n'],d['endgroup_signature'],sol)].append((s,sol,papers,y))
        found=[]
        for key,items in sets.items():
            items.sort(key=lambda x:x[0])
            for i,a in enumerate(items):
                for b in items[i+1:]:
                    if mode=='different_Doi':
                        pa={x for x in a[2].split('|') if x}
                        pb={x for x in b[2].split('|') if x}
                        if not pa or not pb or pa.intersection(pb):continue
                    result=matched(a,b)
                    if result:found.append(result)
        return found
    strict=pairs('same_Doi')
    exploratory=pairs('different_Doi')
    fields=['n','position','meso','smiles_A','smiles_B','R_A','R_B','xi_A','xi_B','delta_xi','solvent','doi_A','doi_B','abs_A_nm','abs_B_nm','observed_delta_nm','predicted_delta_direction','observed_delta_direction','direction_agrees','endgroup_signature']
    dump_csv(out/'stage16_same_DOI_single_substituent_pairs.csv',strict,fields)
    dump_csv(out/'stage16_different_DOI_same_solvent_EXPLORATORY.csv',exploratory,fields)
    def summary(pairs):
        result={}
        for n in LENGTHS:
            subset=[r for r in pairs if r['n']==n]
            # Also count unique molecule pairs; pairs may share molecules and are not independent.
            unique={tuple(sorted((r['smiles_A'],r['smiles_B']))) for r in subset}
            papers=set(r['doi_A'] for r in subset if r['doi_A']==r['doi_B'] and r['doi_A'])
            result[f'Cy-{n}']={
              'matched_comparisons':len(subset),
              'unique_molecular_pairs':len(unique),
              'distinct_same_paper_DOIs':len(papers),
              'direction_agreement_fraction':round(sum(r['direction_agrees'] for r in subset)/len(subset),4) if subset else None,
              'status':'EXPLORATORY_ONLY_not_independent_validation' if subset else 'NO_MATCHED_PAIRS'
            }
        return result
    length_counts=collections.Counter(d['n'] for d in ann.values())
    fully=collections.Counter(d['n'] for d in ann.values() if d['xi'] is not None)
    report={
        'status':'STRUCTURAL_MATCHING_RESULTS_NOT_FORMAL_INDEPENDENT_VALIDATION',
        'cyanine_tagged_rows':stats['cyanine_rows'],
        'cyanine_tagged_unique_smiles':len(source_smiles),
        'cyanine_tagged_unique_smiles_with_absorption':stats['cyanine_molecules_with_absorption'],
        'strict_chain_candidates_by_length':{str(n):length_counts[n] for n in LENGTHS},
        'xi_annotated_candidates_by_length':{str(n):fully[n] for n in LENGTHS},
        'rejection_reasons':dict(reason),
        'strict_same_DOI_single_substitution_test':summary(strict),
        'different_DOI_same_solvent_exploratory_test':summary(exploratory),
        'limitations':[
           'Only strict noncyclic sp2 carbon bridges between two distinct heterocyclic ring systems are recognized.',
           'Many true cyanines may be excluded; some candidate chains may still be incorrectly classified.',
           'No automatic identity transfer from spectral matching with the 2026 SI.',
           'Sigma values are approximate published simple-substituent Hammett proxies, not inferred quantum parameters.',
           'Both pair comparisons and DOI groups can be dependent; no p-values or formal confirmation reported.',
           'Across-paper pairs confounded by protocols, pH, and aggregation; never treat them as strict tests.'
        ]}
    (out/'stage16_direct_test_audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2));print('FINISHED')

if __name__=='__main__':main()
