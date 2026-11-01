#!/usr/bin/env python3
"""Diagnose why strict Cy-3/5/7/9 matched pairs are absent.
This performs NO hypothesis testing and deliberately does not loosen chemical identity rules.
"""
import argparse
import csv
import json
import math
import re
from collections import defaultdict, Counter
from pathlib import Path
from statistics import median
from itertools import combinations
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.error')

DOI = re.compile(r'10\.\d{4,9}/[^\s,;]+',re.I)
def doi(s):
    m=DOI.search(s or '')
    return m.group().rstrip(').,;').lower() if m else ''
def finite(v):
    try:
        x=float(v)
        return x if math.isfinite(x) and 200<=x<=1500 else None
    except (TypeError,ValueError):return None

def load(path):
    with path.open(newline='',encoding='utf-8-sig') as f:
        yield from csv.DictReader(f)

def main():
    a=argparse.ArgumentParser(description=__doc__)
    home=Path.home()/'Downloads'/'XanPy_QM_Atlas'
    a.add_argument('--data',type=Path,default=home/'stage2_all_endpoints.csv')
    a.add_argument('--annotations',type=Path,default=home/'stage16_strict_chain_annotations.csv')
    a.add_argument('--out',type=Path,default=home)
    args=a.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    for p in [args.data,args.annotations]:
        if not p.exists():raise SystemExit(f'File missing: {p}')

    anns={}
    for r in load(args.annotations):
        if r.get('status')!='strict_candidate':continue
        n=int(r['chain_length'])
        if n not in (3,5,7,9):continue
        xi=r.get('xi_SI_simple_groups_only','').strip()
        if not xi:continue
        try: xi=float(xi)
        except ValueError:continue
        R=tuple(r['chain_substituents'].split(' | '))
        if len(R)!=n:continue
        anns[r['canonical_smiles']]={'n':n,'xi':xi,'R':R,'ends':r['terminal_signature']}

    solvent_cache={}
    def solvent(s):
        if s not in solvent_cache:
            m=Chem.MolFromSmiles(s) if s else None
            solvent_cache[s]=Chem.MolToSmiles(m,isomericSmiles=True) if m and len(Chem.GetMolFrags(m))==1 else ''
        return solvent_cache[s]

    measurements=defaultdict(list)
    for r in load(args.data):
        s=r.get('canonical_smiles','')
        if s not in anns or (r.get('source_tag_name') or '').casefold()!='cyanine':continue
        y=finite(r.get('absorption_nm'))
        sol=solvent((r.get('solvent') or '').strip())
        if y is None or not sol:continue
        measurements[(s,sol,doi(r.get('reference_doi')))].append(y)

    bymolsol=defaultdict(list)
    for (s,sol,p),ys in measurements.items():
        bymolsol[(s,sol)].append((p,float(median(ys))))

    structs_by_n_sol=defaultdict(set)
    for s,sol in bymolsol:
        structs_by_n_sol[(anns[s]['n'],sol)].add(s)

    stats=defaultdict(Counter)
    closest=[]
    for (n,sol),structures in sorted(structs_by_n_sol.items()):
        for s1,s2 in combinations(sorted(structures),2):
            one,two=anns[s1],anns[s2]
            different=[i for i,(x,y) in enumerate(zip(one['R'],two['R'])) if x!=y]
            same_end=one['ends']==two['ends']
            one_site=len(different)==1
            dxi=two['xi']-one['xi']
            nonzero_xi=abs(dxi)>0.045
            doi_a={p for p,_ in bymolsol[(s1,sol)] if p}
            doi_b={p for p,_ in bymolsol[(s2,sol)] if p}
            common_doi=doi_a&doi_b
            c=stats[n]
            c['same_length_and_solvent_pairs']+=1
            if common_doi:c['same_length_solvent_DOI_pairs']+=1
            if same_end:c['plus_identical_terminal_fragments']+=1
            if one_site:c['single_chain_substitution_ignore_terminal']+=1
            if same_end and one_site:c['plus_one_changed_chain_substituent']+=1
            if same_end and one_site and nonzero_xi:c['plus_nonzero_xi_change']+=1
            if same_end and one_site and nonzero_xi and common_doi:
                c['strict_same_DOI_possible_pairs']+=1
            if one_site and nonzero_xi and len(closest)<1000:
                closest.append({
                    'chain':f'Cy-{n}','solvent':sol,'same_terminal_fragments':int(same_end),
                    'same_DOI_available':int(bool(common_doi)),
                    'shared_DOIs':';'.join(sorted(common_doi)),
                    'substitution_site':different[0]+1,
                    'substituent_A':one['R'][different[0]],
                    'substituent_B':two['R'][different[0]],
                    'delta_xi':round(dxi,4),'smiles_A':s1,'smiles_B':s2
                })

    result={}
    for n in (3,5,7,9):
        mols={s for (s,sol) in bymolsol if anns[s]['n']==n}
        ends=defaultdict(set)
        for s in mols:ends[anns[s]['ends']].add(s)
        result[f'Cy-{n}']={
            'xi_annotated_molecules_with_absorption':len(mols),
            'distinct_xi_values':len({anns[s]['xi'] for s in mols}),
            'distinct_terminal_signatures':len(ends),
            'max_molecules_same_terminals':max(map(len,ends.values()),default=0),
            'distinct_solvents':len({sol for (s,sol) in bymolsol if anns[s]['n']==n}),
            **{k:stats[n][k] for k in (
                'same_length_and_solvent_pairs',
                'same_length_solvent_DOI_pairs',
                'plus_identical_terminal_fragments',
                'single_chain_substitution_ignore_terminal',
                'plus_one_changed_chain_substituent',
                'plus_nonzero_xi_change',
                'strict_same_DOI_possible_pairs')}
        }
    report={
        'status':'BOTTLENECK_DIAGNOSTIC_ONLY_NOT_A_TEST',
        'note':'Pairs counted at unique-molecule and solvent level; a pair may appear in multiple solvents and is not independent.',
        'by_chain':result,
        'limitations':['Matching ignores isomeric measurement conditions within papers.',
                       'Single-chain-substitution pairs with DIFFERENT terminal fragments are NOT valid tests.',
                       'Approximate Hammett parameters from Stage16 are not independently validated.',
                       'Unrecognized cyanine architectures remain untested.',
                       'No confirmation or rejection of the 2026 law is inferred.']
    }
    (args.out/'stage16_bottleneck_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    cols=['chain','solvent','same_terminal_fragments','same_DOI_available','shared_DOIs','substitution_site','substituent_A','substituent_B','delta_xi','smiles_A','smiles_B']
    with (args.out/'stage16_near_miss_structures_NOT_TESTS.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=cols);writer.writeheader();writer.writerows(closest)
    print(json.dumps(report,ensure_ascii=False,indent=2));print('FINISHED')

if __name__=='__main__':main()
