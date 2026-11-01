import csv, json, subprocess, sys, tempfile
from pathlib import Path
from rdkit import Chem
sys.path.insert(0,str(Path(__file__).resolve().parent))
import stage16_direct_xi_test as m

base='OC(CCCCCN1/C(C(C)(C)C2=C1C=CC(S(O)(=O)=O)=C2)=C/C=C/C3=[N+](CC)C4=CC=C(S([O-])(=O)=O)C=C4C3(C)C)=O'
assert m.characterize(base)[0]['n']==3
with tempfile.TemporaryDirectory() as t:
    rows=[]
    for n in (3,5,7,9):
        smi=base if n==3 else base.replace('=C/C=C/C3', '=C/'+('/'.join(['C=C']*((n-1)//2)))+'/C3')
        mol=Chem.MolFromSmiles(smi)
        assert mol
        info,reason=m.characterize(smi);assert info and info['n']==n,(n,reason,info)
        (_,path,_),_=m.chain_candidates(mol)[0]
        c=path[n//2]
        rw=Chem.RWMol(mol); f=rw.AddAtom(Chem.Atom('F')); rw.AddBond(c,f,Chem.BondType.SINGLE)
        modified=rw.GetMol(); Chem.SanitizeMol(modified)
        sf=Chem.MolToSmiles(modified)
        b,reason=m.characterize(sf)
        assert b and b['xi'] is not None and b['endgroup_signature']==info['endgroup_signature'],(n,reason,b)
        sign=+1 if n in (3,7) else -1
        rows.extend([
            {'source_tag_name':'Cyanine','canonical_smiles':Chem.MolToSmiles(mol),'solvent':'CO','reference_doi':f'10.1234/{n}','absorption_nm':600,'emission_nm':620},
            {'source_tag_name':'Cyanine','canonical_smiles':sf,'solvent':'CO','reference_doi':f'10.1234/{n}','absorption_nm':600+sign*10,'emission_nm':620+sign*10},
        ])
    p=Path(t)/'raw.csv'
    with p.open('w',newline='') as f:
       w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    result=subprocess.run([sys.executable,str(Path(__file__).with_name('stage16_direct_xi_test.py')),'--input',str(p),'--output',t],capture_output=True,text=True)
    assert result.returncode==0, (result.stdout,result.stderr)
    d=json.loads((Path(t)/'stage16_direct_test_audit.json').read_text())
    print(json.dumps(d,indent=2))
    for n in (3,5,7,9):
       summary=d['strict_same_DOI_single_substitution_test'][f'Cy-{n}']
       assert summary['matched_comparisons']==1,summary
       assert summary['direction_agreement_fraction']==1.0,summary
