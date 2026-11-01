import json,os,subprocess,sys,csv
from pathlib import Path

def test_stage13_small_fixture(tmp_path):
    src=tmp_path/'stage12_candidate_records.csv'
    with src.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['chemical_state','smiles','doi','solvent','absorption_nm','emission_nm','source_row','source_tag','source_database'])
        w.writeheader()
        for i,(a,e,s) in enumerate([(450,510,'CO'),(452,530,'CC#N'),(455,520,'CCO')]):
            w.writerow(dict(chemical_state='O_bridge__conjugated_9C',smiles='COc1ccccc1',doi='10.1234/syntheticfixture',solvent=s,absorption_nm=a,emission_nm=e,source_row=i,source_tag='synthetic',source_database='SYNTHETIC'))
    repo=Path(__file__).resolve().parents[1]
    p=subprocess.run([sys.executable,str(repo/'scripts/stage13.py')],env={**os.environ,'XANPY_ROOT':str(tmp_path)},capture_output=True,text=True)
    assert p.returncode==0,p.stderr+p.stdout
    summary=json.loads((tmp_path/'stage13_o_conjugated_audit.json').read_text())
    assert summary['doi_traceable_molecule_paper_groups_with_at_least_2_solvents']==1
    assert summary['within_study_solvent_pairs_not_independent']==3

def test_user_reported_summary_marked(tmp_path):
    repo=Path(__file__).resolve().parents[1]
    x=json.loads((repo/'results/user_reported/summary.json').read_text())
    assert x['stage13']['molecules_in_multisolvent']==37
    assert x['stage13']['ci95_nm'][0] < .5 < x['stage13']['ci95_nm'][1]
