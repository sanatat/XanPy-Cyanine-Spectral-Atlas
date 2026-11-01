#!/usr/bin/env python3
"""Validate experimentally curated ξ--absorption trend from 2026 SI.

CSV columns: molecule_id,chain_length,solvent,doi,xi,sigma_source,absorption_nm
xi must be calculated from the paper's positional weighting, not fabricated.
A chain has a regression only if >= 3 distinct ξ values, >= 4 data points.
This tool does not auto-label cyanines or extract SIs from PDFs/DOCX.
"""
import argparse,csv,math,json
from collections import defaultdict
from pathlib import Path
import numpy as np

def main():
    p=argparse.ArgumentParser()
    p.add_argument('csv_path',help='Manually curated actual SI measurements')
    p.add_argument('--out',default='hammett_results.json')
    args=p.parse_args()
    data=defaultdict(list)
    with open(args.csv_path,newline='',encoding='utf-8-sig') as f:
        r=csv.DictReader(f)
        required={'molecule_id','chain_length','solvent','doi','xi','sigma_source','absorption_nm'}
        if not required.issubset(r.fieldnames or []):raise SystemExit(f'Missing: {required-set(r.fieldnames or [])}')
        for a in r:
            if not a['doi'].strip() or not a['sigma_source'].strip() or not a['molecule_id'].strip():
                raise SystemExit('Each input needs molecule_id, DOI and sigma_source')
            n=int(a['chain_length'])
            x=float(a['xi']); y=float(a['absorption_nm'])
            if not(n in (3,5,7,9) and math.isfinite(x) and math.isfinite(y) and 200<=y<=1500):
                raise SystemExit(f'Bad input {a}')
            data[(n,a['solvent'])].append((x,y,a['molecule_id'],a['doi']))
    if not data:
        raise SystemExit("No real SI data in the input template; collect values from official Supporting Information first.")
    report={}
    for (n,sol),entries in sorted(data.items()):
        x=np.array([e[0] for e in entries]);y=np.array([e[1] for e in entries]);ids=[e[2] for e in entries]
        if len(x)<4 or len(set(x))<3:
            report[f'Cy-{n}|{sol}']={'status':'insufficient','n':len(x)};continue
        if len(set(ids))<4:
            report[f'Cy-{n}|{sol}']={'status':'insufficient_unique_molecules','n':len(x)};continue
        slope,intercept=np.polyfit(x,y,1)
        corr=float(np.corrcoef(x,y)[0,1]);expect=1 if n%4==3 else -1
        report[f'Cy-{n}|{sol}']={'status':'evaluated','n':len(x),'unique_molecules':len(set(ids)), 'slope_nm_per_xi':round(float(slope),6),'pearson_r':round(corr,6),'expected_slope_sign':('positive' if expect==1 else 'negative'),'sign_matches_paper_rule':bool(slope*expect>0),'note':'Descriptive fit; samples may be dependent and solvation/auxochromes may confound'}
    Path(args.out).write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
