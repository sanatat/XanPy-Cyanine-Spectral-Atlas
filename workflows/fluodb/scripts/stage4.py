#!/usr/bin/env python3
"""Paired model-gain bootstrap and within-molecule solvent-shift diagnostics.
Bootstraps evaluate the fixed test split only, not cross-split uncertainty.
"""
import os, csv, json
from collections import defaultdict
from pathlib import Path
import numpy as np
ROOT=Path(os.environ.get('XANPY_ROOT',str(Path.home()/'Downloads/XanPy_QM_Atlas')))
results=[]
for target in ('emission_nm','absorption_nm'):
    src=ROOT/f'stage3_{target}_predictions.csv'
    if not src.is_file():raise SystemExit('Missing '+str(src))
    observations=defaultdict(dict)
    with src.open(newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f):
            if r['model'] not in ('structure','structure_solvent'):continue
            solvent=r.get('solvent',r.get('solvent_raw_smiles',''))
            key=(r['split'],r['smiles'],solvent)
            d=observations[key]
            d['actual']=float(r['actual_nm']);d['scaffold']=r['scaffold']
            d[r['model']]=float(r['predicted_nm'])
    for split in ('molecule','scaffold'):
        clusters=defaultdict(lambda:[0.0,0]);by_mol=defaultdict(list);gains=[]
        for (sp,smi,sol),d in observations.items():
            if sp!=split or not ('structure' in d and 'structure_solvent' in d):continue
            gain=abs(d['actual']-d['structure'])-abs(d['actual']-d['structure_solvent'])
            group=d['scaffold'] if split=='scaffold' else smi
            clusters[group][0]+=gain;clusters[group][1]+=1
            gains.append(gain)
            if sol:by_mol[smi].append((sol,d))
        if not gains:continue
        rng=np.random.default_rng(20261101)
        values=np.asarray(list(clusters.values()),dtype=float)
        boot=np.empty(1000)
        for i in range(len(boot)):
            selected=values[rng.integers(0,len(values),len(values))]
            boot[i]=selected[:,0].sum()/selected[:,1].sum()
        within_base=[];within_solv=[];ranges=[]
        for smi,obs in by_mol.items():
            if len(set(s for s,_ in obs))<2:continue
            y=np.array([d['actual'] for _,d in obs]);p=np.array([d['structure'] for _,d in obs]);q=np.array([d['structure_solvent'] for _,d in obs])
            within_base.append(float(np.mean(np.abs((y-y.mean())-(p-p.mean())))))
            within_solv.append(float(np.mean(np.abs((y-y.mean())-(q-q.mean())))))
            ranges.append(float(y.max()-y.min()))
        result={'target':target,'split':split,'test_pairs':len(gains),'bootstrap_groups':len(values),
          'mae_gain_nm':round(float(np.mean(gains)),4),
          'gain_ci_low_nm':round(float(np.quantile(boot,.025)),4),
          'gain_ci_high_nm':round(float(np.quantile(boot,.975)),4),
          'multi_solvent_test_molecules':len(within_base),
          'centered_mae_structure_nm':round(float(np.mean(within_base)),4) if within_base else None,
          'centered_mae_with_solvent_nm':round(float(np.mean(within_solv)),4) if within_solv else None,
          'median_observed_solvent_range_nm':round(float(np.median(ranges)),4) if ranges else None}
        results.append(result)
(ROOT/'stage4_paired_audit.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
print(json.dumps(results,indent=2))
