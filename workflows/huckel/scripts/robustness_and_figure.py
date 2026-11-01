from pathlib import Path
import sys,csv
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from huckel_test import rows_for_chains, save_csv, CHAINS

ROOT=Path(__file__).resolve().parents[1]
rows=[]
for s in [0.01,0.05,0.10,0.15,0.25,0.40]:
    for r in rows_for_chains(s):
        rows.append({'strength_beta':s, 'chain_length':r['chain_length'],'position':r['position'],
             'electronic_class':r['substituent_type'],'predicted_direction':r['predicted_absorption_shift'],
             'delta_gap':r['exact_model_delta_gap']})
save_csv(ROOT/'results'/'sensitivity_96_cases.csv',rows)
reference={(r['chain_length'],r['position'],r['substituent_type']):r['predicted_absorption_shift'] for r in rows_for_chains(0.15)}
assert all(reference[(r['chain_length'],r['position'],r['electronic_class'])]==r['predicted_direction'] for r in rows)
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(8,4.5))
    for site in ('meso','alpha'):
        data=[next(r for r in rows_for_chains(0.15) if r['chain_length']==n and r['position']==site and r['substituent_type']=='EDG')['exact_model_delta_gap'] for n in CHAINS]
        ax.plot(CHAINS,data,marker='o',label=f'EDG, {site} position')
    ax.axhline(0,linewidth=1,linestyle='--')
    ax.set_xticks(CHAINS)
    ax.set_xlabel('Polymethine chain length (number of carbon sites)')
    ax.set_ylabel(r'Change in HOMO–LUMO gap ($|\beta|$ units)')
    ax.set_title('Numerical Hückel test: substituent-position sign inversion')
    ax.legend()
    fig.tight_layout()
    fig.savefig(ROOT/'figures'/'huckel_gap_sign_inversion.png',dpi=175)
    plt.close(fig)
except ImportError:
    print('Matplotlib not available; figure skipped')
print('Sensitivity test: 96 / 96 signs consistent across strengths')
