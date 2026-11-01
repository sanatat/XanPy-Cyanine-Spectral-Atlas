import csv
import json
import os
import subprocess
import sys
from pathlib import Path

D = Path.home() / 'Downloads'
P = Path(os.environ.get('XANPY_ROOT', str(D / 'XanPy_QM_Atlas')))
source = Path(__file__).resolve().parent / 'stage3.py'
target = sys.argv[1] if len(sys.argv) > 1 else 'emission_nm'
if target not in ('emission_nm', 'absorption_nm'):
    sys.exit('Use emission_nm or absorption_nm')
if not source.is_file():
    sys.exit('Missing original script: ' + str(source))

original_metrics = P / f'stage3_{target}_metrics.csv'
old = {}
if original_metrics.is_file():
    with original_metrics.open(newline='') as f:
        for r in csv.DictReader(f):
            if r.get('split') == 'scaffold' and r.get('model') in ('structure', 'structure_solvent'):
                old[r['model']] = float(r['MAE_nm'])

code = source.read_text(encoding='utf-8')
# Keep the same architecture, features, and RF random_state as Stage 3.
# Change ONLY the grouped test split seed and remove molecule split to reduce work.
groups = ("('molecule', mol_groups),", "('molecule', smiles_arr),")
if not any(x in code for x in groups):
    sys.exit('Could not locate molecule split in Stage 3; not modifying source.')
for x in groups:
    code = code.replace(x, '')
patterns = (
    ('test_size=0.20,\n        random_state=42',
     'test_size=0.20,\n        random_state=int(__import__("os").environ["XANPY_SEED"])'),
    ('test_size=0.20, random_state=42',
     'test_size=0.20, random_state=int(__import__("os").environ["XANPY_SEED"])')
)
changed = 0
for a, b in patterns:
    if a in code:
        code = code.replace(a, b, 1)
        changed += 1
if changed != 1:
    sys.exit('Could not safely locate test split random_state=42')

# Do not overwrite your Stage 3 metrics, predictions, or audit.
code = code.replace('f\'stage3_{TARGET}_metrics.csv\'', 'f\'stage5_tmp_{TARGET}_metrics.csv\'')
code = code.replace('f\'stage3_{TARGET}_predictions.csv\'', 'f\'stage5_tmp_{TARGET}_predictions.csv\'')
code = code.replace('f\'stage3_{TARGET}_audit.json\'', 'f\'stage5_tmp_{TARGET}_audit.json\'')
worker = P / '_xanpy_stage5_worker.py'
worker.write_text(code, encoding='utf-8')

rows = []
for seed in (11, 42, 2026):
    env = dict(os.environ, XANPY_SEED=str(seed))
    print(f'\n=== Scaffold holdout seed {seed} ===', flush=True)
    subprocess.run([sys.executable, str(worker), target], env=env, check=True)
    metrics_path = P / f'stage5_tmp_{target}_metrics.csv'
    with metrics_path.open(newline='', encoding='utf-8') as f:
        m = {r['model']: r for r in csv.DictReader(f) if r['split'] == 'scaffold'}
    if not {'structure', 'structure_solvent'}.issubset(m):
        sys.exit('Expected structure and structure_solvent result missing')
    a, b = float(m['structure']['MAE_nm']), float(m['structure_solvent']['MAE_nm'])
    row = dict(seed=seed, target=target, test_type='scaffold',
               test_pairs=m['structure'].get('test_pairs', m['structure'].get('n_test')),
               structure_MAE_nm=a, structure_solvent_MAE_nm=b,
               improvement_nm=round(a-b, 4))
    rows.append(row)
    print('Improvement:', round(a-b, 4), 'nm', flush=True)

output = P / f'stage5_{target}_repeated_scaffold.csv'
with output.open('w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)

checks = []
if old:
    checks = [{"model": label, "seed42_minus_original_MAE_nm":
               round(next(r for r in rows if r['seed'] == 42)[
                   'structure_MAE_nm' if label == 'structure' else 'structure_solvent_MAE_nm'] - value, 5)}
              for label, value in old.items()]
report = dict(target=target, n_splits=len(rows),
              mean_improvement_nm=round(sum(r['improvement_nm'] for r in rows)/len(rows), 4),
              positive_splits=sum(r['improvement_nm'] > 0 for r in rows),
              seed42_check=checks, per_split=rows,
              caveat='Exploratory repeated scaffold holdouts; only 3 splits, not definitive statistical inference.')
report_path = P / f'stage5_{target}_summary.json'
report_path.write_text(json.dumps(report, indent=2), encoding='utf-8')
print('\n',json.dumps(report, indent=2), flush=True)
print('FINISHED:', report_path)
