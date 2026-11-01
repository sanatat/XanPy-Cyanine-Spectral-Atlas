#!/usr/bin/env python3
"""Stage 17: coverage audit of the 2026 cyanine xi rule in ALL FluoDB tags.

Only exact same terminal fragments, solvent and chain length are compared;
within-paper pairs are separate from across-paper exploratory pairs.
A successful spectral trend is *not* an independent validation until the
structures, compound identities, citations and chemical assumptions are reviewed.
"""
import argparse
import collections
import csv
import importlib.util
import itertools
import json
import math
from pathlib import Path
from statistics import median
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.error')

LENGTHS = (3, 5, 7, 9)


def load_stage16(path):
    if not path.exists():
        raise SystemExit(f'Missing stage16 code: {path}')
    spec = importlib.util.spec_from_file_location('xanpy_stage16_for_stage17', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def nfloat(value):
    try:
        f = float(value)
        return f if math.isfinite(f) and 200 <= f <= 1500 else None
    except (TypeError, ValueError):
        return None


def save_csv(path, rows, columns):
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=columns, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def main():
    home = Path.home() / 'Downloads'
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', default=str(home/'XanPy_QM_Atlas'/'stage2_all_endpoints.csv'))
    ap.add_argument('--out', default=str(home/'XanPy_QM_Atlas'))
    ap.add_argument('--stage16-script', default=str(home/'XanPy-Stage16-Direct-Xi-Test'/'stage16_direct_xi_test.py'))
    opt = ap.parse_args()
    inp, out = Path(opt.input), Path(opt.out)
    if not inp.exists():
        raise SystemExit(f'Data not found: {inp}')
    out.mkdir(parents=True, exist_ok=True)
    engine = load_stage16(Path(opt.stage16_script))
    total = collections.Counter()
    measurements = collections.defaultdict(list)
    source_tags = collections.defaultdict(set)
    sol_cache = {}

    with inp.open(newline='', encoding='utf-8-sig') as f:
        for r in csv.DictReader(f):
            total['rows'] += 1
            smi = (r.get('canonical_smiles') or '').strip()
            if not smi:
                total['no_smiles'] += 1
                continue
            v = nfloat(r.get('absorption_nm'))
            if v is None:
                total['no_valid_absorption'] += 1
                continue
            tag = (r.get('source_tag_name') or 'unclassified').strip()
            solraw = (r.get('solvent') or '').strip()
            if solraw not in sol_cache:
                mol = Chem.MolFromSmiles(solraw) if solraw else None
                sol_cache[solraw] = (Chem.MolToSmiles(mol,isomericSmiles=True)
                                     if mol and len(Chem.GetMolFrags(mol)) == 1 else '')
            sol = sol_cache[solraw]
            if not sol:
                total['invalid_solvent'] += 1
                continue
            source_tags[smi].add(tag)
            doi = engine.doi(r.get('reference_doi'))
            measurements[(smi, sol, doi)].append(v)
            total['usable_measurement_rows'] += 1

    annotations, reasons = {}, collections.Counter()
    for idx, smi in enumerate(sorted(source_tags)):
        info, why = engine.characterize(smi)
        reasons[why] += 1
        if info is not None:
            annotations[smi] = info

    # Only structurally recognized and explicitly xi-annotated compounds enter comparisons.
    med = {key: float(median(vals)) for key, vals in measurements.items()
           if key[0] in annotations and annotations[key[0]]['xi'] is not None}

    # A full, inspectable list of the structurally accepted compounds.
    ann_rows = []
    for smi, info in annotations.items():
        if info['xi'] is None:
            continue
        ann_rows.append({
            'smiles': smi, 'chain_length': info['n'],
            'xi': info['xi'], 'R': ' | '.join(info['R']),
            'terminal_signature': info['endgroup_signature'],
            'source_tags': ';'.join(sorted(source_tags[smi])),
            'with_absorption_measurements': sum(k[0]==smi for k in med),
        })
    save_csv(out/'stage17_alltags_xi_annotated.csv', ann_rows,
             ['smiles','chain_length','xi','R','terminal_signature','source_tags','with_absorption_measurements'])

    def build_pair(a, b, solvent, doi_a, doi_b, y_a, y_b, mode):
        if a == b:
            return None
        ai, bi = annotations[a], annotations[b]
        diffs = [j for j,(x,y) in enumerate(zip(ai['R'],bi['R'])) if x!=y]
        if len(diffs) != 1:
            return None
        delta_xi = bi['xi'] - ai['xi']
        if abs(delta_xi) < .045:
            return None
        observed_delta = y_b-y_a
        predicted_sign = 1 if ai['n'] in (3,7) else -1
        predicted_sign *= 1 if delta_xi > 0 else -1
        agrees = None if abs(observed_delta) < 1e-8 else int(predicted_sign * observed_delta > 0)
        j = diffs[0]
        return dict(chain_length=ai['n'], position=j+1, meso=int(j == ai['n']//2),
                    solvent=solvent, same_doi=int(mode=='same_doi'), doi_a=doi_a,
                    doi_b=doi_b, smiles_a=a, smiles_b=b,
                    source_tags_a=';'.join(sorted(source_tags[a])),
                    source_tags_b=';'.join(sorted(source_tags[b])),
                    substituent_a=ai['R'][j], substituent_b=bi['R'][j],
                    xi_a=ai['xi'], xi_b=bi['xi'], delta_xi=round(delta_xi,4),
                    absorption_a=round(y_a,3), absorption_b=round(y_b,3),
                    delta_abs_nm=round(observed_delta,3), direction_agrees=agrees,
                    terminal_signature=ai['endgroup_signature'])

    def pairs_from_groups(groups, mode):
        pairs=[]; counts=collections.Counter()
        for key, entries in groups.items():
            if len(entries) < 2:
                continue
            counts['same_length_terminals_solvent_groups_with_2plus']+=1
            for aa,bb in itertools.combinations(sorted(entries,key=lambda x:x[0]),2):
                counts['same_length_terminals_solvent_pairs']+=1
                a,y_a,p_a = aa; b,y_b,p_b=bb
                if a==b:continue
                if mode!='same_doi':
                    if not p_a or not p_b or p_a & p_b:
                        continue
                counts['pairs_pass_paper_condition']+=1
                diffs=sum(x!=y for x,y in zip(annotations[a]['R'], annotations[b]['R']))
                if diffs!=1:
                    continue
                counts['one_site_differing']+=1
                if abs(annotations[b]['xi']-annotations[a]['xi'])<.045:
                    continue
                counts['nonzero_xi_difference']+=1
                doi_a=key[-1] if mode=='same_doi' else '|'.join(sorted(p_a))
                doi_b=key[-1] if mode=='same_doi' else '|'.join(sorted(p_b))
                pair=build_pair(a,b,key[2],doi_a,doi_b,y_a,y_b,mode)
                if pair:
                    pairs.append(pair)
        return pairs, dict(counts)

    within_groups=collections.defaultdict(list)
    for (s,sol,d),y in med.items():
        if not d:continue
        meta=annotations[s]
        within_groups[(meta['n'],meta['endgroup_signature'],sol,d)].append((s,y,d))

    strict,strict_counts=pairs_from_groups(within_groups,'same_doi')

    aggregated=collections.defaultdict(list)
    for (s,sol,d),y in med.items():
        if d:aggregated[(s,sol)].append((d,y))

    # Across-publication comparisons are explorations, not controlled tests.
    cross_groups=collections.defaultdict(list)
    for (s,sol),vals in aggregated.items():
        meta=annotations[s]
        cross_groups[(meta['n'],meta['endgroup_signature'],sol)].append(
            (s,float(median(y for _,y in vals)),set(d for d,_ in vals)))
    exploratory,exploratory_counts=pairs_from_groups(cross_groups,'different_doi')

    pair_cols=['chain_length','position','meso','solvent','same_doi','doi_a','doi_b',
               'smiles_a','smiles_b','source_tags_a','source_tags_b','substituent_a',
               'substituent_b','xi_a','xi_b','delta_xi','absorption_a','absorption_b',
               'delta_abs_nm','direction_agrees','terminal_signature']
    save_csv(out/'stage17_same_DOI_pairs_FOR_MANUAL_REVIEW.csv', strict, pair_cols)
    save_csv(out/'stage17_across_DOI_pairs_EXPLORATORY.csv',exploratory,pair_cols)

    by_length={}
    for n in LENGTHS:
        a={s:i for s,i in annotations.items() if i['n']==n and i['xi'] is not None}
        used={s for s,_,_ in med if s in a}
        unique_xi=sorted({a[s]['xi'] for s in used})
        terms=collections.defaultdict(set)
        for s in used:terms[a[s]['endgroup_signature']].add(s)
        sp=[p for p in strict if p['chain_length']==n]
        ep=[p for p in exploratory if p['chain_length']==n]
        by_length[f'Cy-{n}']={
            'identified_molecules_with_absorption': len(used),
            'distinct_xi_values':len(unique_xi),
            'xi_values_examples':unique_xi[:30],
            'distinct_terminal_signatures':len(terms),
            'same_terminal_groups_2plus':sum(len(v)>1 for v in terms.values()),
            'nonzero_xi_molecules':sum(abs(a[s]['xi'])>.045 for s in used),
            'same_DOI_single_substitution_candidates':len(sp),
            'same_DOI_unique_molecular_pairs':len({tuple(sorted((p['smiles_a'],p['smiles_b']))) for p in sp}),
            'same_DOI_distinct_papers':len({p['doi_a'] for p in sp}),
            'across_DOI_exploratory_candidates':len(ep),
            'direction_agreement_same_DOI':(
                round(sum(p['direction_agrees'] for p in sp if p['direction_agrees'] is not None) /
                      sum(p['direction_agrees'] is not None for p in sp),4)
                if any(p['direction_agrees'] is not None for p in sp) else None),
        }
    report={
       'status':'ALL_TAGS_STRUCTURAL_AUDIT_NOT_FORMAL_VALIDATION',
       'input_rows':total['rows'],
       'unique_structures_with_valid_absorption_and_solvent':len(source_tags),
       'all_tags_structurally_accepted_molecules':len(annotations),
       'all_tags_xi_annotated_molecules':len(ann_rows),
       'rejection_reasons':dict(reasons),
       'by_length':by_length,
       'strict_pair_filter_counts':strict_counts,
       'across_DOI_filter_counts':exploratory_counts,
       'remarks':[
           'Source tags were NOT used to restrict the structural scan.',
           'Exact terminal-fragment identity is imposed to control auxochrome changes.',
           'Nonzero xi difference and a single chain substitution are required.',
           'Matched candidates require manual chemical identity and DOI verification.',
           'Published SI and FluoDB may share source papers, so overlap is not independent replication.',
           'No pair or no xi variance means this source cannot test the rule under these conditions.',
           'Xi is an approximate Hammett proxy for a limited set of simple substituents.'
       ]
    }
    (out/'stage17_alltags_audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))
    print('FINISHED')

if __name__=='__main__': main()
