#!/usr/bin/env python3
"""Minimal independent Hückel calculation for 2026 Dunlop--Strada--Sebej ξ-scheme.

This is a tight-binding orbital model, NOT DFT nor a molecular excitation-energy model.
On-site EDG/EWG perturbation is a modeling convention. Energies in units |beta|.
"""
from __future__ import annotations
import argparse
import csv
from pathlib import Path
import numpy as np

CHAINS = (3, 5, 7, 9)
SITE_PERTURBATION = 0.15

def hamiltonian(n: int) -> np.ndarray:
    if n < 3 or n % 2 != 1:
        raise ValueError('n must be an odd integer >= 3')
    h = np.zeros((n, n), dtype=float)
    for i in range(n - 1):
        h[i, i+1] = -1.0
        h[i+1, i] = -1.0
    return h

def orbital_data(n: int):
    e, vectors = np.linalg.eigh(hamiltonian(n))
    # (n+1) electrons: for n odd, occupy (n+1)/2 orbitals doubly.
    homo = (n - 1)//2
    lumo = homo+1
    return e, vectors, homo, lumo

def gap_shift(n: int, position: int, onsite_shift: float) -> float:
    e, _, h, l = orbital_data(n)
    initial_gap = e[l] - e[h]
    hp = hamiltonian(n)
    hp[position, position] += onsite_shift
    ep = np.linalg.eigvalsh(hp)
    return float(ep[l]-ep[h]-initial_gap)

def rows_for_chains(strength: float = SITE_PERTURBATION):
    rows = []
    for n in CHAINS:
        es, vec, h, l = orbital_data(n)
        for name, idx in (('meso', n//2), ('alpha', n//2 - 1)):
            dh = float(vec[idx, h]**2)
            dl = float(vec[idx, l]**2)
            derivative = dl - dh # first-order d(gap)/d(onsite energy)
            for name_electronic, sign in (('EDG', 1), ('EWG', -1)):
                onsite = sign * strength
                delta = gap_shift(n, idx, onsite)
                # The sign of a gap shift predicts the qualitative direction of absorption.
                direction = 'blue' if delta > 0 else 'red' if delta < 0 else 'no_shift'
                rows.append({
                    'chain': f'Cy-{n}', 'chain_length': n,
                    'class': '4n-1' if n%4 == 3 else '4n+1',
                    'position': name, 'site_index_1based': idx+1,
                    'substituent_type': name_electronic,
                    'onsite_shift_beta_units': round(onsite, 5),
                    'homo_density_on_site': round(dh, 7),
                    'lumo_density_on_site': round(dl, 7),
                    'first_order_gap_derivative': round(derivative, 7),
                    'first_order_delta_gap': round(derivative*onsite, 7),
                    'exact_model_delta_gap': round(delta, 7),
                    'predicted_absorption_shift': direction,
                })
    return rows

PAPER_EXAMPLES = [
    # Approximate shifts explicitly stated in the 2026 paper, *not* recomputed from SI.
    {'chain': 'Cy-5', 'position': 'meso', 'substituent': 'OH', 'type': 'EDG', 'reported_absorption_shift_nm_approx': +30},
    {'chain': 'Cy-5', 'position': 'meso', 'substituent': 'NO2', 'type': 'EWG', 'reported_absorption_shift_nm_approx': -35},
    {'chain': 'Cy-7', 'position': 'meso', 'substituent': 'OEt', 'type': 'EDG', 'reported_absorption_shift_nm_approx': -45},
    {'chain': 'Cy-7', 'position': 'meso', 'substituent': 'CN', 'type': 'EWG', 'reported_absorption_shift_nm_approx': +65},
]

def compare_examples(calculated):
    key_to_row = {(r['chain'],r['position'],r['substituent_type']):r for r in calculated}
    out=[]
    for item in PAPER_EXAMPLES:
        model = key_to_row[(item['chain'],item['position'],item['type'])]
        paper_direction = 'red' if item['reported_absorption_shift_nm_approx'] > 0 else 'blue'
        out.append({**item,
            'paper_direction': paper_direction,
            'huckel_direction': model['predicted_absorption_shift'],
            'qualitative_match': paper_direction==model['predicted_absorption_shift'],
            'note':'Reported approximate contrast from 2026 paper, not an independently measured datapoint'})
    return out

def save_csv(path:Path, rows:list):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)

def main():
    pa=argparse.ArgumentParser()
    pa.add_argument('--out',default=None,help='Output directory (default: project results)')
    pa.add_argument('--strength',type=float,default=SITE_PERTURBATION)
    args=pa.parse_args()
    if args.strength<=0 or args.strength > 0.4:pa.error('use 0<strength<=0.4 for a small perturbation')
    root=Path(__file__).resolve().parent.parent
    out=Path(args.out) if args.out else root/'results'
    out.mkdir(parents=True,exist_ok=True)
    result=rows_for_chains(args.strength)
    compare=compare_examples(result)
    save_csv(out/'huckel_16_sign_tests.csv',result)
    save_csv(out/'paper_four_examples_sign_check.csv',compare)
    print('HÜCKEL MODEL: 4 lengths x 2 positions x 2 substituent polarities = 16 tests')
    print('Chain  meso-EDG  alpha-EDG  meso-EWG  alpha-EWG')
    lookup={(r['chain_length'],r['position'],r['substituent_type']):r['predicted_absorption_shift'] for r in result}
    for n in CHAINS:
        print(f'{n:>5}  {lookup[n,"meso","EDG"]:>8}  {lookup[n,"alpha","EDG"]:>9}  {lookup[n,"meso","EWG"]:>8}  {lookup[n,"alpha","EWG"]:>9}')
    print(f'Paper excerpt sign matches: {sum(x["qualitative_match"] for x in compare)}/4 (NOT full experimental validation)')
    print('Saved:',out)

if __name__=='__main__':main()
