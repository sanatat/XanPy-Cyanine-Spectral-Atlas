import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from huckel_test import rows_for_chains,compare_examples,gap_shift,orbital_data

def test_four_lengths_sixteen_directions():
    rows=rows_for_chains()
    assert len(rows)==16
    lookup={(r['chain_length'],r['position'],r['substituent_type']):r['predicted_absorption_shift'] for r in rows}
    expected={
        3:('blue','red','red','blue'),
        5:('red','blue','blue','red'),
        7:('blue','red','red','blue'),
        9:('red','blue','blue','red')}
    for n, directions in expected.items():
        got=tuple(lookup[n,p,t] for p,t in [('meso','EDG'),('alpha','EDG'),('meso','EWG'),('alpha','EWG')])
        assert got==directions

def test_first_order_and_exact_agree_on_direction():
    for row in rows_for_chains():
        assert row['exact_model_delta_gap']*row['first_order_delta_gap']>0

def test_molecular_symmetry():
    for n in (3,5,7,9):
        center=n//2
        left=gap_shift(n,center-1,0.15)
        right=gap_shift(n,center+1,0.15)
        assert abs(left-right)<1e-12

def test_four_paper_examples_match_signs():
    assert all(x['qualitative_match'] for x in compare_examples(rows_for_chains()))
