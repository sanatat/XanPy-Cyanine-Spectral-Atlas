"""Synthetic numerical unit tests and archived-data integration if available."""
import json
from pathlib import Path
import numpy as np
import pytest
from xanpy_cyanine.analysis import (parabolic_peak, width_at_fraction, sg,
                                   shape_comparison, summarize, run, default_zip_dir)

def test_quadratic_gaussian_peak():
    x=np.linspace(-20,20,1601)
    y=np.exp(-.5*((x-1.2)/3)**2)
    assert abs(parabolic_peak(x,y)-1.2)<.01

def test_fwhm_gaussian():
    x=np.linspace(-20,20,1601)
    y=sg(np.exp(-.5*((x-1.2)/3)**2),11)
    got=width_at_fraction(x,y,int(np.argmax(y)),.5)
    expected=2*np.sqrt(2*np.log(2))*3
    assert got is not None
    assert abs(got['width_nm']-expected)<.02

def test_missing_wing_no_width():
    x=np.arange(1,31,dtype=float)
    y=x.copy()
    assert width_at_fraction(x,y,len(y)-1,.5) is None

def test_subgrid_shift():
    x=np.linspace(-20,20,1601)
    y=np.exp(-.5*((x-1.2)/3)**2)
    z=np.exp(-.5*((x-.2)/3)**2)
    assert abs(parabolic_peak(x,z)-parabolic_peak(x,y)+1)<.01

def test_shape_identical():
    x=np.linspace(600,750,751)
    y=np.exp(-.5*((x-680)/22)**2)
    d=shape_comparison(x,y,x,y,680,680,'absorption')
    assert d['peak_aligned_mean_absolute_normalized_signal_difference']<1e-12
    assert d['peak_aligned_area_intersection_over_union']==pytest.approx(1)

def test_summaries_of_synthetic_trace():
    x=np.linspace(550,750,1001)
    y=np.exp(-.5*((x-650)/15)**2)
    result=summarize(x,y,650,'absorption')
    assert result['raw_grid_peak_nm']==pytest.approx(650)
    assert result['width_at_fraction_of_exported_peak_height']['0.5']['width_nm']==pytest.approx(35.32,abs=.03)

@pytest.mark.integration
def test_against_original_archives_when_present(tmp_path):
    source=default_zip_dir()
    if not all((source/f'{k}.zip').exists() for k in ['Cy5-7','Cy5-8','Cy7-3','Cy7-8']):
        pytest.skip('Original archived 4 ZIPs not present: numerical tests only')
    result=run(source,tmp_path)
    pairs={p['pair']:p for p in result['pairs']}
    cy5=pairs['Cy5-7 -> Cy5-8']['absorption']
    cy7=pairs['Cy7-3 -> Cy7-8']['absorption']
    assert cy5['sg11_subgrid_difference_nm']==pytest.approx(-.3824,abs=.03)
    assert cy5['apparent_halfmax_width_difference_nm']==pytest.approx(-33.448,abs=.05)
    assert cy7['sg11_subgrid_difference_nm']==pytest.approx(-4.6068,abs=.03)
