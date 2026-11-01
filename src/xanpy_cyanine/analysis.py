#!/usr/bin/env python3
"""XanPy stage 20: reproducible, descriptive review of original spectral traces.

No uncertainty intervals are inferred from single traces. Wavelength-domain widths are
apparent widths at fractions of exported peak height, NOT deconvolved electronic bands.
"""
from __future__ import annotations
import argparse, csv, hashlib, io, json, math, zipfile
from pathlib import Path
import numpy as np
from scipy.signal import savgol_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HASHES = {
    'Cy5-7':'15d35122836b39ef7e475c3904ffceaf',
    'Cy5-8':'de045dcd49ff0eb86ad80c65e98b382f',
    'Cy7-3':'5c199383545d8fb62cc6adff7011508c',
    'Cy7-8':'5cab22162f49e71a1275d3db5b877f4e'
}
PUBLISHED={'Cy5-7':(649.,693.),'Cy5-8':(649.,681.),'Cy7-3':(745.,806.),'Cy7-8':(741.,788.)}
PAIRS=[('Cy5-7','Cy5-8'),('Cy7-3','Cy7-8')]
COLORS={'Cy5-7':'#176d9c','Cy5-8':'#d65f31','Cy7-3':'#176d9c','Cy7-8':'#d65f31'}


def default_zip_dir()->Path:
    options=[Path(__file__).resolve().parents[2]/'data/raw', Path.home()/'Downloads/XanPy-Stage19-Independent-Validation/raw_zips', Path.home()/'Downloads']
    return next((p for p in options if all((p/f'{k}.zip').exists() for k in HASHES)),options[0])


def read_trace(z:zipfile.ZipFile, kind:str):
    names=[n for n in z.namelist() if ('absorption' if kind=='absorption' else 'emission') in n.lower() and n.lower().endswith('.txt') and ('spectr' in n.lower() or 'spetra' in n.lower())]
    if len(names)!=1:raise ValueError(f'Expected exactly one {kind} spectrum; got {names}')
    d=np.genfromtxt(io.BytesIO(z.read(names[0])),comments='#',invalid_raise=False)
    if d.ndim!=2 or d.shape[1]<2:raise ValueError('Data not 2-column')
    d=d[:, :2];d=d[np.isfinite(d).all(axis=1)];d=d[np.argsort(d[:,0])]
    x,y=d[:,0],d[:,1]
    if len(x)<40 or np.any(np.diff(x)<=0):raise ValueError('Invalid wavelength axis')
    return x,y,names[0]


def sg(y:np.ndarray, window:int):
    w=min(window,len(y)//2*2-1)
    return savgol_filter(y,w,3,mode='interp') if w>=5 else y.copy()


def parabolic_peak(x,y):
    k=int(np.argmax(y));
    if k==0 or k==len(y)-1:return float(x[k])
    # parabola around the three measured samples, vertex constrained locally
    a,b,c=np.polyfit(x[k-1:k+2]-x[k],y[k-1:k+2],2)
    offset= -b/(2*a) if a<0 else 0
    return float(x[k]+offset) if abs(offset)<1.5*min(x[k]-x[k-1],x[k+1]-x[k]) else float(x[k])


def width_at_fraction(x,y,peak_idx,fraction=.5):
    # Local crossings adjacent to main peak; a value is returned only if both wings cross.
    ylevel=float(y[peak_idx]*fraction)
    inds=np.flatnonzero(y[:peak_idx]<ylevel)
    jinds=np.flatnonzero(y[peak_idx+1:]<ylevel)
    if not len(inds) or not len(jinds):return None
    i=int(inds[-1]);j=int(peak_idx+1+jinds[0])
    def cross(k):
        x1,x2=x[k],x[k+1]; y1,y2=y[k],y[k+1]
        return float(x1+(ylevel-y1)*(x2-x1)/(y2-y1)) if y1!=y2 else float((x1+x2)/2)
    L=cross(i);R=cross(j-1)
    return {'width_nm':R-L,'left_nm':L,'right_nm':R,'left_halfwidth_nm':x[peak_idx]-L,'right_halfwidth_nm':R-x[peak_idx]}


def summarize(x,y,published,kind):
    # Actual maximum must be near the literature peak; suppress unrelated far-UV features.
    span=30 if kind=='absorption' else 45
    mask=(x>=published-span)&(x<=published+span)
    if mask.sum()<20: raise ValueError('not enough data around published peak')
    xx=x[mask];yy=y[mask]
    subgrid={}
    for method in ['raw','sg7','sg11','sg21','sg31','sg51']:
        vals=yy if method=='raw' else sg(yy,int(method[2:]))
        subgrid[method]=float(parabolic_peak(xx,vals))
    fullsm=sg(y,11)
    k=int(np.argmin(np.abs(x-subgrid['sg11'])))
    band={}
    width_smoothing={}
    for window in [7,11,21,31,51]:
        vv=sg(y,window)
        kk=int(np.argmin(np.abs(x-subgrid['sg'+str(window)])))
        crossing=width_at_fraction(x,vv,kk,.5)
        width_smoothing[str(window)]=round(crossing['width_nm'],3) if crossing else None
    for fraction in [.4,.5,.6,.7,.8]:
        v=width_at_fraction(x,fullsm,k,fraction)
        band[str(fraction)]=v
    return {'published_peak_nm':published,'raw_grid_peak_nm':float(xx[np.argmax(yy)]),
            'peak_estimates_nm':{k:round(v,4) for k,v in subgrid.items()},
            'peak_smoothing_range_nm':[round(min(subgrid.values()),4),round(max(subgrid.values()),4)],
            'wavelength_grid_step_nm':float(np.median(np.diff(x))),
            'n_points':int(len(x)), 'measured_wavelength_range_nm':[float(x.min()),float(x.max())],
            'width_at_fraction_of_exported_peak_height':band,
            'apparent_halfmax_width_by_savgol_window_nm':width_smoothing,
            'bandwidth_caution':'These are descriptive widths of overlapping exported bands, not baseline-corrected/deconvolved FWHM or measurement uncertainty.'}


def shape_comparison(xa,ya,xb,yb,peak_a,peak_b,kind):
    # Compare normalized curves on common wavelength domain; overlap after shifting to same peak.
    # For Cy5, absorption absorption near 650 has a broad blue wing; use a shared ±75 nm region.
    # For Cy7, use ±95 nm. Emission use asymmetric ±(35/75) nm.
    left,right=(85,85) if kind=='absorption' else (35,75)
    if kind=='absorption' and peak_a>700:left,right=95,85
    # shift each trace to a common peak-centred coordinate u
    u0=max(-left,xa.min()-peak_a,xb.min()-peak_b)
    u1=min(right,xa.max()-peak_a,xb.max()-peak_b)
    u=np.arange(math.ceil(u0),math.floor(u1)+.01,.5)
    A=np.interp(peak_a+u,xa,ya);B=np.interp(peak_b+u,xb,yb)
    # compare shapes (relative heights) after peak alignment
    A=A/max(np.max(A),1e-12); B=B/max(np.max(B),1e-12)
    mean_abs=float(np.mean(np.abs(A-B)))
    rms=float(np.sqrt(np.mean((A-B)**2)))
    overlap=float(np.trapezoid(np.minimum(A,B),u)/np.trapezoid(np.maximum(A,B),u))
    return {'peak_aligned_mean_absolute_normalized_signal_difference':mean_abs,
            'peak_aligned_RMS_normalized_signal_difference':rms,
            'peak_aligned_area_intersection_over_union':overlap,
            'relative_window_nm':[float(u[0]),float(u[-1])],
            'caution':'Descriptive shape metrics on normalized spectra; not uncertainty or quantum yield.'}


def run(zipdir,outdir):
    outdir.mkdir(parents=True,exist_ok=True)
    traces={};report={'status':'RAW_SPECTRA_VERIFIED_DESCRIPTIVE_NO_REPLICATES',
        'source_record':'https://zenodo.org/records/20119205',
        'source_doi':'10.1021/acs.joc.6c00573',
        'zip_md5_verified':{},'traces':{},'pairs':[],
        'limitations':['Single exported trace per condition: no experimental CI, significance or true wavelength calibration uncertainty.',
                       'Absorption sampling is 1 nm, so shifts of 1 nm or less are not robustly resolved.',
                       'The apparent width describes the exported band envelope, which may contain overlapping transitions.',
                       'Original signal appears normalized; intensity magnitude/area cannot imply brightness, quantum yield or extinction coefficient.',
                       'Sample pH/concentration/aggregation, repeat scans and instrument metadata were not independently verified.',
                       'The Hückel–Hammett rule is tested only for qualitative direction in selected pairs, not proven or rejected by these spectra.']}
    for compound,expected in HASHES.items():
        zp=zipdir/f'{compound}.zip'
        if not zp.exists():raise FileNotFoundError(zp)
        h=hashlib.md5(zp.read_bytes()).hexdigest()
        if h!=expected:raise ValueError(f'{zp}: MD5 mismatch expected {expected} got {h}')
        report['zip_md5_verified'][compound]=h
        with zipfile.ZipFile(zp) as z:
            report['traces'][compound]={}
            for kind,pub in zip(['absorption','emission'],PUBLISHED[compound]):
                x,y,name=read_trace(z,kind)
                traces[compound,kind]=(x,y)
                d=summarize(x,y,pub,kind);d['source_file_in_zip']=name
                report['traces'][compound][kind]=d
    for a,b in PAIRS:
        row={'pair':f'{a} -> {b}','data_source':'same published paper and solvent; not replicate scans'}
        for kind in ('absorption','emission'):
            da=report['traces'][a][kind];db=report['traces'][b][kind]
            row[kind]={
                'reported_table_difference_nm':db['published_peak_nm']-da['published_peak_nm'],
                'raw_grid_peak_difference_nm':db['raw_grid_peak_nm']-da['raw_grid_peak_nm'],
                'sg11_subgrid_difference_nm':round(db['peak_estimates_nm']['sg11']-da['peak_estimates_nm']['sg11'],4),
                'smooth_window_differences_nm':{w:round(db['peak_estimates_nm'][w]-da['peak_estimates_nm'][w],4) for w in ('raw','sg7','sg11','sg21','sg31','sg51')},
                'apparent_width_at_halfmax_nm':{
                    a:round(da['width_at_fraction_of_exported_peak_height']['0.5']['width_nm'],3) if da['width_at_fraction_of_exported_peak_height']['0.5'] else None,
                    b:round(db['width_at_fraction_of_exported_peak_height']['0.5']['width_nm'],3) if db['width_at_fraction_of_exported_peak_height']['0.5'] else None},
                'apparent_halfmax_width_difference_nm':round(db['width_at_fraction_of_exported_peak_height']['0.5']['width_nm']-da['width_at_fraction_of_exported_peak_height']['0.5']['width_nm'],3) if da['width_at_fraction_of_exported_peak_height']['0.5'] and db['width_at_fraction_of_exported_peak_height']['0.5'] else None,
                'shape_after_peak_alignment':shape_comparison(*traces[a,kind],*traces[b,kind],da['peak_estimates_nm']['sg11'],db['peak_estimates_nm']['sg11'],kind),
            }
        report['pairs'].append(row)
    (outdir/'stage20_spectral_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    rows=[]
    for c,d in report['traces'].items():
        for typ,r in d.items():
            w=r['width_at_fraction_of_exported_peak_height']['0.5']
            rows.append({'compound':c,'spectrum':typ,'published_peak_nm':r['published_peak_nm'],
                         'raw_grid_peak_nm':r['raw_grid_peak_nm'],
                         'sg11_subgrid_peak_nm':r['peak_estimates_nm']['sg11'],
                         'peak_method_min_nm':r['peak_smoothing_range_nm'][0],
                         'peak_method_max_nm':r['peak_smoothing_range_nm'][1],
                         'apparent_halfmax_bandwidth_nm':round(w['width_nm'],3) if w else '',
                         'apparent_blue_halfwidth_nm':round(w['left_halfwidth_nm'],3) if w else '',
                         'apparent_red_halfwidth_nm':round(w['right_halfwidth_nm'],3) if w else '',
                         'source_zip':c+'.zip','md5':report['zip_md5_verified'][c]})
    with (outdir/'spectral_features.csv').open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
    # Independent graphs: one plot per graph, no subplots.
    for kind in ['absorption','emission']:
        fig,ax=plt.subplots(figsize=(10.2,5.6))
        for c in ('Cy5-7','Cy5-8','Cy7-3','Cy7-8'):
            x,y=traces[c,kind]
            domain=(555,825) if kind=='absorption' else (625,910)
            mask=(x>=domain[0])&(x<=domain[1])
            ax.plot(x[mask],y[mask],label=c,lw=1.55,color=COLORS[c],ls='-' if c[-1] in ('7','3') else '--')
        ax.set_xlim(domain);ax.set_ylim(bottom=-.02)
        ax.set_ylabel('Original exported signal (approx. normalized)');ax.set_xlabel('Wavelength (nm)')
        ax.set_title(('Absorption' if kind=='absorption' else 'Emission')+' — archived instrument traces')
        ax.grid(alpha=.2);ax.legend(ncol=4);fig.tight_layout()
        fig.savefig(outdir/f'{kind}_four_raw_overlays.png',dpi=170);plt.close(fig)
    for a,b in PAIRS:
        for kind in ('absorption','emission'):
            fig,ax=plt.subplots(figsize=(9.3,5.3))
            for c in (a,b):
                x,y=traces[c,kind]
                peak=report['traces'][c][kind]['peak_estimates_nm']['sg11']
                r=100 if kind=='absorption' else 85
                mask=(x>=peak-r)&(x<=peak+r)
                ax.plot(x[mask],y[mask],label=f'{c} (raw)',lw=1.3,color=COLORS[c]);
                ax.axvline(peak,ls=':',lw=1.2,color=COLORS[c],alpha=.8)
            ax.set_xlabel('Wavelength (nm)');ax.set_ylabel('Exported instrument signal (approx. normalized)')
            ax.set_title(f'{a} vs {b}: {kind} (dotted = SG11 sub-grid maximum)')
            ax.grid(alpha=.2);ax.legend();fig.tight_layout()
            fig.savefig(outdir/f'{a}_vs_{b}_{kind}_detail.png',dpi=180);plt.close(fig)
    print(json.dumps({'zip_md5_verified':list(report['zip_md5_verified']),
                      'features':rows,'pairs':report['pairs']},ensure_ascii=False,indent=2))
    return report

if __name__=='__main__':
    parser=argparse.ArgumentParser(description='Reanalyse source-verified 2026 cyanine absorption/emission traces')
    parser.add_argument('--zip-dir',type=Path,default=default_zip_dir())
    parser.add_argument('--out-dir',type=Path,default=Path(__file__).resolve().parents[2]/'results')
    args=parser.parse_args()
    run(args.zip_dir,args.out_dir)
