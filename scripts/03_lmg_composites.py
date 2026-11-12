#!/usr/bin/env python3
"""v7 step 7: significant-only LMG attribution (Fig 3A data) + composite construct validation.
Attribution: per index, OLS y ~ significant category PC1s (FDR<0.10, MEM10 primary screen)
+ MEM1-10 spatial block; LMG (exact, all orderings) decomposes R² into blocks.
Indices with no significant category: MEM-only model (bar all grey, marked n.s.).
Composites: spatial partial Spearman (rank, MEM1-10 residualized) for
USAI~xenobiotic_degrad, IPI~EMHI, EMHI~pathogen_relab.
Outputs: stats_ms/attribution_v7.csv, composites_validation_v7.csv
"""
import os, warnings, itertools; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd, statsmodels.api as sm
from scipy import stats as sst


from pathlib import Path
import os
PACKAGE = Path(__file__).resolve().parents[1]
BASE = str(PACKAGE / 'data')
MS = str(Path(os.environ['URBAN_RESULTS_DIR']).resolve())
IDX = str(PACKAGE / 'data' / 'indices')
SW = str(PACKAGE / 'data' / 'sewage') + '/'

def wcsv(df, name, **kw):
    p = f'{MS}/{name}'
    if os.path.exists(p): os.remove(p)
    df.to_csv(p, **kw)
    print(f'  wrote {name} {df.shape}', flush=True)

def haversine(lat1, lon1, lat2, lon2):
    R = 6371.0; p1, p2 = np.radians(lat1), np.radians(lat2)
    dp, dl = np.radians(lat2-lat1), np.radians(lon2-lon1)
    a = np.sin(dp/2)**2 + np.cos(p1)*np.cos(p2)*np.sin(dl/2)**2
    return 2*R*np.arcsin(np.sqrt(a))

def dbmem(coords, n_axes=10):
    n = len(coords); D = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            D[i, j] = D[j, i] = haversine(coords[i,0], coords[i,1], coords[j,0], coords[j,1])
    A = -0.5*D**2; Bm = A - A.mean(0) - A.mean(1)[:, None] + A.mean()
    ev, ec = np.linalg.eigh(Bm); o = np.argsort(ev)[::-1]
    ev, ec = ev[o], ec[:, o]
    return ec[:, ev > 0][:, :n_axes]

city = pd.read_csv(f'{IDX}/indices_city_v7.csv').set_index('city')
ms6 = pd.read_csv(f'{MS}/indices_city_ms.csv'); ms6['city'] = ms6['city'].astype(str).str.strip(); ms6 = ms6.set_index('city')
for extra in ['func_redundancy', 'func_tax_decouple']:
    city[extra] = ms6[extra].reindex(city.index)
co = pd.read_csv(f'{BASE}/city_coords_metasub.csv').set_index('city')
common = [c for c in city.index if c in co.index]
idxc = city.loc[common]
MEM = pd.DataFrame(dbmem(co.loc[common, ['lat','lon']].values, 10), index=common,
                   columns=[f'MEM{i+1}' for i in range(10)])
PC1 = pd.read_csv(f'{MS}/category_pc1_v7.csv', index_col=0)
cs = pd.read_csv(f'{MS}/category_screen_v7.csv')  # MEM10 primary (city-level)
_cts = pd.read_csv(f'{MS}/category_twostage_v7.csv')
_confirmed = set(zip(_cts.loc[(_cts.fdr_sample<0.10)&(_cts.direction_match),'index'],
                     _cts.loc[(_cts.fdr_sample<0.10)&(_cts.direction_match),'category']))

FUNC16 = ['pathogen_relab','pathway_richness','pathway_evenness','amr_total_depth_perM',
          'amr_class_diversity','vf_total_depth_perM','VPI','allergen_richness',
          'allergen_load_perM','IPI','xenobiotic_degrad','EMHI','USAI',
          'func_redundancy','func_tax_decouple','taxa_core_ubiquity']
LOG_IDX = {'pathogen_relab','amr_total_depth_perM','vf_total_depth_perM',
           'allergen_load_perM','xenobiotic_degrad','allergen_richness'}

def r2_of(y, Xcols_df):
    Xd = sm.add_constant(Xcols_df.values)
    m = sm.OLS(y, Xd).fit()
    return m.rsquared

rows = []
for ix in FUNC16:
    sig = sorted(c for (i2, c) in _confirmed if i2 == ix)  # two-stage-confirmed categories
    y0 = idxc[ix].copy()
    if ix in LOG_IDX: y0 = np.log1p(y0.clip(lower=0))
    blocks = {c: PC1[c] for c in sig}
    blocks['Spatial structure'] = None  # placeholder; MEM block
    # complete cases
    data = pd.DataFrame({'y': y0})
    for c in sig: data[c] = PC1[c]
    data = data.dropna()
    if len(data) < 20:
        rows.append({'index': ix, 'block': 'ALL', 'lmg': np.nan, 'R2_total': np.nan,
                     'n': int(len(data)), 'n_sig_categories': len(sig), 'sig_categories': '|'.join(sig)})
        continue
    yz = (data['y'] - data['y'].mean())/data['y'].std(ddof=0)
    mem = MEM.loc[data.index]
    names = sig + ['Spatial structure']
    def block_df(name):
        return mem if name == 'Spatial structure' else data[[name]]
    # LMG exact over block orderings
    lmg = {b: 0.0 for b in names}
    perms = list(itertools.permutations(names))
    for perm in perms:
        r2_prev = 0.0; cur = []
        for b in perm:
            cur.append(b)
            Xcur = pd.concat([block_df(x) for x in cur], axis=1)
            r2 = r2_of(yz.values, Xcur)
            lmg[b] += (r2 - r2_prev)/len(perms)
            r2_prev = r2
    r2_tot = r2_of(yz.values, pd.concat([block_df(x) for x in names], axis=1))
    for b in names:
        rows.append({'index': ix, 'block': b, 'lmg': lmg[b], 'R2_total': r2_tot,
                     'n': int(len(data)), 'n_sig_categories': len(sig), 'sig_categories': '|'.join(sig)})
att = pd.DataFrame(rows)
wcsv(att, 'attribution_v7.csv', index=False)
print(att[att['block'] != 'Spatial structure'].to_string(index=False), flush=True)
print('\nspatial block LMG per index:', flush=True)
print(att[att['block'] == 'Spatial structure'][['index','lmg','R2_total','n','n_sig_categories']].to_string(index=False), flush=True)

# ---- composites spatial partial Spearman ----
print('\n=== composites ===', flush=True)
def partial_city(xcol, ycol):
    dd = idxc[[xcol, ycol]].dropna()
    M = sm.add_constant(MEM.loc[dd.index].values)
    rx = sst.rankdata(dd[xcol]); ry = sst.rankdata(dd[ycol])
    bx, *_ = np.linalg.lstsq(M, rx, rcond=None); resx = rx - M @ bx
    by, *_ = np.linalg.lstsq(M, ry, rcond=None); resy = ry - M @ by
    r, p = sst.pearsonr(resx, resy)
    r0, p0 = sst.spearmanr(dd[xcol], dd[ycol])
    return len(dd), r0, p0, r, p
crows = []
for x, y in [('USAI','xenobiotic_degrad'), ('IPI','EMHI'), ('EMHI','pathogen_relab')]:
    n, r0, p0, rp, pp = partial_city(x, y)
    crows.append({'x': x, 'y': y, 'n': n, 'rho_raw': r0, 'p_raw': p0, 'rho_spatial': rp, 'p_spatial': pp})
    print(f'  {x} ~ {y}: raw rho={r0:+.3f} (p={p0:.1e}) | spatial rho={rp:+.3f} (p={pp:.1e}), n={n}', flush=True)
wcsv(pd.DataFrame(crows), 'composites_validation_v7.csv', index=False)
print('DONE step7', flush=True)
