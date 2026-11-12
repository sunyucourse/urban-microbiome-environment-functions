#!/usr/bin/env python3
"""v7 step 6: final screens with user-approved specs.
Primary spatial spec: MEM1-10 (city-level factor screen + category screen + sample-level).
Sensitivity: MEM1-5.
Outputs (stats_ms/):
  spatial_full_v7.csv (MEM1-10 primary), spatial_full_v7_mem5.csv (sensitivity)
  category_screen_v7.csv (MEM1-10), category_screen_v7_mem5.csv
  twostage_v7.csv (MEM1-10 nomination + sample-level confirmation), twostage_v7_mem5.csv
"""
import os, warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd, statsmodels.api as sm
from scipy import stats as sst
from statsmodels.stats.multitest import multipletests


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
env = pd.read_csv(f'{BASE}/master_env_city.csv'); env['city'] = env['city'].astype(str).str.strip()
envi = env.set_index('city')
# data fix: zurich station_lat_med erroneously 8.54 (city centroid 47.0 N)
envi.loc['zurich', 'station_lat_med'] = 47.0
co = pd.read_csv(f'{BASE}/city_coords_metasub.csv').set_index('city')
common = [c for c in city.index if c in envi.index and c in co.index]
envc, idxc, coc = envi.loc[common], city.loc[common], co.loc[common]
MEM = pd.DataFrame(dbmem(coc[['lat','lon']].values, 10), index=common,
                   columns=[f'MEM{i+1}' for i in range(10)])

cl = pd.read_csv(f'{MS}/env_clusters_v7.csv')
kept = set(cl.loc[cl['kept'], 'factor']) | set(cl.loc[cl['protected'], 'factor'])
X = envc.select_dtypes(include=[np.number])
X = X.loc[:, X.notna().mean() >= 0.6]; X = X.loc[:, X.std() > 0]
X = X.drop(columns=['n_samples'], errors='ignore')
singletons = [c for c in X.columns if c not in set(cl['factor'])]
allkeep = sorted(kept | set(singletons))
Xz = (X[allkeep] - X[allkeep].mean())/X[allkeep].std(ddof=0)

FUNC16 = ['pathogen_relab','pathway_richness','pathway_evenness','amr_total_depth_perM',
          'amr_class_diversity','vf_total_depth_perM','VPI','allergen_richness',
          'allergen_load_perM','IPI','xenobiotic_degrad','EMHI','USAI',
          'func_redundancy','func_tax_decouple','taxa_core_ubiquity']
TAX = ['species_richness','shannon_tax','invsimpson_tax','pielou_evenness','genus_richness']
LOG_IDX = {'pathogen_relab','amr_total_depth_perM','vf_total_depth_perM',
           'allergen_load_perM','xenobiotic_degrad','allergen_richness'}
# UK pilot batch: mock-community contamination; excluded from pathogen-index analyses only
UK_PILOT = {'lands_end','bury','belfast','eastbourne','kensington','jaywick',
            'newcastle','swansea','edinburgh'}

def mask_uk(y, ix):
    """Mask UK pilot cities to NaN for pathogen-index analyses."""
    if ix == 'pathogen_relab':
        y = y.copy()
        y[y.index.isin(UK_PILOT)] = np.nan
    return y

def city_screen(nmem):
    rows = []
    for ix in FUNC16 + [t for t in TAX if t in idxc.columns]:
        y0 = idxc[ix].copy()
        if ix in LOG_IDX: y0 = np.log1p(y0.clip(lower=0))
        for e in Xz.columns:
            ok = y0.notna() & Xz[e].notna()
            if ok.sum() < 20: continue
            yz = (y0[ok]-y0[ok].mean())/y0[ok].std(ddof=0)
            xe = Xz.loc[ok, e]; mem = MEM.loc[ok].iloc[:, :nmem]
            rho, p = sst.spearmanr(xe, yz)
            Xs = sm.add_constant(np.column_stack([xe.values, mem.values]))
            try:
                m = sm.OLS(yz.values, Xs).fit()
                rows.append({'index': ix, 'env': e, 'n': int(ok.sum()), 'rho': rho, 'p': p,
                             'beta_spatial': m.params[1], 'p_spatial': m.pvalues[1]})
            except Exception: pass
    spa = pd.DataFrame(rows)
    fmask = spa['index'].isin(FUNC16)
    spa['fdr'] = np.nan
    spa.loc[fmask, 'fdr'] = multipletests(spa.loc[fmask, 'p'], method='fdr_bh')[1]
    spa['fdr_spatial'] = np.nan
    spa.loc[fmask, 'fdr_spatial'] = multipletests(spa.loc[fmask, 'p_spatial'], method='fdr_bh')[1]
    return spa

for nmem, tag in [(10, ''), (5, '_mem5')]:
    spa = city_screen(nmem)
    wcsv(spa, f'spatial_full_v7{tag}.csv', index=False)
    fu = spa[spa['index'].isin(FUNC16)]
    print(f'factor screen MEM{nmem}: {len(fu)} functional tests; '
          f'FDR<0.10: {(fu.fdr_spatial<0.10).sum()}; p_spatial<0.05: {(fu.p_spatial<0.05).sum()}', flush=True)

# ---- category screen, both specs ----
PC1 = pd.read_csv(f'{MS}/category_pc1_v7.csv', index_col=0)
def cat_screen(nmem):
    rows = []
    for cat in PC1.columns:
        x = PC1[cat]
        for ix in FUNC16:
            y0 = idxc[ix].copy()
            if ix in LOG_IDX: y0 = np.log1p(y0.clip(lower=0))
            ok = y0.notna() & x.notna()
            if ok.sum() < 20: continue
            yz = (y0[ok]-y0[ok].mean())/y0[ok].std(ddof=0)
            xe = x[ok]; mem = MEM.loc[ok].iloc[:, :nmem]
            rho, p = sst.spearmanr(xe, yz)
            Xs = sm.add_constant(np.column_stack([xe.values, mem.values]))
            try:
                m = sm.OLS(yz.values, Xs).fit()
                rows.append({'category': cat, 'index': ix, 'n': int(ok.sum()), 'rho': rho, 'p': p,
                             'beta_spatial': m.params[1], 'p_spatial': m.pvalues[1]})
            except Exception: pass
    cs = pd.DataFrame(rows)
    cs['fdr_spatial'] = multipletests(cs['p_spatial'].fillna(1), method='fdr_bh')[1]
    cs['fdr_raw'] = multipletests(cs['p'], method='fdr_bh')[1]
    return cs.sort_values('p_spatial')

for nmem, tag in [(10, ''), (5, '_mem5')]:
    cs = cat_screen(nmem)
    wcsv(cs, f'category_screen_v7{tag}.csv', index=False)
    print(f'category screen MEM{nmem}: FDR<0.10: {(cs.fdr_spatial<0.10).sum()}; FDR<0.20: {(cs.fdr_spatial<0.20).sum()}', flush=True)

# ---- two-stage, both specs ----
samp = pd.read_csv(f'{IDX}/indices_sample_v7.csv', index_col=0)
s2 = samp[(samp['source'] == 'metasub') & (~samp['qc_fail'])].copy()
s2 = s2[s2['city'].isin(common)]
for mcol in MEM.columns:
    s2[mcol] = s2['city'].map(MEM[mcol])

def twostage(nmem):
    spa = pd.read_csv(f'{MS}/spatial_full_v7{"" if nmem==10 else "_mem5"}.csv')
    nom = spa[(spa.p_spatial < 0.05) & (spa.env != 'n') &
              spa['index'].isin([f for f in FUNC16 if f in s2.columns])]
    memcols = [f'MEM{i+1}' for i in range(nmem)]
    rows = []
    for _, r in nom.iterrows():
        ix, e = r['index'], r['env']
        s2['xe'] = s2['city'].map(envi[e].astype(float))
        d = s2[[ix, 'xe', 'city'] + memcols].dropna()
        if len(d) < 100: continue
        y = d[ix].copy()
        if ix in LOG_IDX: y = np.log1p(y.clip(lower=0))
        yz = (y - y.mean())/y.std(ddof=0)
        xz = (d['xe'] - d['xe'].mean())/d['xe'].std(ddof=0)
        Xd = sm.add_constant(np.column_stack([xz.values, d[memcols].values]))
        m = sm.OLS(yz.values, Xd).fit(cov_type='cluster', cov_kwds={'groups': d['city'].values})
        rows.append({'index': ix, 'env': e, 'n_samples': len(d), 'n_cities': d['city'].nunique(),
                     'rho_city': r['rho'], 'beta_city': r['beta_spatial'], 'p_city': r['p_spatial'],
                     'fdr_city': r['fdr_spatial'], 'beta_sample': m.params[1], 'p_sample': m.pvalues[1],
                     'direction_match': np.sign(m.params[1]) == np.sign(r['beta_spatial'])})
    sl = pd.DataFrame(rows)
    sl['fdr_sample'] = multipletests(sl['p_sample'], method='fdr_bh')[1]
    return sl.sort_values('p_sample')

for nmem, tag in [(10, ''), (5, '_mem5')]:
    sl = twostage(nmem)
    wcsv(sl, f'twostage_v7{tag}.csv', index=False)
    conf = sl[(sl.fdr_sample < 0.10) & (sl.direction_match)]
    sugg = sl[(sl.fdr_sample >= 0.10) & (sl.fdr_sample < 0.20) & (sl.direction_match)]
    print(f'two-stage MEM{nmem}: nominated {len(sl)}, confirmed {len(conf)}, suggestive {len(sugg)}', flush=True)
print('DONE step6', flush=True)
