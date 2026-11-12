#!/usr/bin/env python3
"""v7 step6b: category-level two-stage (nominate from MEM10 category screen p_spatial<0.05;
confirm at sample level with city-clustered OLS + MEM1-10, FDR<0.10 + direction).
Writes category_twostage_v7.csv. Run AFTER step4 (category PC1s) and step6 (category screen).
"""
import pandas as pd, numpy as np
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from scipy import stats as sst


from pathlib import Path
import os
PACKAGE = Path(__file__).resolve().parents[1]
BASE = str(PACKAGE / 'data')
MS = str(Path(os.environ['URBAN_RESULTS_DIR']).resolve())
IDX = str(PACKAGE / 'data' / 'indices')
SW = str(PACKAGE / 'data' / 'sewage') + '/'

def wcsv(df, name, **kw):
    import os
    p = f'{MS}/{name}'
    if os.path.exists(p):
        os.remove(p)
    df.to_csv(p, **kw)

FUNC16 = ['pathogen_relab','pathway_richness','pathway_evenness','amr_total_depth_perM',
          'amr_class_diversity','vf_total_depth_perM','VPI','allergen_richness',
          'allergen_load_perM','IPI','xenobiotic_degrad','EMHI','USAI',
          'func_redundancy','func_tax_decouple','taxa_core_ubiquity']
LOG_IDX = {'pathogen_relab','amr_total_depth_perM','vf_total_depth_perM',
           'allergen_load_perM','xenobiotic_degrad','allergen_richness'}

# city-level MEMs (same construction as step6)
city = pd.read_csv(f'{IDX}/indices_city_v7.csv').set_index('city')
co = pd.read_csv(f'{BASE}/city_coords_metasub.csv').set_index('city')
common = [c for c in city.index if c in co.index]

def haversine(lat1, lon1, lat2, lon2):
    R = 6371.0; p1, p2 = np.radians(lat1), np.radians(lat2)
    dp, dl = np.radians(lat2-lat1), np.radians(lon2-lon1)
    a = np.sin(dp/2)**2 + np.cos(p1)*np.cos(p2)*np.sin(dl/2)**2
    return 2*R*np.arcsin(np.sqrt(a))

def dbmem(coords, n_axes=10):
    # identical to v7_step6_finalscreens.py: haversine PCNM, no truncation
    n = len(coords); D = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            D[i, j] = D[j, i] = haversine(coords[i,0], coords[i,1], coords[j,0], coords[j,1])
    A = -0.5*D**2; Bm = A - A.mean(0) - A.mean(1)[:, None] + A.mean()
    ev, ec = np.linalg.eigh(Bm); o = np.argsort(ev)[::-1]
    ev, ec = ev[o], ec[:, o]
    return ec[:, ev > 0][:, :n_axes]

MEM = pd.DataFrame(dbmem(co.loc[common, ['lat','lon']].values, 10), index=common,
                   columns=[f'MEM{i+1}' for i in range(10)])

PC1 = pd.read_csv(f'{MS}/category_pc1_v7.csv', index_col=0)
cs = pd.read_csv(f'{MS}/category_screen_v7.csv')  # MEM10 primary

samp = pd.read_csv(f'{IDX}/indices_sample_v7.csv', index_col=0)
s2 = samp[(samp['source'] == 'metasub') & (~samp['qc_fail'])].copy()
s2 = s2[s2['city'].isin(common)]
for mcol in MEM.columns:
    s2[mcol] = s2['city'].map(MEM[mcol])

nom = cs[(cs.p_spatial < 0.05) & (cs['index'].isin([f for f in FUNC16 if f in s2.columns]))]
memcols = [f'MEM{i+1}' for i in range(10)]
rows = []
for _, r in nom.iterrows():
    ix, cat = r['index'], r['category']
    s2['xe'] = s2['city'].map(PC1[cat].astype(float))
    d = s2[[ix, 'xe', 'city'] + memcols].dropna()
    if len(d) < 100:
        continue
    y = d[ix].copy()
    if ix in LOG_IDX:
        y = np.log1p(y.clip(lower=0))
    yz = (y - y.mean()) / y.std(ddof=0)
    xz = (d['xe'] - d['xe'].mean()) / d['xe'].std(ddof=0)
    Xd = sm.add_constant(np.column_stack([xz.values, d[memcols].values]))
    m = sm.OLS(yz.values, Xd).fit(cov_type='cluster', cov_kwds={'groups': d['city'].values})
    rows.append({'category': cat, 'index': ix, 'n_samples': len(d),
                 'beta_city': r['beta_spatial'], 'p_city': r['p_spatial'],
                 'beta_sample': m.params[1], 'p_sample': m.pvalues[1],
                 'direction_match': np.sign(m.params[1]) == np.sign(r['beta_spatial'])})
sl = pd.DataFrame(rows)
sl['fdr_sample'] = multipletests(sl['p_sample'], method='fdr_bh')[1]
sl = sl.sort_values('p_sample')
wcsv(sl, 'category_twostage_v7.csv', index=False)
conf = sl[(sl.fdr_sample < 0.10) & sl.direction_match]
sugg = sl[(sl.fdr_sample >= 0.10) & (sl.fdr_sample < 0.20) & sl.direction_match]
print(f'category two-stage MEM10: nominated {len(sl)}, confirmed {len(conf)}, suggestive {len(sugg)}')
print(sl.to_string())
print('DONE step6b')
