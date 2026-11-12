#!/usr/bin/env python
"""Historical statistical follow-up; see README for construct caveats.
v7 step 11: EXTENDED cross-cohort replication (after env extraction extension).
Sewage-side env coverage now: metals(16) + pharma + pm25 + abs_lat + pesticides(3) +
gdp/uhc/pop + GHSL(4) + UHII(2) + WorldClim(3) + CAMS tcno2/pm25 + SoilGrids(2).
Both sides: raw Spearman (rank-based, scale-free). Replication = direction-prespecified
one-sided p < 0.05. Outputs replication_{factor,category,enrichment,concordance}_v8.csv.
"""
import os, warnings; warnings.filterwarnings('ignore')
import pandas as pd, numpy as np
from scipy import stats


from pathlib import Path
import os
PACKAGE = Path(__file__).resolve().parents[1]
BASE = str(PACKAGE / 'data') + '/'
MS = str(Path(os.environ['URBAN_RESULTS_DIR']).resolve()) + '/'
IDX = str(PACKAGE / 'data' / 'indices')
SW = str(PACKAGE / 'data' / 'sewage') + '/'

def wcsv(df, name):
    p = f'{MS}/{name}'
    if os.path.exists(p): os.remove(p)
    df.to_csv(p, index=False)
    print(f'  wrote {name} {df.shape}', flush=True)

# ---------- load discovery side ----------
ts = pd.read_csv(MS + 'twostage_v7.csv')
ts['confirmed'] = (ts.fdr_sample < 0.10) & ts.direction_match
conf = ts[ts.confirmed].copy()
print('confirmed factor pairs:', len(conf))
cat = pd.read_csv(MS + 'category_twostage_v7.csv')
cat_conf = cat[(cat.fdr_sample < 0.10) & (cat.direction_match == True)].copy()
print('confirmed category pairs:', len(cat_conf))

denv = pd.read_csv(BASE + 'master_env_city.csv')
didx = pd.read_csv(BASE + 'indices/indices_city_v7.csv')
print('discovery env cols:', denv.shape, '| indices:', didx.shape)

# ---------- load sewage side ----------
swi = pd.read_csv(SW + 'sewage_city_indices.csv')[
    ['city', 'sw_pathogen_relab', 'sw_amr_depth_pm', 'sw_amr_richness', 'sw_shannon', 'sw_genus_richness']]
met = pd.read_csv(SW + 'gssi_city_metals.csv')
METALS7 = ['As', 'Cd', 'Co', 'Cr', 'Cu', 'Ni', 'Pb']
met_cols = [f'{m}_{s}' for m in METALS7 for s in ('hhet', 'at')] + ['hm_hhet_mean', 'hm_at_mean']
met = met[['city'] + met_cols]
ph = pd.read_csv(SW + 'gssi_city_pharma.csv')
pcol = [c for c in ph.columns if 'log10' in c][0]
ph = ph[['city', pcol]].rename(columns={pcol: 'pharma_log10'})
pm = pd.read_csv(SW + 'gssi_city_pm25.csv')[['city', 'pm25_vd_mean']]
ext = pd.read_csv(SW + 'gssi_city_env_ext.csv')
sw = swi.merge(met, on='city').merge(ph, on='city').merge(pm, on='city').merge(ext, on='city')
sw['abs_lat'] = sw.lat.abs()
sw = sw.replace([np.inf, -np.inf], np.nan)  # log10(0) -> -inf in pharma_cum_log10
print('sewage merged:', sw.shape, flush=True)

IDX_MAP = {'pathogen_relab': 'sw_pathogen_relab',
           'amr_total_depth_perM': 'sw_amr_depth_pm',
           'amr_class_diversity': 'sw_amr_richness'}  # proxy, flagged
ENV_MAP = {'pharma_cum_log10': 'pharma_log10', 'cams_pm25_ugm3_1617': 'pm25_vd_mean',
           'city_lat': 'abs_lat'}
NEW_ENV = ['pest_pesticides', 'pest_fungicides_bactericides', 'pest_rodenticides',
           'gdp_pc_2017', 'uhc_index_2017', 'pop', 'ghs_pop_total', 'ghs_pop_dens',
           'smod_urbancentre_frac', 'smod_urban_frac', 'mean_diurnal_range', 'ann_precip',
           'precip_warmest_q', 'uhii_day', 'uhii_nig', 'cams_tcno2_1617', 'phh2o', 'nitrogen']
for e in NEW_ENV:
    ENV_MAP[e] = e
for m in met_cols:
    ENV_MAP[m] = m
ENV_MAP['pm25_vd_mean'] = 'pm25_vd_mean'

def sw_test(idx, env):
    """Spearman between discovery-index-mapped sewage index and sewage env; one-sided p
    prespecified by discovery direction. Returns dict or None."""
    si, se = IDX_MAP.get(idx), ENV_MAP.get(env)
    if si is None or se is None or se not in sw.columns:
        return None
    d = sw[[si, se]].dropna()
    if len(d) < 20:
        return None
    rho, p2 = stats.spearmanr(d[si], d[se])
    return dict(sw_index=si, sw_env=se, rho_sewage=rho, p_sewage_2sided=p2, n_sewage=len(d))

# ---------- factor-level replication ----------
rows = []
for _, r in conf.iterrows():
    d = dict(index=r['index'], env=r.env, rho_city=r.rho_city, beta_sample=r.beta_sample,
             fdr_sample=r.fdr_sample)
    t = sw_test(r['index'], r.env)
    if t is None:
        d.update(testable=False, replicated=False)
    else:
        disc_sign = np.sign(r.beta_sample)
        match = bool(np.sign(t['rho_sewage']) == disc_sign)
        p1 = t['p_sewage_2sided'] / 2 if match else 1 - t['p_sewage_2sided'] / 2
        d.update(testable=True, **t, direction_match_sewage=match,
                 p_sewage_1sided=p1, replicated=bool(p1 < 0.05))
    rows.append(d)
rep = pd.DataFrame(rows)
wcsv(rep, 'replication_factor_v8.csv')
tb = rep[rep.testable]
print(f'\n=== FACTOR: {len(tb)}/{len(rep)} testable; {tb.replicated.sum()} replicated (1-sided p<0.05) ===')
print(tb[['index', 'env', 'beta_sample', 'fdr_sample', 'rho_sewage', 'p_sewage_1sided', 'n_sewage', 'replicated']].to_string())

# ---------- category-level replication with sewage-side PC1s ----------
def pc1(dfm, members):
    Z = dfm[members].apply(pd.to_numeric, errors='coerce')
    Z = (Z - Z.mean()) / Z.std(ddof=0)
    cc = Z.dropna()
    if len(cc) < 20:
        return None
    U, S, Vt = np.linalg.svd(cc.values - cc.values.mean(0), full_matrices=False)
    sc = pd.Series(U[:, 0] * S[0], index=cc.index)
    if np.corrcoef(sc, cc.mean(axis=1))[0, 1] < 0:
        sc = -sc
    return sc

CAT_MEM_SW = {
    'Heavy metals': met_cols,
    'Pharmaceuticals': ['pharma_log10'],
    'PM2.5 (satellite)': ['pm25_vd_mean'],
    'Pesticides': ['pest_pesticides', 'pest_fungicides_bactericides', 'pest_rodenticides'],
    'Socioeconomic': ['gdp_pc_2017', 'pop', 'uhc_index_2017'],  # popdens unavailable; validated below
    'Built environment': ['ghs_pop_dens', 'smod_urban_frac'],
    'Urban heat island': ['uhii_day', 'uhii_nig'],
    'Geographic': ['abs_lat'],
}
sw_scores = {}
for cname, mem in CAT_MEM_SW.items():
    if len(mem) == 1:
        s = sw.set_index('city')[mem[0]]
        sw_scores[cname] = (s - s.mean()) / s.std(ddof=0)
    else:
        sw_scores[cname] = pc1(sw.set_index('city'), mem)
    n = sw_scores[cname].notna().sum() if sw_scores[cname] is not None else 0
    print(f'sewage PC1 {cname}: n={n}', flush=True)

# discovery-side validation: Socioeconomic 3-member vs full 4-member PC1
full = pc1(denv.set_index('city'), ['gdp_pc_2017', 'pop', 'popdens', 'uhc_index_2017'])
red = pc1(denv.set_index('city'), ['gdp_pc_2017', 'pop', 'uhc_index_2017'])
common = full.index.intersection(red.index)
r_v, p_v = stats.spearmanr(full[common], red[common])
print(f'Socioeconomic PC1: full(4) vs reduced(3, no popdens): rho={r_v:.3f} p={p_v:.1e} n={len(common)}', flush=True)

# PC1 loading diagnostics (both cohorts) for Socioeconomic / BuiltEnv / UHI / Pesticides
def pc1_loadings(dfm, members):
    Z = dfm[members].apply(pd.to_numeric, errors='coerce')
    Z = (Z - Z.mean()) / Z.std(ddof=0)
    cc = Z.dropna()
    U, S, Vt = np.linalg.svd(cc.values - cc.values.mean(0), full_matrices=False)
    sc = pd.Series(U[:, 0] * S[0], index=cc.index)
    ld = pd.Series(Vt[0], index=members)
    if np.corrcoef(sc, cc.mean(axis=1))[0, 1] < 0:
        ld = -ld
    return ld
for cname, mem_d, mem_s in [
    ('Socioeconomic', ['gdp_pc_2017', 'pop', 'popdens', 'uhc_index_2017'], ['gdp_pc_2017', 'pop', 'uhc_index_2017']),
    ('Built environment', ['ghs_pop_dens', 'smod_urban_frac'], ['ghs_pop_dens', 'smod_urban_frac']),
    ('Urban heat island', ['uhii_day', 'uhii_nig'], ['uhii_day', 'uhii_nig']),
    ('Pesticides', ['pest_pesticides', 'pest_fungicides_bactericides', 'pest_rodenticides'],
     ['pest_pesticides', 'pest_fungicides_bactericides', 'pest_rodenticides'])]:
    ld = pc1_loadings(denv.set_index('city'), mem_d)
    ls = pc1_loadings(sw.set_index('city'), mem_s)
    print(f'  loadings {cname}: discovery={ld.round(2).to_dict()} | sewage={ls.round(2).to_dict()}', flush=True)

crows = []
for _, r in cat_conf.iterrows():
    cname = r['category']
    d = dict(category=cname, index=r['index'], beta_sample=r.beta_sample, fdr_sample=r.fdr_sample)
    si = IDX_MAP.get(r['index'])
    sc = sw_scores.get(cname)
    if si is None or sc is None:
        d.update(testable=False, replicated=False)
    else:
        tmp = pd.DataFrame({'score': sc}).join(sw.set_index('city')[[si]]).dropna()
        if len(tmp) < 20:
            d.update(testable=False, replicated=False)
        else:
            rho, p2 = stats.spearmanr(tmp.score, tmp[si])
            match = bool(np.sign(rho) == np.sign(r.beta_sample))
            p1 = p2 / 2 if match else 1 - p2 / 2
            d.update(testable=True, sw_index=si, rho_sewage=rho, p_sewage_2sided=p2,
                     n_sewage=len(tmp), direction_match_sewage=match,
                     p_sewage_1sided=p1, replicated=bool(p1 < 0.05))
    crows.append(d)
crep = pd.DataFrame(crows)
wcsv(crep, 'replication_category_v8.csv')
ct = crep[crep.testable]
print(f'\n=== CATEGORY: {len(ct)}/{len(crep)} testable; {ct.replicated.sum()} replicated ===')
print(ct[['category', 'index', 'beta_sample', 'fdr_sample', 'rho_sewage', 'p_sewage_1sided', 'n_sewage', 'replicated']].to_string())

# ---------- expanded concordance: all sewage-available factors x 5 shared indices ----------
CONC_IDX = {'pathogen_relab': 'sw_pathogen_relab', 'amr_total_depth_perM': 'sw_amr_depth_pm',
            'amr_gene_richness': 'sw_amr_richness', 'shannon_tax': 'sw_shannon',
            'genus_richness': 'sw_genus_richness'}
# discovery env factors available in sewage (same or mapped names)
D2S = {'pharma_cum_log10': 'pharma_log10', 'cams_pm25_ugm3_1617': 'pm25_vd_mean',
       'city_lat': 'abs_lat', 'pm25_vd_mean': 'pm25_vd_mean'}
for e in NEW_ENV + met_cols:
    D2S[e] = e
didx_i = didx.set_index('city')
denv_i = denv.set_index('city')
denv_i['abs_lat'] = denv_i.city_lat.abs()

crows2 = []
for di, si in CONC_IDX.items():
    if di not in didx_i.columns:
        print('missing discovery index:', di); continue
    for denv_name, sw_env in D2S.items():
        if denv_name not in denv_i.columns or sw_env not in sw.columns:
            continue
        dd = pd.DataFrame({'x': denv_i[denv_name], 'y': didx_i[di]}).dropna()
        if len(dd) < 20: continue
        r_d, p_d = stats.spearmanr(dd.x, dd.y)
        sd = sw[[si, sw_env]].dropna()
        if len(sd) < 20: continue
        r_s, p_s = stats.spearmanr(sd[si], sd[sw_env])
        crows2.append(dict(index_discovery=di, index=si, env_discovery=denv_name, env=sw_env,
                           n_discovery=len(dd), rho_discovery=r_d, p_discovery=p_d,
                           n_sewage=len(sd), rho_sewage=r_s, p_sewage=p_s))
conc = pd.DataFrame(crows2)
# drop duplicate sewage-variable mapping: cams_pm25_ugm3_1617 and pm25_vd_mean both map to
# sewage pm25_vd_mean; keep the direct satellite-PM2.5 construct (as in the v7 95-pair set)
conc = conc[~((conc.env_discovery == 'cams_pm25_ugm3_1617') & (conc.env == 'pm25_vd_mean'))]
wcsv(conc, 'replication_concordance_v8.csv')
r_all, p_all = stats.spearmanr(conc.rho_discovery, conc.rho_sewage)
print(f'\n=== CONCORDANCE: {len(conc)} pairs ({conc.env.nunique()} factors x {conc["index"].nunique()} indices); rho={r_all:.3f}, p={p_all:.1e} ===')
# old 19-factor subset for continuity with v7
OLD19 = ['pm25_vd_mean', 'pharma_log10', 'abs_lat'] + met_cols
sub = conc[conc.env.isin(OLD19)]
r_s, p_s2 = stats.spearmanr(sub.rho_discovery, sub.rho_sewage)
print(f'  old-19-factor subset ({len(sub)} pairs): rho={r_s:.3f}, p={p_s2:.1e}')

# AMR-depth axis direction consistency
amr = conc[conc['index'] == 'sw_amr_depth_pm']
same = (np.sign(amr.rho_sewage) == np.sign(amr.rho_discovery)).sum()
pb = stats.binomtest(int(same), len(amr), 0.5).pvalue
print(f'AMR-depth axis: {same}/{len(amr)} direction-consistent, binomial p={pb:.1e}')
# metals subset
metmask = amr.env.isin(met_cols)
msame = (np.sign(amr[metmask].rho_sewage) == np.sign(amr[metmask].rho_discovery)).sum()
msig = ((amr[metmask].p_sewage < 0.001) & (np.sign(amr[metmask].rho_sewage) == np.sign(amr[metmask].rho_discovery))).sum()
print(f'  metals->AMR depth: {msame}/{metmask.sum()} direction-consistent; {msig} with sewage p<0.001 same direction')
# discovery-sig pairs replication (uniform raw p<0.05 both sides)
conc['disc_sig'] = conc.p_discovery < 0.05
conc['sw_rep'] = (conc.p_sewage < 0.05) & (np.sign(conc.rho_sewage) == np.sign(conc.rho_discovery))
tab = pd.crosstab(conc.disc_sig, conc.sw_rep)
print('\ndiscovery-sig x sewage-replicated:'); print(tab)
if tab.shape == (2, 2):
    OR, pf = stats.fisher_exact(tab)
    print(f'Fisher OR={OR:.2f}, p={pf:.4f}')
ds = conc[conc.disc_sig]
sd_same = (np.sign(ds.rho_sewage) == np.sign(ds.rho_discovery)).sum()
print(f'discovery-sig pairs: {len(ds)}; same direction in sewage: {sd_same} '
      f'(binom p={stats.binomtest(int(sd_same), len(ds), 0.5).pvalue:.4f})')
wcsv(conc, 'replication_enrichment_v8.csv')
