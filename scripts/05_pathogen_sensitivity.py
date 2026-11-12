#!/usr/bin/env python3
"""Pathogen dominance sensitivity: new supp table S14 + renumber S14-S18 -> S15-S19 + manuscript text."""
import pandas as pd, numpy as np
from scipy import stats

from pathlib import Path
import os
PACKAGE = Path(__file__).resolve().parents[1]
BASE = str(PACKAGE / 'data')
MS = str(Path(os.environ['URBAN_RESULTS_DIR']).resolve())
IDX = str(PACKAGE / 'data' / 'indices')
SW = str(PACKAGE / 'data' / 'sewage') + '/'

B = BASE
SP22 = ['Stenotrophomonas_maltophilia','Klebsiella_pneumoniae','Acinetobacter_baumannii',
 'Enterobacter_cloacae','Klebsiella_oxytoca','Proteus_mirabilis','Citrobacter_freundii',
 'Citrobacter_koseri','Moraxella_catarrhalis','Enterococcus_faecium','Haemophilus_influenzae',
 'Morganella_morganii','Streptococcus_pneumoniae','Serratia_marcescens','Burkholderia_multivorans',
 'Vibrio_cholerae','Klebsiella_aerogenes','Burkholderia_cenocepacia','Aeromonas_hydrophila',
 'Vibrio_natriegens','Campylobacter_jejuni','Clostridium_perfringens']

# ---------- 1. build S14 table ----------
v7s = pd.read_csv(f'{B}/indices/indices_sample_v7.csv')
v7c = pd.read_csv(f'{B}/indices/indices_city_v7.csv')
d = v7s[v7s['city'].isin(set(v7c.city)) & (v7s.qc_fail == 0) & (v7s.source == 'metasub')].copy()
sp = pd.read_csv(f'{B}/pathogen_species22.csv', index_col=0)
X = sp.loc[d['Unnamed: 0'], SP22]; X.index = d.index
d['path_rich'] = (X.values > 0).sum(axis=1)
d['path_logsum'] = np.log1p(X.values).sum(axis=1)
d['path_sqrtsum'] = np.sqrt(X.values).sum(axis=1)
d['smalto'] = X['Stenotrophomonas_maltophilia'].values
city = d.groupby('city')[['pathogen_relab','path_rich','path_logsum','path_sqrtsum','smalto']].mean()
env = pd.read_csv(f'{B}/master_env_city.csv')
m = city.merge(env.set_index('city'), left_index=True, right_index=True)

ts = pd.read_csv(f'{MS}/twostage_v7.csv')
tsp = ts[ts['index'] == 'pathogen_relab'].set_index('env')
factors = ['city_lat','phh2o','pest_fungicides_bactericides','pop','ghs_pop_dens',
           'mean_diurnal_range','cams_tcno2_1617','precip_warmest_q','smod_urban_frac',
           'ann_precip','nitrogen']
variants = [('abundance_sum','pathogen_relab'), ('richness_presence_absence','path_rich'),
            ('log_sum','path_logsum'), ('sqrt_sum','path_sqrtsum'), ('s_maltophilia_only','smalto')]
rows = []
for f in factors:
    row = {'factor': f,
           'screen_beta_city': round(tsp.loc[f, 'beta_city'], 3),
           'screen_p_city': f"{tsp.loc[f, 'p_city']:.2e}",
           'screen_FDR_sample': f"{tsp.loc[f, 'fdr_sample']:.2e}"}
    for vname, col in variants:
        dd = m.dropna(subset=[col, f])
        r = stats.spearmanr(dd[f], dd[col])
        row[f'rho_{vname}'] = round(r.statistic, 3)
        row[f'p_{vname}'] = f'{r.pvalue:.2e}'
        row['n_cities'] = len(dd)
    rows.append(row)
s14 = pd.DataFrame(rows)
s14.to_csv(f'{MS}/S14_pathogen_dominance.csv', index=False)
print('S14 table:'); print(s14[['factor','rho_abundance_sum','rho_richness_presence_absence','rho_s_maltophilia_only']].to_string(index=False))

