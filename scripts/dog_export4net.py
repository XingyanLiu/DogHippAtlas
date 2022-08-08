# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 22:11:24 2020

@author: xyliu
"""

import sys
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

sys.path.append('../')

from src import funx as fx, funx_dog as fxd
from src import build as B
from src import preproc as pp

# settings load annotated data (normalized)

resdir = fxd.DATADIR_NEW / 'hipp_subset'
figdir = resdir / 'figs'
pp.check_dirs(figdir)

adata = fxd.DataHippFormal()
hvgs = adata.var_names

# 0.1 select genes (DEGs)
n_top = 50
marker_unique = pp.top_markers_from_adata(adata, n_top)

# expression averages and proportions
expr_avgs = B.GroupMean(adata, 'leiden', binary=False, use_raw=True)
expr_props = B.GroupMean(adata, 'leiden', binary=True, use_raw=True)

expr_avgs.to_csv(fxd.HippDirFormal / 'expr_avgs.csv')
expr_props.to_csv(fxd.HippDirFormal / 'expr_props.csv')


# 0.2 select genes (screening by max_proportion)
# max proportion for each gene
max_prop = expr_props.max(axis=1)
plt.hist(max_prop, bins=20)

# filter genes
cut_prop = 0.25
passed = max_prop >= cut_prop
(passed).sum()
genes_pass = max_prop[passed].index.tolist()


# inspect intersection
fx.VennPlot([marker_unique, hvgs, genes_pass])

# save genes used for WGCAN
fx.save_namelist(genes_pass, 
                 fxd.GENEDIR / f'screened_genes_cut{cut_prop}.csv')


# take out subset gene-expression matrix
adata_reset = adata.raw.to_adata()
adata_sub = adata_reset[:, genes_pass]


## save
#fx.saveMtx2df(adata_sub, resdir / f'sc_expr_cutprop{cut_prop}.csv')

# Downsampling of single cells

adata.obs['leiden_anno'].value_counts()
sample_frac = 0.2
adata_sub_sub = sc.pp.subsample(adata_sub, fraction=sample_frac, copy=True)
adata_sub_sub.obs['leiden_anno'].value_counts()

fx.saveMtx2df(adata_sub_sub, resdir / f'sc_expr_sub_cutprop{cut_prop}.csv')


# recurrent clustering - to break the whole dataset into many small clusters
key_orig = 'leiden'
key_new = f'{key_orig}_lv2'
max_cnt = 100
reso = 1
adata.obs[key_new] = adata.obs[key_orig]
grp_cnts = adata.obs[key_new].value_counts()
while True:
    candi_grp = grp_cnts.idxmax()
    _n_cells = grp_cnts[candi_grp]
    print(_n_cells)
    if _n_cells > max_cnt:
        print(f'breaking cluster `{candi_grp}`...')
        sc.tl.leiden(adata, resolution=reso,
                     restrict_to=(key_new, [candi_grp]),
                     key_added=key_new)
        grp_cnts = adata.obs[key_new].value_counts()
    else:
        print(f'no cluster contain more than {max_cnt} cells, stop.')
        break
    
len(grp_cnts)


# In[] ###################
tt = 'number of cells in each aggregation'
grp_cnts.hist()
plt.title(tt)
plt.xlabel('number of cells')
plt.ylabel('frequency')
plt.savefig(figdir/ f'{tt}.pdf', bbox_inches='tight')


# change the names of the categories
new_cats = list(map(lambda x: x.replace(',', '_'), 
                    adata.obs[key_new].cat.categories))
adata.obs[key_new].cat.categories = new_cats

# aggragate cells in small clusters - remove batch effect by the way :)
expr_agg = B.GroupMean(adata, groupby=key_new, use_raw=True)


expr_agg.to_csv(resdir / f'agg_expr_all.csv', header=True, index=True)

expr_agg.loc[genes_pass, :].to_csv(
        resdir / f'agg_expr_cutprop{cut_prop}.csv',
        header=True, index=True)

# make metadata for aggregated expression
# cell types and colors
df_color_match = pd.read_csv(fxd.HippDirFormal / 'type_color_match.csv',
                             header=None, )
df_color_match.columns = ['leiden_anno', 'color']
df_color_match['leiden'] = df_color_match['leiden_anno'].apply(
              lambda x:x.split(' ')[0])

### leiden cluster number
lbs_leiden = list(map(lambda x: x.split('_')[0], expr_agg.columns))
pd.value_counts(lbs_leiden)

### merge those 2 annotations
meta_agg = pd.DataFrame({'agg_name': expr_agg.columns, 'leiden': lbs_leiden}, )
meta_agg = pd.merge(meta_agg, df_color_match, 
                    how='left', on='leiden')
meta_agg.head()
meta_agg['leiden_anno'].value_counts()

meta_agg.to_csv(resdir / f'agg_expr_cutprop{cut_prop}_metadata.csv',
        header=True, index=False)

