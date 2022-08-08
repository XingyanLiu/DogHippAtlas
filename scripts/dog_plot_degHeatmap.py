# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 22:03:53 2022

@author: xyliu
"""

import os
from pathlib import Path
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

WORKDIR = Path(r'D:\Users\xyliu\003')
os.chdir(WORKDIR)

from src import funx as fx, funx_dog as fxd
from src import build as B

'''     Settings
==========================
'''
# sub = 'ppca'
# resdir = fxd.HippDir / f'Analysis_{sub}'
resdir = fxd.HippDirFormal
figdir = resdir / 'figs'  # Path('figures/dog')
fx.check_dirs(figdir)

sc.settings.figdir = figdir
sc.set_figure_params(fontsize=14, )

''' Loading dataset
'''
adata = fxd.DataHippFormal()
adata

# In[]

_gmap = {'ENSCAFG00000030472': 'NRXN3'}
adata.var_names = [_gmap.get(x, x) for x in adata.var_names]
# adata.raw.var_names = [_gmap.get(x, x) for x in adata.raw.var_names]

# In[] 

cate_plt = [
    'Glutamatergic neurons', 'GABAergic neurons', 'Oligodendrocytes',
    #       'Unknown',
    'Myelinating oligodendrocytes', 'Ependymal cells',
    'Endothelial cells', 'Microglials', 'Astrocytes', 'Cajal-Retzius']

key_major = 'major_class'

avgs_major = B.GroupMean(adata, key_major, use_raw=True)

avgs_major.index = [_gmap.get(x, x) for x in avgs_major.index]
# avgs_major_maxnorm = avgs_major.apply(lambda x: x / x.max(), axis=1)
avgs_major.to_csv(resdir / 'mean_expressions-major_class.csv')

# In[]
"""all degs"""
dir_de = Path(r'E:\Users\xyliu\data003\dog\DE\20210902')

df = pd.read_csv(dir_de / 'deg_intersects-p0.001.tsv', sep='\t',
                 header=None, index_col=0)
dedict = df.iloc[:, 0].apply(lambda x: x.split(',')).to_dict()
unique_degs = set()
for k, lst in dedict.items():
    unique_degs.update(lst)

print(len(unique_degs))
tag = 'withUNK'

# In[]
''' heatmap
'''
# genes_hmap = fx.load_namelist(fxd.GENEDIR / 'genes_heatmap.txt')
genes_hmap = unique_degs

df_hmap = avgs_major.loc[genes_hmap, cate_plt]
df_hmap = df_hmap.apply(lambda x: x / x.max(), axis=1)

df_hmap1 = fx.order_contingency_df(df_hmap, axis=1)

df_hmap1[df_hmap1.isna().any(1)]

df_hmap1.to_csv(resdir / 'heatmap_values-mean-major_class-maxNorm.csv')

# In[]
sc.set_figure_params(dpi_save=180, fontsize=11, )
cmap_heat = ['RdBu_r', 'magma_r'][0]

fig, ax = plt.subplots(figsize=(4, 12))
sns.heatmap(df_hmap1,
            yticklabels=False,
            #            cbar=False,
            cbar_kws={'shrink': 0.3},
            ax=ax,
            cmap=cmap_heat,
            #            col_colors = cl_color_match,
            )
fig.savefig(figdir / f'heatmap_plain-{cmap_heat}-220314.pdf',
            bbox_inches='tight')

# In[]
''' heatmap with UNK (from scratch)
'''
dir_de = Path(r'E:\Users\xyliu\data003\dog\DE\20210902')

df = pd.read_csv(dir_de / 'deg_intersects-p0.001.tsv', sep='\t',
                 header=None, index_col=0)
dedict = df.iloc[:, 0].apply(lambda x: x.split(',')).to_dict()
unique_degs = set()
for k, lst in dedict.items():
    unique_degs.update(lst)
print(len(unique_degs))
tag = 'withUNK'

# genes_hmap = fx.load_namelist(fxd.GENEDIR / 'genes_heatmap.txt')
genes_hmap = unique_degs
cate_all = adata.obs[key_major].cat.categories
df_hmap = avgs_major.loc[genes_hmap, cate_all]
df_hmap = df_hmap.apply(lambda x: x / x.max(), axis=1)

df_hmap1 = fx.order_contingency_df(df_hmap, axis=1)

df_hmap1[df_hmap1.isna().any(1)]

df_hmap1.to_csv(resdir / 'heatmap_values-mean-major_class-maxNorm-{tag}.csv')

# 
sc.set_figure_params(dpi_save=180, fontsize=11, )
cmap_heat = ['RdBu_r', 'magma_r'][0]

fig, ax = plt.subplots(figsize=(4, 12))
sns.heatmap(df_hmap1,
            yticklabels=False,
            #            cbar=False,
            cbar_kws={'shrink': 0.3},
            ax=ax,
            cmap=cmap_heat,
            #            col_colors = cl_color_match,
            )
fig.savefig(figdir / f'heatmap_plain-{cmap_heat}-{tag}-220314.pdf',
            bbox_inches='tight')
