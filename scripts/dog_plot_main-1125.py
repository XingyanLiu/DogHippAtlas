# -*- coding: utf-8 -*-
'''
Created on Tue May  5 18:02:06 2020

@author: xyliu
'''

import os
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

WORKDIR = Path(r'D:\Users\xyliu\003')
os.chdir(WORKDIR)

from src import funx as fx, funx_dog as fxd
from src import build as B

# In[]
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
#### or load formal data (annotated)
# adata = sc.read_h5ad(fxd.HippDir / 'Merged_formal.h5ad')
###########[ export metadata ]###############
# metadata = adata.obs[['batch', 'n_genes', 'n_counts', 'hvg_counts',
#                      'primer', 'leiden', 'leiden_anno', 'major_class']].copy()
# X_umap = adata.obsm['X_umap']
# metadata['UMAP1'] = X_umap[:, 0]
# metadata['UMAP2'] = X_umap[:, 1]
# metadata.to_csv(resdir / f'metadata.csv', index_label='barcode')

# In[] 
''' Data overview '''

adata.obs['primer'] = adata.obs['RNAcap']
adata.obs['primer'].cat.categories = ['polyT', 'random']

summ_g = adata.obs.groupby('primer')['n_genes'].describe()
summ_c = adata.obs.groupby('primer')['n_counts'].describe()
# summ_g.to_csv(fxd.DATADIR_NEW / 'summ_gene.csv', index=True, header=True)
# summ_c.to_csv(fxd.DATADIR_NEW / 'summ_count.csv', index=True, header=True)
# summ_g.plot.hist()
# summ = pd.concat


sc.pl.violin(adata, ['n_genes', 'n_counts'], groupby='primer', log=True)

# add annotations

AnnoDF = pd.read_csv(fxd.HippDirFormal / 'cluster_annotations.csv')
AnnoDF.columns
AnnoDF['leiden'] = AnnoDF['leiden'].astype(str)
leiden_anno = AnnoDF[['leiden', 'cell_type']].apply(
    lambda x: ' '.join(x),
    axis=1)

# clusterNames = adata.obs['leiden'].astype(int).map(leiden_anno.to_dict())
adata.obs['leiden'].cat.categories
adata.obs['leiden_anno'] = adata.obs['leiden'].copy()
adata.obs['leiden_anno'].cat.categories = leiden_anno.values

# take out cluster colors
cl_colors = adata.uns['leiden_colors']
cl_color_match = pd.Series(cl_colors, index=leiden_anno)
print(cl_color_match)

''' save annotated data
'''
formal_name = 'analyzed'
# adata.obs.to_csv(fxd.HippDirFormal / f'{formal_name}_metadata.csv')
# cl_color_match.to_csv(fxd.HippDirFormal / 'type_color_match.csv', index=True, header=False)
##adata.write(fxd.HippDirFormal / f'{formal_name}.h5ad')

# In[]
''' annotation of the major types
'''
dct_major = {
    'Glutamatergic neurons': [0, 17, 11, 14, 3, 5, 6, 7, 10, 13, 15, 20],
    'GABAergic neurons': [4, 8, 12, 21],
    'Oligodendrocytes': [9],
    'Unknown': [16],  # 16 is removed
    'Myelinating oligodendrocytes': [2, 18],
    'Ependymal cells': [19, 22],
    'Endothelial cells': [23],
    'Microglials': [24],
    'Astrocytes': [1],
    'Cajal-Retzius': [25]}

ordered_major_types = dct_major.keys()
dct_major_rev = fx.reverse_dict_of_list(dct_major)

major_colors = fx.get_colors('tab20b', n=len(dct_major))
fx.view_color_map(major_colors)

key_major = 'major_class'

major_colors_srs = pd.Series(major_colors, index=ordered_major_types,
                             name='color')
major_colors_srs.to_csv(
    fxd.HippDirFormal / 'major_class_colors.csv', index_label=key_major,
    header=True)

adata.obs[key_major] = adata.obs['leiden'].astype(str).apply(
    lambda x: dct_major_rev[int(x)])
adata.obs[key_major] = pd.Categorical(
    adata.obs[key_major], categories=ordered_major_types)
adata.uns[f'{key_major}_colors'] = major_colors

sc.pl.umap(adata, color=key_major)  # test for label and colors

# In[]
''' UMAP vis cell types '''

n_cells = adata.shape[0]
tt = f'Hippocampus ({n_cells} cells)'
sc.set_figure_params(dpi_save=150, fontsize=14, )

sc.pl.umap(adata, color='leiden_anno', title=tt, save=f'_.pdf', )
sc.pl.umap(adata, color='leiden', legend_loc='on data',
           legend_fontsize=12,
           title=tt,
           save=f'_ondata.pdf',
           )

# In[]

''' inspect some marker genes on UMAP
'''
figdir = resdir / 'figs-gene_umap_all'  # Path('figures/dog')
fx.check_dirs(figdir)
sc.settings.figdir = figdir
sc.set_figure_params(dpi_save=150, fontsize=14, )
# plt.hist(adata.raw.X.data, bins=50)
# genes_all = adata.raw.var_names.tolist()
# 'CUX2' in genes_all

# _gene = ['MEGF11', 'ARHGAP24', 'PTPRZ1'][-1]
# genes4umap = pd.read_csv(fxd.GENEDIR / 'genes4umap.csv', header=None)
genes4umap = pd.read_csv(fxd.GENEDIR / '1105-gene_list.txt', header=None)
g4umap = genes4umap.iloc[:, 0]

# _genes = ['CUX2', 'CCBE1', 'NDNF', 'MBP']
# _genes = fxd.GenesDomastic + ['NDNF']
# genes_plot = [g for g in g4umap if g not in _genes]
genes_plot = adata.raw.var_names.tolist()
######################################################
cmap_n = ['PuRd', 'Reds', 'YlGnBu', 'GnBu'][-1]
cmap_n = 'diy'
if cmap_n != 'diy':
    cmap_g = cmap_n
else:
    candi_colors = ['#d9d9d9'] * 15 + fx.get_colors('RdPu', 100)[15:]
    #    candi_colors = ['#d9d9d9']*15 + fx.get_colors('PuRd', 100)[20: -15]
    cmap_g = mpl.colors.ListedColormap(candi_colors)
#    nodes = [0., 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
#    list(zip(nodes, candi_colors))
#    camp_g = mpl.colors.LinearSegmentedColormap.from_list('diy', list(zip(nodes, candi_colors)))
######################################################
name_dict = {'ENSCAFG00000030472': 'NRXN3'}
cut = False
for _g in genes_plot:
    gexpr = adata.raw[:, _g].X.flatten()
    _vmax = np.quantile(gexpr[gexpr > 0], 0.99) if cut else None
    gn = name_dict.get(_g, _g)

    sc.pl.umap(adata, color=_g,
               #               vmax=_vmax,
               size=5,
               title=gn,
               cmap=cmap_g,
               save=f'_{cmap_n}_{gn}.pdf')
    print(_g, ':', _vmax)

##vmax = 1.6
# sc.pl.umap(adata, color=genes_plot[:],
##           vmax=vmax, 
#           ncols=4,
#           cmap = cmap_g,
#           save=f'_{cmap_n}_{genes_plot}T.pdf')


# In[] 
''' just for testing some gene un UMAP (not run) '''
_tmp_g = 'RGS5'  # 'EMCN'
sc.pl.umap(adata, color=_tmp_g,
           #           vmax=vmax,
           ncols=2,
           cmap=cmap_g,
           save=f'_{cmap_g}_{_tmp_g}.pdf')

# In[]
''' correlation matrix
'''
n_top = 20
genes_on = B.TopMarkers(adata, n_top)
groupby = 'leiden_anno'
corr_mtd = 'correlation'

gmeans = B.GroupMean(adata, groupby, features=genes_on, use_raw=True)

# zscore
_gmeans_z = B.Zscore(gmeans.T).T
gmeans_z = pd.DataFrame(_gmeans_z, index=gmeans.index, columns=gmeans.columns)
sims = gmeans_z.corr()
# gmeans.to_csv(resdir / f'gmeans_top{n_top}mk.csv', index=True, header=True)
# gmeans_z.to_csv(resdir / f'gmeans_z_top{n_top}mk.csv', index=True, header=True)
# sims.to_csv(resdir / f'sims_z_top{n_top}mk.csv', index=True, header=True)

# sims = B.GroupSimilarityMean(adata, groupby, use_genes=genes_on,
#                        metric=corr_mtd, 
#                        use_raw=True, tags=None,)
# In[]
''' correlation matrix (Scanpy build-in) 
'''
# import seaborn as sns
cmap_corr = ['RdBu_r', 'vlag', 'magma_r'][-1]
sc.set_figure_params(dpi_save=150, fontsize=14, )

cluster_grid = sns.clustermap(
    sims, cmap=cmap_corr,  # ax=ax,
    method=['single', 'average'][1],
    metric="correlation",
    xticklabels=False,
    linewidths=.5,
    square=True,
    cbar_pos=None,
    figsize=(15, 10),
    row_colors=cl_color_match,
    row_linkage=None, col_linkage=None, )
cluster_grid.ax_col_dendrogram.remove()

# plt.savefig(figdir / f'corr_{cmap_corr}_top{n_top}mk-1122-1.pdf', bbox_inches='tight')

# rowids = cluster_grid.dendrogram_row.reordered_ind
# type_ordered = sims.index.take(rowids)
# fx.save_namelist(type_ordered, resdir / f'celltype_ordered-{n_top}.csv')

# In[]
''' heatmap
'''
genes_hmap = fx.load_namelist(fxd.GENEDIR / 'genes_heatmap.txt')
avgs_major = B.GroupMean(adata, key_major, features=genes_hmap, use_raw=True)
tag = ['_with16', ''][-1]

if tag == '_with16':
    ordered_major_types_use = ordered_major_types
else:
    ordered_major_types_use = [c for c in ordered_major_types if
                               c not in ['Unknown']]
avgs_major = avgs_major[ordered_major_types_use]  # order columns
avgs_major.columns

df_hmap = avgs_major.apply(lambda x: x / x.max(), axis=1)
df_hmap

df_hmap1 = fx.order_contingency_df(df_hmap, axis=1)
df_hmap1.to_csv(resdir / f'mean_expressions_maxNormed_ordered-1125{tag}.csv',
                index_label='gene')
# df_hmap = B.Zscore(df_hmap)

# In[]
sc.set_figure_params(dpi_save=180, fontsize=11, )
cmap_heat = ['RdBu_r', 'magma_r'][0]

fig, ax = plt.subplots(figsize=(4, 12))
sns.heatmap(df_hmap1,
            yticklabels=False,
            # cbar=False,
            cbar_kws={'shrink': 0.3},
            ax=ax,
            cmap=cmap_heat,
            # col_colors = cl_color_match,
            )
fig.savefig(figdir / f'heatmap_plain{tag}-{cmap_heat}-1125.pdf',
            bbox_inches='tight')

# In[]
sc.set_figure_params(dpi_save=180, fontsize=11, )
col_colors = pd.DataFrame(
    {'major class': major_colors_srs[ordered_major_types_use], },
    index=df_hmap1.columns)

heatmap_grid = sns.clustermap(
    df_hmap1, figsize=(4, 12),
    cmap=cmap_heat,
    col_cluster=False,
    row_cluster=False,
    #        metric = "cosine",
    col_colors=col_colors,
    #        vmax = 2.5,
    yticklabels=False,
    cbar_pos=None,
    #        cbar_kws = {'cbar': False}
)
heatmap_grid.ax_col_dendrogram.remove()
heatmap_grid.ax_row_dendrogram.remove()
plt.show()
heatmap_grid.savefig(
    figdir / f'heatmap_anno{tag}-{cmap_heat}-1125.pdf', bbox_inches='tight')

# col_colors.to_csv(resdir / f'col_colors_heatmap.csv', )


# In[]
''' order of the clusters in the heatmap'''
dct0 = {'Glutamatergic neurons': [0, 17, 11, 14, 3, 5, 6, 7, 10, 13, 15, 20],
        'GABAergic neurons': [4, 8, 12, 21],
        'Oligodendrocytes': [9, 16, 2, 18],
        'Ependymal cells': [19, 22],
        'Endothelial cells': [23],
        'Microglials': [24],
        'Astrocytes': [1],
        'Cajal-Retzius': [25]}

candi_colors = fx.get_colors('tab20b', n=len(dct0))
fx.view_color_map(candi_colors)

# ordered_clusters = []
main_type_colors = []
ordered_ids = []
i = 0
for nm, ids in dct0.items():
    # ordered_clusters += [f'{i} {nm}' for i in ids]
    ordered_ids += ids
    main_type_colors += [candi_colors[i]] * len(ids)
    i += 1

leiden_anno
# ordered_clusters = leiden_anno[ordered_ids]
#
# print(ordered_clusters)
# fx.save_namelist(ordered_clusters, fxd.HippDirFormal / 'cluster_order_heatmap.txt')


# dotplot
_n = 2
degs = B.TopMarkers(adata, _n, )
sc.set_figure_params(dpi_save=150, fontsize=12, )

vmax_dot = 3
ax = sc.pl.rank_genes_groups_dotplot(
    adata, n_genes=_n, vmax=vmax_dot,
    figsize=(6 * _n + 0.5, 6),
    save=f'_markers_top{_n}.pdf')

# In[]
''' UMAP highlight only lineage
'''
tt1 = 'putative trajectory'
lin_gr0 = ['2', '9', '16', ]
adata.obs['isLineage'] = adata.obs['leiden'].apply(lambda x: x in lin_gr0)
sc.pl.umap(adata, color='leiden', groups=lin_gr0,
           legend_loc='on data',
           title=tt1,
           save='_highlight_lin_ondata')

sc.pl.umap(adata, color='isLineage',  # groups = [True],
           title=tt1,
           palette=['lightgrey', '#ce1256'],
           save='_highlight_lin1')

lin_gr = leiden_anno[list(map(int, lin_gr0))].values

# In[]
''' Take out the lineage and save (the normalized data)
'''

dir_formal_lin = fxd.DATADIR / '_formal_traj'

# adata_l = sc.read_h5ad(dir_formal_lin / 'lineage0.h5ad')
# adata_l.obs['leiden_anno'] = adata.obs['leiden_anno'][adata_l.obs_names].values
# sc.pl.umap(adata_l, color='leiden_anno', groups=lin_gr, save='_0_lin', title='re-embedded trajectory')
#
## re-analysis
# sc.tl.leiden(adata_l, resolution=0.5)
# sc.pl.umap(adata_l, color='leiden', legend_loc='on data', save='_lin')
#
# sc.tl.rank_genes_groups(adata_l, 'leiden')
# sc.pl.rank_genes_groups_dotplot(adata_l, n_genes=5, )
# B.marker_table(adata).to_csv(resdir / 'lineage_markers.csv', index=False)
##sc.pl.diffmap(adata_l, color='leiden', legend_loc='on data', save='_lin')
##sc.pl.pca(adata_l, color='leiden', legend_loc='on data', save='_lin')
##adata_l.write(dir_formal_lin / 'lineage0.h5ad')
##fx.saveNamedMtx(adata_l, dir_formal_lin / 'lineage_norm_mtx')
##B.save_embeddings(adata_l, dir_formal_lin, 'X_umap', tail='lineage0')
# adata_l.obs.to_csv(dir_formal_lin / 'metadata_lineage.csv', index=True)
# genes_candi = B.TopMarkers(adata_l, 50,)
# fx.save_namelist(genes_candi, dir_formal_lin / 'candidate_genes.csv')
# adata_l

###################################################
''' or directly load from `formal dir`
=========================================================
'''
dir_formal_lin = fxd.DATADIR_NEW / 'formal-traj'
adata_l = sc.read_h5ad(dir_formal_lin / 'lineage0.h5ad')

key_time = 'sling_pseudotime'

sudo_time0 = pd.read_csv(dir_formal_lin / 'pseudotime_slingshot.csv',
                         header=None, )
sudo_time0.head()
sudo_time = sudo_time0.iloc[:, 1].values
sudo_time = sudo_time / max(sudo_time)  # scale to (0, 1)
adata_l.obs['dpt_pseudotime'] = sudo_time

adata_l.obs[key_time] = sudo_time
adata_l.obs[key_time].head()

cmap_lin = [None, 'Spectral_r', 'magma_r'][-1]
frame = False
sc.pl.umap(adata_l, color=key_time, title='pseudotime',
           cmap=cmap_lin, size=25,
           frameon=frame,
           save=f'_time_lineage_{cmap_lin}_{frame}')

# In[]

path_name = '245310'
path_ord = list(path_name)
top_n_gene = 6
genes_lin = B.TopMarkers(adata_l, top_n_gene, groups=path_ord)
use_origin = True

if use_origin:
    sc.tl.paga(adata_l, 'leiden_anno')
    path_name = 'orig'
    path_ord = np.take(lin_gr, [1, 2, 0])

figsize = (8, 7)
cmap_hm_lin = ['magma_r', 'viridis'][-1]
fig, axs = plt.subplots(figsize=figsize, dpi=150)
sc.pl.paga_path(adata_l,
                nodes=path_ord,
                keys=genes_lin,
                ax=axs,
                n_avg=50,
                save='%s_top%d_%s.pdf' % (path_name, top_n_gene, cmap_hm_lin),
                show_node_names=False,
                ytick_fontsize=9,
                legend_fontsize=8,
                color_map=cmap_hm_lin,
                )
