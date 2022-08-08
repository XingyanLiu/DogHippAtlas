# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 12:13:50 2019

@author: xyliu
"""

from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import pyplot as plt
from _submit.scripts.src import build as B

'''
                * Settings
==========================================
'''

DATADIR = Path(r'E:\Users\xyliu\data\hippocampus')
DATADIR_NEW = Path(r'E:\Users\xyliu\data003\dog')
HippDir0 = DATADIR / 'All' / 'concat'
HippDir = DATADIR_NEW / 'hipp_gHVGs_ppca_50'
HippDirFormal = DATADIR_NEW / 'formal-hipp'

GENEDIR = DATADIR_NEW / 'formal-genes'

NSDir = DATADIR / 'NS' / 'concat'

Primers = ('random', 'polyT')
_name_to_primer = {'random': 'EF', 'polyT': 'AB'}

# ==========================================
#         1: For cell type annotations
# ==========================================

GenesDomastic = 'AGAP1 ANKS1B CCBE1 \
    CUX2 ENSCAFG00000022711 GAD2\
    MBP OPCML PCSK5  ROR1'.split()


'''
==========================================
        0-0: Getting datasets
==========================================
'''


def DataHippFormal():
    datadir = HippDirFormal
    fname = 'analyzed.h5ad'
    return _load_data(datadir, fname)


def DataHipp(sub='ppca', datadir=HippDir, fname=None):
    '''
    sub: either be 'ppca' or 'plssvd'
    if `datadir` and `fname` are not NoneType, the `sub` will be ignored.

    '''
    if datadir is None:
        datadir = HippDir0 / f'Merged_gHVGs_{sub}_50'
    fname = 'Merged_cosine_1.h5ad' if fname is None else fname
    return _load_data(datadir, fname)


def DataHippLineage(sub='ppca'):
    datadir = HippDir / f'Analysis_{sub}'
    fname = 'lineage_9-16-2.h5ad'
    return _load_data(datadir, fname)


def RawDataHipp():
    datadir = HippDirFormal
    fname = 'afterQC.h5ad'
    return _load_data(datadir, fname)


def GeneCorr(corr='spearman', detailed=False, bulk=True):
    datadir = NSDir / 'GeneStageCorr_bulk'
    fname = f'corr_{corr}_qc150_filtered.csv'
    df = _load_data(datadir, fname, pd.read_csv, index_col=0)
    df = df.sort_values('correlation', )
    print(df.head())
    if detailed:
        return df
    else:
        return df.index.tolist()


def _load_data(datadir, fname, _method=sc.read_h5ad, **kwds):
    fpath = datadir / fname
    print(f'Loading data from file:\n {fpath}')
    adata = _method(fpath, **kwds)
    return adata


'''
==========================================
        2: Lineage analysis
==========================================
Step 0:
    Separate groups
    
Step 1:
    Re-embedding: neighbors -> UMAP
    Re-clustering
    Purify (optional)
    DE (marker finding)
    
Step 2:
    
    pseudotime inference


'''


def PH_analysis(ph, n_pcs=50, nneigh=10, metric='cosine',
                min_dist=0.25, res=0.6, save=True):
    ph.Neighbors(nneigh=nneigh, n_pcs=n_pcs, metric=metric)
    ph.Umap(min_dist=min_dist, plot_color='leiden')
    ph.Clustering(res=res, keep_results=True)
    ph.DE(plot=True, save=save)  # , method='wilcoxon')


def PH_pseudotime(ph, root_group, groupby='leiden', diff_comps=15, n_dcs=10, ):
    sc.tl.diffmap(ph.adata, n_comps=diff_comps)
    sc.pl.diffmap(ph.adata, color=groupby, components=['1,2', '3,4'])

    ph.adata.uns['iroot'] = group_center(ph.adata, root_group, groupby)
    sc.tl.dpt(ph.adata, n_dcs=n_dcs)
    PH_vis_group_dpt(ph, groupby='leiden')


#    ph.vis(color = [groupby, 'dpt_pseudotime'],
#           legend_loc='on data',
#           save=f'_{groupby}_dpt_{ph.name}.png')

def PH_paga(ph, groups='leiden'):
    sc.tl.paga(ph.adata, groups=groups)
    sc.pl.paga(ph.adata)
    # ph.Umap(init_pos='paga')
    sc.pl.paga_compare(ph.adata, save=f'_{ph.name}.png')


def PH_dynamic(ph, path_ord, top_n_gene=5,
               figsize=(8, 6), save_genes=True, saveph=True, **kw):
    if isinstance(path_ord, str):
        path_name = path_ord
        path_ord = list(path_ord)
    elif isinstance(path_ord, list):
        path_name = '-'.join(path_ord)
    else:
        print('`path_ord` should be either string or a list!')

    genes = ph.markers_top(top_n_gene, groups=path_ord)
    if save_genes:
        pd.Series(genes).to_csv(
            ph.resdir / ('dynamic_genes_path%s.csv' % path_ord),
            header=False, index=False)
    #    len(genes)
    fig, axs = plt.subplots(figsize=figsize, dpi=100)
    sc.pl.paga_path(ph.adata,
                    nodes=list(path_ord),
                    keys=genes,
                    ax=axs,
                    n_avg=50,
                    save='%s_top%d.png' % (path_name, top_n_gene),
                    #                    save = '%s_top%d_avg50.png'%(path_ord, n),
                    ytick_fontsize=8,
                    legend_fontsize=8,
                    #                    annotations=['dpt_pseudotime'],
                    **kw)
    if saveph: ph.save()


def PH_vis_group_dpt(ph, groupby='leiden'):
    ph.vis(color=[groupby, 'dpt_pseudotime'],
           legend_loc='on data',
           save=f'_{groupby}_dpt_{ph.name}.png')


def PH_vis_paga_compare(ph):
    sc.pl.paga_compare(ph.adata, save=f'_{ph.name}.png')


def PH_vis_dpt(ph, ):
    ph.vis(color='dpt_pseudotime',
           legend_loc='on data',
           save=f'_dpt_{ph.name}.png')


def group_center(adata, group, groupby='leiden'):
    """
    find the center of a group, return the index
    temperally take the node with max degree.
    """
    #    indices_ = ph.adata.obs[groupby]  == group
    indices = np.flatnonzero(adata.obs[groupby] == group)
    A = adata.uns['neighbors']['connectivities'][indices, :]

    # compute the degrees
    d = A.sum(axis=1).A1
    center_id = indices[np.argmax(d)]
    print(f'Auto-selected node {center_id}, with max degree {d.max()}')
    return center_id


def dynamic_paga(adata, path_ord, top_n_gene=5, genes=None,
                 figsize=(8, 6), **kwds):
    if isinstance(path_ord, str):
        path_name = path_ord
        path_ord = list(path_ord)
    elif isinstance(path_ord, list):
        path_name = '-'.join(path_ord)
    else:
        print('`path_ord` should be either string or a list!')

    genes = B.TopMarkers(adata, top_n_gene,
                         groups=path_ord) if genes is None else genes
    #    len(genes)
    fig, axs = plt.subplots(figsize=figsize, dpi=100)
    sc.pl.paga_path(adata,
                    nodes=list(path_ord),
                    keys=genes,
                    ax=axs,
                    n_avg=50,
                    save='%s_top%d.png' % (path_name, top_n_gene),
                    #                    save = '%s_top%d_avg50.png'%(path_ord, n),
                    ytick_fontsize=8,
                    legend_fontsize=8,
                    **kwds)


''' 
================================================================
                            Plottings
================================================================
'''


def umap_dpt_wrapper(adata, figsize=None, tt=None, save=True):
    tt = 'Diffusion pseudotime' if tt is None else tt
    figsize = (4, 3.5) if figsize is None else figsize
    fig, ax = plt.subplots(figsize=figsize)
    sc.pl.umap(adata, color='dpt_pseudotime', ax=ax,  # frameon=False,
               title=tt, show=False,
               legend_loc='on data', )
    d = _relocate_ax_horiz(ax, 2)
    _save_fig_or_nothing(ax.figure, save=save)

    return ax


##################################################################
# Other helper functions
def _save_fig_or_nothing(fig, save=True):
    if save:
        fname = save if isinstance(save, (str, Path)) else 'temp_figure.png'
        print(f'figure was saved into file:\n\t{fname}')
        fig.savefig(fname, bbox_inches='tight')


def _relocate_ax_horiz(ax, nd=1, d=None, right=False, ):
    """
    """
    pos0 = ax.get_position()
    x, y, w, h = pos0.x0, pos0.y0, pos0.width, pos0.height
    d = x / 10 if d is None else d
    if right: d *= -1  # move left by default
    x = x - nd * d
    pos1 = [x, y, w, h]
    ax.set_position(pos1)
    #    if squre: ax.set_aspect('equal', adjustable='datalim')
    return d
