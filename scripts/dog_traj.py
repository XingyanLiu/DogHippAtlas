# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 12:12:51 2019

@author: xyliu
"""

import os
from pathlib import Path
import scanpy as sc

import matplotlib as mpl
#matplotlib.use('TkAgg')
print(mpl.get_backend()) # module://ipykernel.pylab.backend_inline
print(mpl.matplotlib_fname())

WORKDIR = Path(r'D:\Users\xyliu\003')
os.chdir(WORKDIR)

from src.pipe_handler import PipeHandler
from _submit.scripts.src import funx as fx, funx_dog as fxd
from src import build as B

'''     Settings
==========================
'''
sub = 'ppca'
resdir = fxd.HippDir / f'Analysis_{sub}'
figdir = resdir / 'figs'

sc.settings.figdir = figdir
sc.set_figure_params(dpi=100, fontsize=8,)

''' Loading dataset
'''


adata = fxd.DataHipp('ppca')

#sc.pl.umap(adata, color='leiden', legend_loc='on data')

''' Cell type annotaions
==========================
'''


adata.obs['leiden'].cat.categories
adata.obs['leiden_anno'] = adata.obs['leiden']
adata.obs['leiden_anno'].cat.categories = [
        '0', '1', '2 Oligo', '3 Glu Sst4_14', '4 GABA', '5 Glu neuron', 
        '6 Glu Sst4_14', '7 Glu neuron', '8 GABA', '9 Progenitor', '10', '11', '12',
       '13', '14', '15', '16 Oligo', '17 neuron', '18', '19', '20', '21', '22', '23', '24',
       '25'
        ]
sc.pl.umap(adata, color='leiden_anno', legend_loc='on data', save=f'_{sub}_ondata.png')
sc.pl.umap(adata, color='leiden_anno', save=f'_{sub}.png', )

# cluster 9 is the Progenitor cell

''' 
            Lineage Analysis
=============================================
9 --> 16
'''

groups0 = ['9', '16', '2']
_tail = '-'.join(groups0)

## have a look for validation
tt = 'Putative lineage from the hippocampus'
sc.pl.umap(adata, color='leiden', groups=groups0, title=tt,
           legend_loc='on data', save=f'_fromAll_{_tail}.png')


adt0 = B.TakeGroups(adata, groups0, 'leiden', copy=True)
ph0 = PipeHandler(adt0, f'lineage_{_tail}', resdir=resdir)



## re-embedding
ph = ph0
fxd.PH_analysis(ph, n_pcs = 50, nneigh=10, metric='cosine', 
                min_dist=0.25, res=0.6, save=True)

B.MergeGroups(ph.adata, 'leiden', ['0', '1'], new_key='leiden')

# purify: seems not necessary to re-analyze after puritication!
adt00 = B.RemoveGroups(ph.adata, list('678'), 'leiden')
ph1 = PipeHandler(adt00, f'lineage_{_tail}', resdir=resdir)
# seems not necessary to re-analyze after puritication!
#fxd.PH_analysis(ph, n_pcs = 30, nneigh=10, metric='cosine', 
#                min_dist=0.25, res=0.6, save=True)
ph1.DE(save=True)

'''     PAGA and Pseudotime inference 
'''
fxd.PH_paga(ph1, 'leiden')
#fxd.group_center(ph.adata, '3')

fxd.PH_pseudotime(ph, '3')

## just for visualization, not necessary
fxd.PH_vis_group_dpt(ph1)
fxd.PH_vis_paga_compare(ph1)
fxd.PH_vis_dpt(ph1)
ph1.save()

'''expression dynamic
'''

n = 8
path_ord = list('3452') + ['0_1']
fxd.PH_dynamic(ph1, path_ord, top_n_gene=n)

''' 9 --> 5
'''

groups1 = ['9', '5']
adt1 = B.TakeGroups(adata, groups1, 'leiden', copy=True)


fx.close_figs()


#############################################################################
## GENE plot on Hippocampus UMAP-Embedding
corr = 'spearman'
ns_genes = fxd.GeneCorr(corr)
ns_genes

for i, g in enumerate(ns_genes):
    nb = str(i + 1)
    if len(nb) < 2: nb = '0' + nb
    save = f'_{corr}_gene{nb}({g}).png'
    sc.pl.umap(adata, color=g, save=save)




