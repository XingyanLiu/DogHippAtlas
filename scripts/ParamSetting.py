# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 22:36:18 2019

@author: xyliu
"""


_parameters0 = dict(
    All=dict(qc=True, min_genes=350, max_genes=5000,
             rmv_mito=True, mito_perc=0.02,
             counts_level=None,
             plot_hvgs=True, min_mean=0.035, min_disp=0.25,
             do_regress=False, batch_key=None,
             n_comps=100,
             metric='cosine',
             nneigh=30, n_pcs=60,
             min_dist=0.5,
             de=True, plot_de=True,
             cluster=True,
             save_middle=True,
             save_by_default=True
             ),
)


def get_parameters(key):
    return _parameters0[key].copy()

