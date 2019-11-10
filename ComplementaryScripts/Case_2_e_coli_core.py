#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-09-12

"""Step1_iML1515_test.py
:description : script
:param : 
:returns: 
:rtype: 
"""
import os

import ConvexHull_yield
import cobra
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull

os.chdir('../ComplementaryData/')

# %% <method EFMs > can not calculate
print('\n---------- Loading EFMs ... ---------- ')
core_EFMs_z = np.genfromtxt('Case2_1_ecoli_core/e_coli_core_EFMs_standardized.csv', delimiter = ',')
EFMs_all_points = core_EFMs_z.T
EFMs_all_points = EFMs_all_points[:,1:]   # modes x reas
experiment_datas = []
print('EFMs:', EFMs_all_points.shape)

# %% <method EFVs > can not calculate
print('\n---------- Loading EFVs ... ---------- ')
core_EFVs_z = np.genfromtxt('Case2_1_ecoli_core/e_coli_core_EFVs_standardized.csv', delimiter = ',')
EFVs_all_points = core_EFVs_z.T
EFVs_all_points = EFVs_all_points[:,1:]   # modes x reas
experiment_datas = []
print('EFVs:', EFVs_all_points.shape)

#%% <FBA mode>
print('\n---------- Loading FBA modes ... ---------- ')
core_FBAMs_z = np.genfromtxt('Case2_1_ecoli_core/ecoli_core_FBAMs_standardized.csv', delimiter = ',')

FBAMs_all_points = core_FBAMs_z.T  # modes x rea
FBAMs_all_points = FBAMs_all_points[:, 1:]
experiment_datas = []
print('FBA models:', FBAMs_all_points.shape)

#%% <Our method>
print('\n---------- Caculating Yield space by this method ... ---------- ')
e_coli_core = cobra.io.read_sbml_model('../ComplementaryData/Case2_1_ecoli_core/e_coli_core.xml')
#e_coli_core.reactions.get_by_id('EX_o2_e').bounds = (0.0,1000.0)
model = e_coli_core.copy()

production_rea_ids_x = ['BIOMASS_Ecoli_core_w_GAM', 'EX_ac_e', 'EX_for_e', 'EX_etoh_e', 'EX_lac__D_e', 'EX_succ_e', ]
production_rea_ids_y = ['BIOMASS_Ecoli_core_w_GAM', 'EX_ac_e', 'EX_for_e', 'EX_etoh_e', 'EX_lac__D_e', 'EX_succ_e', ]
carbon_source_rea_id = 'EX_glc__D_e'

steps = 10
carbon_uptake_direction = -1
model.reactions.get_by_id(carbon_source_rea_id).bounds = (-10, -0.001)

# all modes
# all modes
# yield_normalized_df, fluxes_all, hull_index_all = GEM2pathways.get_yield_space_multi(model, production_rea_ids_x,
#                                                                                      production_rea_ids_y,
#                                                                                      carbon_source_rea_id,
#                                                                                      steps=steps,
#                                                                                      carbon_uptake_direction=carbon_uptake_direction,
#                                                                                      draw=True)
# yield_normalized_df.to_csv('Case2_1_ecoli_core/ecoli_core_our_yield_normalized_df.csv', index=True, sep=',')
# yield_normalized_df_ = pd.read_csv('Case2_1_ecoli_core/ecoli_core_our_yield_normalized_df.csv', sep=',',
#                                    index_col=0, )
# yield_normalized_df_hull = yield_normalized_df_[yield_normalized_df_.columns[hull_index_all]]
# yield_normalized_df_hull.to_csv('Case2_1_ecoli_core/ecoli_core_our_yield_normalized_df_hull.csv', index=True,
#                                 sep=',')
yield_normalized_df_hull_ = pd.read_csv('Case2_1_ecoli_core/ecoli_core_our_yield_normalized_df_hull.csv', sep=',',
                                        index_col=0, )

our_all_points = yield_normalized_df_hull_.values.T
our_all_points = our_all_points[ abs(our_all_points[:,0]) > 1e-10,:]
our_all_points = our_all_points[:,1:]       #exclude the carbon colume
print('Our method:', our_all_points.shape)

# %% <MYA >
# experiment data

# %%
experiment_datas = []
qhull_options = 'QJ Qx A0.9999999'  # 'Qt QJ Pp Qw Qx'      #'QG0'mean expect point 0(index) TODO Qx will lose 7 points check it
cutoff_persent = 1  # can't reduce much points because too many dimasions....
# Fixme check the hull_all ??? 169points??? to many !!! why ??? Decimals??
# Note not Decimals problem try other reason
# decimals = 3
# EFMs_all_points_ = np.around(EFMs_all_points,decimals,)
# FBAMs_all_points_ = np.around(FBAMs_all_points,decimals,)
# our_all_points_ = np.around(our_all_points,decimals,)

# EFMs_indexes, EFMs_weights, EFMs_estimated_datas, EFMs_in_hulls = ConvexHull_yield.pipeline_mya(EFMs_all_points,
#                                                                                                 experiment_datas,
#                                                                                                 qhull_options=qhull_options,
#                                                                                                 method=1,
#                                                                                                 cutoff_persent=cutoff_persent)
#
# EFVs_indexes, EFVs_weights, EFVs_estimated_datas, EFVs_in_hulls = ConvexHull_yield.pipeline_mya(EFVs_all_points,
#                                                                                                 experiment_datas,
#                                                                                                 qhull_options=qhull_options,
#                                                                                                 method=1,
#                                                                                                 cutoff_persent=cutoff_persent)

# FBAMs_indexes, FBAMs_weights, FBAMs_estimated_datas, FBAMs_in_hulls = ConvexHull_yield.pipeline_mya(FBAMs_all_points,
#                                                                                                     experiment_datas,
#                                                                                                     qhull_options=qhull_options,
#                                                                                                     method=1,
#                                                                                                     cutoff_persent=cutoff_persent)

our_indexes, our_weights, our_estimated_datas, our_in_hulls = ConvexHull_yield.pipeline_mya(our_all_points,
                                                                                            experiment_datas,
                                                                                            qhull_options=qhull_options,
                                                                                            method=1,
                                                                                            cutoff_persent=cutoff_persent)


#%% <plot initial yield>: TODO not divided by carbon, so not yield

# yield_rea_ids_name = yield_rea_ids
yield_rea_ids_name = ['Acetate', 'Formate', 'Ethanol', 'Lacate', 'Succinate']


def trim_axs(axs, N):
    """little helper to massage the axs list to have correct length..."""
    axs = axs.flat
    for ax in axs[N:]:
        ax.remove()
    return axs[:N]

cols = 3
rows = len(yield_rea_ids_name) // cols + 1
figsize = (10, 8)
fig, axs = plt.subplots(cols, rows,figsize=figsize)
axs = trim_axs(axs, len(yield_rea_ids_name))

for ax , exmet_reaction in zip(axs,yield_rea_ids_name):
    index = yield_rea_ids_name.index(exmet_reaction)+1

    xy_EFMs = EFMs_all_points[:, [0,index]]
    xy_EFVs = EFVs_all_points[:, [0,index]]
    xy_FBAMs = FBAMs_all_points[:, [0, index]]
    xy_our = our_all_points[:, [0,index]]

    colors_list = ['blue','teal','tab:red','tab:orange']

    points_EFMs = ax.plot(xy_EFMs[:,0], xy_EFMs[:,1],       '.', markerfacecolor='none', color=colors_list[0],
                    alpha = 0.5, markersize=3)
    points_EFVs = ax.plot(xy_EFVs[:,0], xy_EFVs[:,1],       '.', markerfacecolor='none', color=colors_list[1],
                    alpha = 0.5, markersize=3)

    points_our = ax.plot(xy_our[:,0], xy_our[:,1],          '^', markerfacecolor='none', color=colors_list[2],
                         alpha=0.8, markersize=6)
    points_FBAMs = ax.plot(xy_FBAMs[:, 0], xy_FBAMs[:, 1], 's', color=colors_list[3],
                           alpha=1, markersize=8)

    draw_hull = True

    if draw_hull:
        hull_EFMs = ConvexHull(xy_EFMs,qhull_options = qhull_options )
        hull_EFVs = ConvexHull(xy_EFVs, qhull_options=qhull_options)
        hull_FBAMs = ConvexHull(xy_FBAMs, qhull_options=qhull_options)
        hull_our = ConvexHull(xy_our,qhull_options = qhull_options )

        for simplex in hull_EFMs.simplices:
            line_EFMs = ax.plot(xy_EFMs[simplex, 0], xy_EFMs[simplex, 1], 'x-', markerfacecolor='none', color=colors_list[0],
                            alpha= 0.5 ,  markersize=3,lw=2)

        for simplex in hull_EFVs.simplices:
            line_EFVs = ax.plot(xy_EFVs[simplex, 0], xy_EFVs[simplex, 1], 'v-', markerfacecolor='none', color=colors_list[1],
                            alpha = 0.5, markersize=3 ,lw=2)

        for simplex in hull_our.simplices:
            line_our = ax.plot(xy_our[simplex, 0], xy_our[simplex, 1], '^-',markerfacecolor='none', color=colors_list[2],
                            alpha = 0.5, markersize=5,lw=2)

        for simplex in hull_FBAMs.simplices:
            line_FBAMs = ax.plot(xy_FBAMs[simplex, 0], xy_FBAMs[simplex, 1], 'o-', markerfacecolor='none',
                                 color=colors_list[3],
                                 alpha=0.5, markersize=5, lw=2)

    else:
        line_EFMs = points_EFMs
        line_EFVs = points_EFVs
        line_our = points_our
        line_FBAMs = points_FBAMs

    ax.set_ylabel(yield_rea_ids_name[index-1]+' / Glc',fontsize = 18)

ax.set_xlabel('Yield Biomass/Glucose', fontsize=18)
fig.legend((line_EFMs[0], line_EFVs[0], line_our[0], line_FBAMs[0]), ('EFMs', 'EFVs', 'This study', 'FBA modes'),
           bbox_to_anchor=(0.6, 0.3), loc='upper left', borderaxespad=0., prop={'size': 18})
fig.show()

fig.savefig('Case2_1_ecoli_core/4.png',format ='png')


# %% ac

yield_rea_ids_name = ['Acetate', 'Formate', 'Ethanol', 'Lacate', 'Succinate']
yield_rea_ids_name = ['Acetate']

fig = plt.figure()
ax = fig.add_subplot(111)

cols = 1
rows = len(yield_rea_ids_name) // cols + 1
figsize = (10, 8)
# ax= fig.subplots(111)

index = 1


xy_EFMs = EFMs_all_points[:, [0,index]]
xy_EFVs = EFVs_all_points[:, [0,index]]
xy_FBAMs = FBAMs_all_points[:, [0, index]]
xy_our = our_all_points[:, [0,index]]

colors_list = ['blue','teal','tab:red','tab:orange']

points_EFMs = ax.plot(xy_EFMs[:,0], xy_EFMs[:,1],       '.', markerfacecolor='none', color=colors_list[0],
                alpha = 0.5, markersize=3)
points_EFVs = ax.plot(xy_EFVs[:,0], xy_EFVs[:,1],       '.', markerfacecolor='none', color=colors_list[1],
                alpha = 0.5, markersize=3)

points_our = ax.plot(xy_our[:,0], xy_our[:,1],          '^', markerfacecolor='none', color=colors_list[2],
                     alpha=0.8, markersize=6)
points_FBAMs = ax.plot(xy_FBAMs[:, 0], xy_FBAMs[:, 1], 's', color=colors_list[3],
                       alpha=1, markersize=8)

draw_hull = True

if draw_hull:
    hull_EFMs = ConvexHull(xy_EFMs,qhull_options = qhull_options )
    hull_EFVs = ConvexHull(xy_EFVs, qhull_options=qhull_options)
    hull_FBAMs = ConvexHull(xy_FBAMs, qhull_options=qhull_options)
    hull_our = ConvexHull(xy_our,qhull_options = qhull_options )

    for simplex in hull_EFMs.simplices:
        line_EFMs = ax.plot(xy_EFMs[simplex, 0], xy_EFMs[simplex, 1], 'x-', markerfacecolor='none', color=colors_list[0],
                        alpha= 0.5 ,  markersize=3,lw=2)

    for simplex in hull_EFVs.simplices:
        line_EFVs = ax.plot(xy_EFVs[simplex, 0], xy_EFVs[simplex, 1], 'v-', markerfacecolor='none', color=colors_list[1],
                        alpha = 0.5, markersize=3 ,lw=2)

    for simplex in hull_our.simplices:
        line_our = ax.plot(xy_our[simplex, 0], xy_our[simplex, 1], '^-',markerfacecolor='none', color=colors_list[2],
                        alpha = 0.5, markersize=5,lw=2)

    for simplex in hull_FBAMs.simplices:
        line_FBAMs = ax.plot(xy_FBAMs[simplex, 0], xy_FBAMs[simplex, 1], 'o-', markerfacecolor='none',
                             color=colors_list[3],
                        alpha = 0.5, markersize=10,lw=2)

else:
    line_EFMs = points_EFMs
    line_EFVs = points_EFVs
    line_our = points_our
    line_FBAMs = points_FBAMs

ax.set_ylabel(yield_rea_ids_name[index-1]+' / Glc',fontsize = 18)

ax.set_xlabel('Yield Biomass/Glucose', fontsize=18)
fig.legend((line_EFMs[0], line_EFVs[0], line_our[0], line_FBAMs[0]), ('EFMs', 'EFVs', 'This study', 'FBA modes'),
           bbox_to_anchor=(0.6, 0.95), loc='upper left', borderaxespad=0., prop={'size': 16})
fig.show()

fig.savefig('Case2_1_ecoli_core/0_ac.png',format ='png')