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
import GEM2pathways
import Cybernetic_Functions
import cobra
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull

os.chdir('../Data/')

# %% <method EFMs > from efmtools Matlab
print('\n---------- Loading EFMs ... ---------- ')
core_EFMs_z = np.genfromtxt('Case2_1_ecoli_core/e_coli_core_EFMs_standardized.csv', delimiter=',')
EFMs_all_points = core_EFMs_z.T
EFMs_all_points = EFMs_all_points[:, 1:]  # modes x reas
experiment_datas = []
print('EFMs:', EFMs_all_points.shape)

# %% <method EFVs >  by Matbab
print('\n---------- Loading EFVs ... ---------- ')
core_EFVs_z = np.genfromtxt('Case2_1_ecoli_core/e_coli_core_EFVs_standardized.csv', delimiter=',')
EFVs_all_points = core_EFVs_z.T
EFVs_all_points = EFVs_all_points[:, 1:]  # modes x reas
experiment_datas = []
print('EFVs:', EFVs_all_points.shape)

# %% <FBA mode>
print('\n---------- Loading FBA modes ... ---------- ')
core_FBAMs_z = np.genfromtxt('Case2_1_ecoli_core/ecoli_core_FBAMs_standardized.csv', delimiter=',')

FBAMs_all_points = core_FBAMs_z.T  # modes x rea
FBAMs_all_points = FBAMs_all_points[:, 1:]
experiment_datas = []
print('FBA models:', FBAMs_all_points.shape)

# %% <Our method>
print('\n---------- Caculating Yield space by this method ... ---------- ')
e_coli_core = cobra.io.read_sbml_model('../Data/Case2_1_ecoli_core/e_coli_core.xml')
# e_coli_core.reactions.get_by_id('EX_o2_e').bounds = (0.0,1000.0)
# glc as sub
model = e_coli_core.copy()

production_rea_ids_x = ['BIOMASS_Ecoli_core_w_GAM', ]
production_rea_ids_y = ['EX_ac_e', 'EX_for_e', 'EX_etoh_e', 'EX_lac__D_e', 'EX_succ_e', ]
carbon_source_rea_id = 'EX_glc__D_e'

steps = 10
carbon_uptake_direction = -1
model.reactions.get_by_id(carbon_source_rea_id).bounds = (-10, -0.001)

# all modes
yield_normalized_df, fluxes_all, hull_index_all = GEM2pathways.get_yield_space_multi(model, production_rea_ids_x,
                                                                                     production_rea_ids_y,
                                                                                     carbon_source_rea_id,
                                                                                     steps=steps,
                                                                                     carbon_uptake_direction=carbon_uptake_direction,
                                                                                     draw=True)
# yield_normalized_df.to_csv('Case2_1_ecoli_core/ecoli_core_our_yield_normalized_df.csv', index=True, sep=',')
yield_normalized_df_ = pd.read_csv('Case2_1_ecoli_core/ecoli_core_our_yield_normalized_df.csv', sep=',',
                                   index_col=0, )
yield_normalized_df_hull = yield_normalized_df_[yield_normalized_df_.columns[hull_index_all]]
# yield_normalized_df_hull.to_csv('Case2_1_ecoli_core/ecoli_core_our_yield_normalized_df_hull.csv', index=True,
#                                 sep=',')
yield_normalized_df_hull_ = pd.read_csv('Case2_1_ecoli_core/ecoli_core_our_yield_normalized_df_hull.csv', sep=',',
                                        index_col=0, )

our_all_points = yield_normalized_df_hull_.values.T
our_all_points = our_all_points[abs(our_all_points[:, 0]) > 1e-10, :]
our_all_points = our_all_points[:, 1:]  # exclude the carbon colume
print('Our method:', our_all_points.shape)

# ac as sub: Note limtied by experiment data
model = e_coli_core.copy()

model.reactions.get_by_id('EX_glc__D_e').bounds = (0, 1000)
model.reactions.get_by_id('EX_ac_e').bounds = (-10, -0.001)
p1_yield_value_max, fluxes_opt_max = GEM2pathways.get_yield_opt(model, 'BIOMASS_Ecoli_core_w_GAM', 'EX_ac_e',
                                                                max_or_min='max',
                                                                carbon_uptake_direction=carbon_uptake_direction)
p1_yield_value_min, fluxes_opt_min = GEM2pathways.get_yield_opt(model, 'BIOMASS_Ecoli_core_w_GAM', 'EX_ac_e',
                                                                max_or_min='min',
                                                                carbon_uptake_direction=-1)

additional_points1 = fluxes_opt_max[['EX_glc__D_e', 'BIOMASS_Ecoli_core_w_GAM', 'EX_ac_e']].values
additional_points2 = fluxes_opt_min[['EX_glc__D_e', 'BIOMASS_Ecoli_core_w_GAM', 'EX_ac_e']].values
additional_points = np.insert(additional_points1, 3, additional_points2, axis=0)
additional_points = additional_points.reshape((2, 3)).T
additional_points = additional_points[:, :] / abs(additional_points[-1, :])

our_all_points = yield_normalized_df_hull_.values.T
our_all_points = our_all_points[abs(our_all_points[:, 0]) > 1e-10, :]
our_all_points = our_all_points[:, 1:]  # exclude the carbon colume
print('Our method:', our_all_points.shape)

# %% <MYA >

print('\n---------- Loading experiment data ... ---------- ')
experiment_data_df = pd.read_csv('Case2_1_ecoli_core/Ecoli_aerobic_experiment_data.txt', delimiter='\t', header=0)

data_cloum_name = ['glc', 'biomass', 'ac']
experiment_data_df_trimed = experiment_data_df[data_cloum_name]
experiment_data_df_trimed_values = experiment_data_df_trimed.values[:, :] - experiment_data_df_trimed.values[0, :]
experiment_data_df_trimed_values = experiment_data_df_trimed_values[1:, :]
experiment_data_df_trimed_values = experiment_data_df_trimed_values.T
experiment_data_df_trimed_values = experiment_data_df_trimed_values[:, :] / abs(experiment_data_df_trimed_values[0, :])
experiment_data_df_trimed_values = experiment_data_df_trimed_values.T
experiment_data_df_trimed_values = experiment_data_df_trimed_values[:, 1:]

experiment_datas = []  # TODO experiment data!!!
for i in range(0, experiment_data_df_trimed_values.shape[0]):
    experiment_datas.append(list(experiment_data_df_trimed_values[i, :]) + [''] * 4)

qhull_options = 'QJ Qx A0.9999999'  # 'Qt QJ Pp Qw Qx'      #'QG0'mean expect point 0(index)
cutoff_persent = 1  # can't reduce much points because too many dimasions....
# Fixed (Fixme check the hull_all ??? 169 points??? to many !!! why ??? Decimals??)
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

# %% fig 4a

figsize = (3.25, 3)

fig, ax = plt.subplots(figsize=figsize)

index = 1
xy_EFMs = EFMs_all_points[:, [0, index]]
xy_EFVs = EFVs_all_points[:, [0, index]]
# fig test
# xy_EFMs = xy_EFMs[0:1000, :]
# xy_EFVs = xy_EFVs[0:1000, :]
our_all_points = yield_normalized_df_.values.T
our_all_points = our_all_points[:, 1:]

xy_our = our_all_points[:, [0, index]]
colors_list = ['blue', 'teal', 'tab:orange', 'tab:red', ]

points_EFMs = ax.plot(xy_EFMs[:, 0], xy_EFMs[:, 1], '.', markerfacecolor='none', color=colors_list[0],
                      alpha=0.3, markersize=2)
points_EFVs = ax.plot(xy_EFVs[:, 0], xy_EFVs[:, 1], '.', markerfacecolor='none', color=colors_list[1],
                      alpha=0.3, markersize=2)

points_our = ax.plot(xy_our[:, 0], xy_our[:, 1], '+', markerfacecolor='none', color=colors_list[2],
                     alpha=1, markersize=5)

hull_EFMs = ConvexHull(xy_EFMs, qhull_options=qhull_options)
hull_EFVs = ConvexHull(xy_EFVs, qhull_options=qhull_options)
hull_our = ConvexHull(xy_our, qhull_options=qhull_options)

for simplex in hull_EFMs.simplices:
    line_EFMs = ax.plot(xy_EFMs[simplex, 0], xy_EFMs[simplex, 1], 'x-', markerfacecolor='none',
                        color=colors_list[0],
                        alpha=0.5, markersize=5, lw=1)

for simplex in hull_EFVs.simplices:
    line_EFVs = ax.plot(xy_EFVs[simplex, 0], xy_EFVs[simplex, 1], '*-', markerfacecolor='none',
                        color=colors_list[1],
                        alpha=0.8, markersize=6, lw=2)

for simplex in hull_our.simplices:
    line_our = ax.plot(xy_our[simplex, 0], xy_our[simplex, 1], '+-', markerfacecolor='none',
                       color=colors_list[2],
                       alpha=0.8, markersize=7, lw=1)

# ax.tick_params(axis="x", labelsize=8, pad=1)
# ax.tick_params(axis="y", labelsize=8, pad=1),

# ax.set_xlim((-0.01, 0.1))
ax.set_ylim((-0.1, 2.2))
ax.set_yticks(np.arange(0, 2.1, 1))
ax.set_yticklabels([0, 1, 2, ])
plt.yticks(fontname="Arial", fontsize=8)
plt.xticks(fontname="Arial", fontsize=8)
ax.tick_params(axis="x", labelsize=8, pad=4)
ax.tick_params(axis="y", labelsize=8, pad=1),
ax.set_ylabel('Yield: Acetate/Glucose', fontsize=10, family='Arial', labelpad=1)
ax.set_xlabel('Yield: Biomass/Glucose', fontsize=10, family='Arial', labelpad=4)

# fig.legend((line_EFMs[0], line_EFVs[0], line_our[0],), ('EFMs', 'EFVs', 'This study'),
#            bbox_to_anchor=(0.6, 0.3), loc='upper left', borderaxespad=0., prop={'family': 'Arial', 'size': 8})
#

fig.legend((line_EFMs[0], line_EFVs[0], line_our[0],),
           ('EFMs', 'EFVs', 'This study'),

           # bbox_to_anchor=(0.65, -0.12),ncol=2, loc=8,
           loc=4, bbox_to_anchor=(0.91, 0.7),
           # loc=2, bbox_to_anchor=(0.15, 0.93),
           prop={'family': 'Arial', 'size': 8})

fig.tight_layout()
# fig.savefig('Case2_1_ecoli_core/fig4a_ac_biomass.pdf', bbox_inches='tight')
# fig.savefig('Case2_1_ecoli_core/fig4a_ac_biomass.svg', bbox_inches='tight')
fig.savefig('Case2_1_ecoli_core/fig4a_ac_biomass.jpg', bbox_inches='tight', dpi=600)

fig.show()

yield_rea_ids_name = ['Acetate', 'Formate', 'Ethanol', 'Lacate', 'Succinate']  # ,
yield_rea_ids_name = ['Ethanol', 'Lactate', 'Succinate', 'Formate', ]  # ,
yield_rea_ids_lables = ['Ethanol', 'Lactate', 'Succinate', 'Formate', ]


def trim_axs(axs, N):
    """little helper to massage the axs list to have correct length..."""
    axs = axs.flat
    for ax in axs[N:]:
        ax.remove()
    return axs[:N]


cols = 2
rows = 2
figsize = (4.5, 3.5)
fig, axs = plt.subplots(cols, rows, figsize=figsize)
axs = trim_axs(axs, len(yield_rea_ids_name))

for ax, exmet_reaction in zip(axs, yield_rea_ids_name):
    index = yield_rea_ids_name.index(exmet_reaction) + 1
    lable = yield_rea_ids_lables[index - 1]
    if lable == 'Ethanol':
        index = 3
    if lable == 'Acetate':
        index = 1
    elif lable == 'Lactate':
        index = 4
    elif lable == 'Succinate':
        index = 5
    elif lable == 'Formate':
        index = 2

    xy_EFMs = EFMs_all_points[:, [0, index]]
    xy_EFVs = EFVs_all_points[:, [0, index]]
    # fig test
    # xy_EFMs = xy_EFMs[0:1000, :]
    # xy_EFVs = xy_EFVs[0:1000, :]
    our_all_points = yield_normalized_df_.values.T
    our_all_points = our_all_points[:, 1:]

    xy_our = our_all_points[:, [0, index]]

    colors_list = ['blue', 'teal', 'tab:orange', 'tab:red', ]

    points_EFMs = ax.plot(xy_EFMs[:, 0], xy_EFMs[:, 1], '.', markerfacecolor='none', color=colors_list[0],
                          alpha=0.3, markersize=1)
    points_EFVs = ax.plot(xy_EFVs[:, 0], xy_EFVs[:, 1], '.', markerfacecolor='none', color=colors_list[1],
                          alpha=0.3, markersize=1)

    points_our = ax.plot(xy_our[:, 0], xy_our[:, 1], '+', markerfacecolor='none', color=colors_list[2],
                         alpha=1, markersize=3)

    draw_hull = True

    if draw_hull:
        hull_EFMs = ConvexHull(xy_EFMs, qhull_options=qhull_options)
        hull_EFVs = ConvexHull(xy_EFVs, qhull_options=qhull_options)
        hull_our = ConvexHull(xy_our, qhull_options=qhull_options)

        for simplex in hull_EFMs.simplices:
            line_EFMs = ax.plot(xy_EFMs[simplex, 0], xy_EFMs[simplex, 1], 'x-', markerfacecolor='none',
                                color=colors_list[0],
                                alpha=0.5, markersize=3, lw=1)

        for simplex in hull_EFVs.simplices:
            line_EFVs = ax.plot(xy_EFVs[simplex, 0], xy_EFVs[simplex, 1], '*-', markerfacecolor='none',
                                color=colors_list[1],
                                alpha=0.8, markersize=5, lw=2)

        for simplex in hull_our.simplices:
            line_our = ax.plot(xy_our[simplex, 0], xy_our[simplex, 1], '+-', markerfacecolor='none',
                               color=colors_list[2],
                               alpha=0.8, markersize=6, lw=1)

    else:
        line_EFMs = points_EFMs
        line_EFVs = points_EFVs
        line_our = points_our
        line_FBAMs = points_FBAMs

    # ax.set_xlim((-0.01, 0.1))
    ax.set_ylim((-0.2, 2.2))
    ax.set_yticks(np.arange(0, 2.1, 1))
    ax.set_yticklabels([0, 1, 2, ])
    if index == 2:
        ax.set_ylim((-0.2, 4.4))
        ax.set_yticks(np.arange(0, 4.1, 1))
        ax.set_yticklabels([0, 1, 2, 3, 4])
        # ax.set_yticks(np.arange(0, 1.2, 0.2))
        # ax.set_yticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.tick_params(axis="x", labelsize=8, pad=1)
    ax.tick_params(axis="y", labelsize=8, pad=1),

    ax.set_ylabel(lable + '/Glucose', fontsize=10, family='Arial', labelpad=1)

# ax.set_xlabel('Yield Biomass/Glucose', fontsize=10, family='Arial')
# fig.legend((line_EFMs[0], line_EFVs[0], line_our[0],), ('EFMs', 'EFVs', 'This study'),
#            bbox_to_anchor=(0.6, 0.3), loc='upper left', borderaxespad=0., prop={'family': 'Arial', 'size': 8})

# fig.legend((line_EFMs[0], line_EFVs[0], line_our[0], ), ('EFMs', 'EFVs', 'This study',),
#            ('All points', 'Key points', 'EFMs', 'ConvexHull'),
#            bbox_to_anchor=(0.8, -0.22),
#            loc=8, ncol=1, prop={'family': 'Arial', 'size': 8})


fig.tight_layout(pad=1, h_pad=0.5, w_pad=0.5)
# fig.savefig('Case2_1_ecoli_core/fig4a_all_mets_biomass.pdf', bbox_inches='tight')
fig.savefig('Case2_1_ecoli_core/fig4a_all_mets_biomass.jpg', bbox_inches='tight', dpi=600)

# fig.savefig('Case2_1_ecoli_core/fig4a_all_mets_biomass.svg', bbox_inches='tight', transparent=True)

fig.show()

# %% <plot initial yield>: TODO not divided by carbon, so not yield


yield_rea_ids_name = ['Acetate', 'Formate', 'Ethanol', 'Lacate', 'Succinate']  # ,


def trim_axs(axs, N):
    """little helper to massage the axs list to have correct length..."""
    axs = axs.flat
    for ax in axs[N:]:
        ax.remove()
    return axs[:N]


cols = 3
rows = len(yield_rea_ids_name) // cols + 1
figsize = (10, 8)
fig, axs = plt.subplots(cols, rows, figsize=figsize)
axs = trim_axs(axs, len(yield_rea_ids_name))

for ax, exmet_reaction in zip(axs, yield_rea_ids_name):
    index = yield_rea_ids_name.index(exmet_reaction) + 1

    xy_EFMs = EFMs_all_points[:, [0, index]]
    xy_EFVs = EFVs_all_points[:, [0, index]]
    xy_FBAMs = FBAMs_all_points[:, [0, index]]
    xy_our = our_all_points[:, [0, index]]

    colors_list = ['blue', 'teal', 'tab:red', 'tab:orange']

    points_EFMs = ax.plot(xy_EFMs[:, 0], xy_EFMs[:, 1], '.', markerfacecolor='none', color=colors_list[0],
                          alpha=0.5, markersize=3)
    points_EFVs = ax.plot(xy_EFVs[:, 0], xy_EFVs[:, 1], '.', markerfacecolor='none', color=colors_list[1],
                          alpha=0.5, markersize=3)

    points_our = ax.plot(xy_our[:, 0], xy_our[:, 1], '^', markerfacecolor='none', color=colors_list[2],
                         alpha=0.8, markersize=6)
    points_FBAMs = ax.plot(xy_FBAMs[:, 0], xy_FBAMs[:, 1], 's', color=colors_list[3],
                           alpha=1, markersize=8)

    draw_hull = True

    if draw_hull:
        hull_EFMs = ConvexHull(xy_EFMs, qhull_options=qhull_options)
        hull_EFVs = ConvexHull(xy_EFVs, qhull_options=qhull_options)
        hull_FBAMs = ConvexHull(xy_FBAMs, qhull_options=qhull_options)
        hull_our = ConvexHull(xy_our, qhull_options=qhull_options)

        for simplex in hull_EFMs.simplices:
            line_EFMs = ax.plot(xy_EFMs[simplex, 0], xy_EFMs[simplex, 1], 'x-', markerfacecolor='none',
                                color=colors_list[0],
                                alpha=0.5, markersize=3, lw=2)

        for simplex in hull_EFVs.simplices:
            line_EFVs = ax.plot(xy_EFVs[simplex, 0], xy_EFVs[simplex, 1], 'v-', markerfacecolor='none',
                                color=colors_list[1],
                                alpha=0.5, markersize=3, lw=2)

        for simplex in hull_our.simplices:
            line_our = ax.plot(xy_our[simplex, 0], xy_our[simplex, 1], '^-', markerfacecolor='none',
                               color=colors_list[2],
                               alpha=0.5, markersize=5, lw=2)

        for simplex in hull_FBAMs.simplices:
            line_FBAMs = ax.plot(xy_FBAMs[simplex, 0], xy_FBAMs[simplex, 1], 'o-', markerfacecolor='none',
                                 color=colors_list[3],
                                 alpha=0.5, markersize=5, lw=2)

    else:
        line_EFMs = points_EFMs
        line_EFVs = points_EFVs
        line_our = points_our
        line_FBAMs = points_FBAMs

    ax.set_ylabel(yield_rea_ids_name[index - 1] + ' / Glc', fontsize=18)

ax.set_xlabel('Yield Biomass/Glucose', fontsize=18)
fig.legend((line_EFMs[0], line_EFVs[0], line_our[0], line_FBAMs[0]), ('EFMs', 'EFVs', 'This study', 'FBA modes'),
           bbox_to_anchor=(0.6, 0.3), loc='upper left', borderaxespad=0., prop={'size': 18})
fig.show()

# fig.savefig('Case2_1_ecoli_core/4.png',format ='png')


# %% plot ac

yield_rea_ids_name = ['Acetate', 'Formate', 'Ethanol', 'Lacate', 'Succinate']
yield_rea_ids_name = ['Acetate']

fig = plt.figure()
ax = fig.add_subplot(111)

cols = 1
rows = len(yield_rea_ids_name) // cols + 1
figsize = (10, 8)
# ax= fig.subplots(111)

index = 1

xy_EFMs = EFMs_all_points[:, [0, index]]
xy_EFVs = EFVs_all_points[:, [0, index]]
xy_FBAMs = FBAMs_all_points[:, [0, index]]
xy_our = our_all_points[:, [0, index]]

colors_list = ['blue', 'teal', 'tab:red', 'tab:orange']

points_EFMs = ax.plot(xy_EFMs[:, 0], xy_EFMs[:, 1], '.', markerfacecolor='none', color=colors_list[0],
                      alpha=0.5, markersize=3)
points_EFVs = ax.plot(xy_EFVs[:, 0], xy_EFVs[:, 1], '.', markerfacecolor='none', color=colors_list[1],
                      alpha=0.5, markersize=3)

points_our = ax.plot(xy_our[:, 0], xy_our[:, 1], '^', markerfacecolor='none', color=colors_list[2],
                     alpha=0.8, markersize=6)
# points_FBAMs = ax.plot(xy_FBAMs[:, 0], xy_FBAMs[:, 1], 's', color=colors_list[3],
#                        alpha=1, markersize=8)

draw_hull = True

if draw_hull:
    hull_EFMs = ConvexHull(xy_EFMs, qhull_options=qhull_options)
    hull_EFVs = ConvexHull(xy_EFVs, qhull_options=qhull_options)
    hull_FBAMs = ConvexHull(xy_FBAMs, qhull_options=qhull_options)
    hull_our = ConvexHull(xy_our, qhull_options=qhull_options)

    for simplex in hull_EFMs.simplices:
        line_EFMs = ax.plot(xy_EFMs[simplex, 0], xy_EFMs[simplex, 1], 'x-', markerfacecolor='none',
                            color=colors_list[0],
                            alpha=0.5, markersize=3, lw=2)

    for simplex in hull_EFVs.simplices:
        line_EFVs = ax.plot(xy_EFVs[simplex, 0], xy_EFVs[simplex, 1], 'v-', markerfacecolor='none',
                            color=colors_list[1],
                            alpha=0.5, markersize=3, lw=2)

    for simplex in hull_our.simplices:
        line_our = ax.plot(xy_our[simplex, 0], xy_our[simplex, 1], '^-', markerfacecolor='none', color=colors_list[2],
                           alpha=0.5, markersize=5, lw=2)

    # for simplex in hull_FBAMs.simplices:
    #     line_FBAMs = ax.plot(xy_FBAMs[simplex, 0], xy_FBAMs[simplex, 1], 'o-', markerfacecolor='none',
    #                          color=colors_list[3],
    #                     alpha = 0.5, markersize=10,lw=2)

else:
    line_EFMs = points_EFMs
    line_EFVs = points_EFVs
    line_our = points_our
    # line_FBAMs = points_FBAMs

ax.set_ylabel(yield_rea_ids_name[index - 1] + ' / Glc', fontsize=18)

ax.set_xlabel('Yield Biomass/Glucose', fontsize=18)
fig.legend((line_EFMs[0], line_EFVs[0], line_our[0], line_FBAMs[0]), ('EFMs', 'EFVs', 'This study', 'FBA modes'),
           bbox_to_anchor=(0.6, 0.95), loc='upper left', borderaxespad=0., prop={'size': 16})
fig.legend((line_EFMs[0], line_EFVs[0], line_our[0]), ('EFMs', 'EFVs', 'This study',),
           bbox_to_anchor=(0.6, 0.95), loc='upper left', borderaxespad=0., prop={'size': 16})
fig.show()

# fig.savefig('Case2_1_ecoli_core/0_ac.png',format ='png')


# %%  <Cybernetic model simulations>:


tStart = 0.0  # DefineTime
tStop = 10
tStep = 0.1
# tspan = np.linspace(tStart, tStop, (tStop - tStart) / tStep)
tspan = np.linspace(tStart, tStop, int(((tStop - tStart) / tStep) + 1))

# matrix Z: reactions x pathways; and Smz metabolites x pathways , S metabolites x reactions : Smz = Sm @ Z
final_index = our_indexes[-1][-1]
Smz = yield_normalized_df_hull_.values[:, final_index]
Smz = Smz[[0, 1, 2], :]
# Smz[0, :] = -Smz[0, :]
# path_ac = np.array([0, 0, -1])  # Note !!! this pathway is nessary  to simulate the for experimentdata
# Smz = np.column_stack((Smz,path_for))
Smz = np.insert(Smz, 3, values=additional_points.T, axis=1)
# metabolites and pathways number
(n_mets, n_path) = Smz.shape

# experiment data
metabObj = ['glc', 'biomass', 'ac']  # , 'for', 'etoh', 'lac', 'succ'

# initial metabolites at tome 0, t0: initial_mets.shape = (n_mets,)
initial_mets = experiment_data_df[metabObj].values[0, :]
# initial_mets = np.insert(initial_mets, 3, values=np.array([0.01,0.01,0.01,0.01]), axis=0)
# initial:Enzyme: initial_enzyme.shape = (n_path,)
initial_enzyme = np.array([0.9] * n_path)

# initial data x0 :initial_x0.shape  = (n_mets + n_path,)
initial_x0 = np.concatenate((initial_mets, initial_enzyme))

# Enzyme Rate Parameters: alpha,beta,ke : de/dt =  alpha + rE(ke) * u - (beta + mu) * e
alpha = np.array([0.04] * n_path)
beta = np.array([0.05] * n_path)
ke = np.array([0.620342] * n_path)  # or 0.5

# Metabolites rate Parameters kmax , Ki : dm/dt =  Smz @ V @ rM(kmax,K) * c
# kmax : n_path
kmax = np.array([

    # 1.6641e+02,
    # 1.6801e+01,
    # 2.7765e+01,
    # 7.3644e+00,
    # 4.0808e-02,
    163.04025717699113, 18.431402942510104, 30.353025531495078, 10.835179478216705, 0.029443545394048323,
])

# K : n_path
K = np.array([
    # 3.2575e-03,
    # 2.4931e+00,
    # 3.5281e+00,
    # 5.1913e-03,
    # 3.8855e-03,
    0.0014519693251603275, 1.2204425112563406, 6.387152105506235, 0.003588519284764188, 0.004125375519668958
])

# carbon number for each pathways
n_carbon = np.array([6, 6, 6, 2, 2])
Biomass_index = 1
sub_index = 0


class Cybernetic_Model_basic(Cybernetic_Functions.Cybernetic_Model):
    pass


CB_model = Cybernetic_Model_basic('CB model for Ecoli core matrix')

# CB_model = Cybernetic_Functions.Cybernetic_Model('CB model for Ecoli core matrix ')
CB_model.Smz = Smz
CB_model.x0 = initial_x0
CB_model.kmax = kmax
CB_model.K = K
CB_model.ke = ke
CB_model.alpha = alpha
CB_model.beta = beta
CB_model.n_carbon = n_carbon
CB_model['sub_index'] = sub_index
CB_model['Biomass_index'] = Biomass_index
CB_model.mets_name = metabObj
CB_model = CB_model


def rate_def(self, x):
    name = self.name
    Smz = self.Smz
    kmax = self.kmax
    K = self.K
    ke = self.ke
    Biomass_index = self.Biomass_index
    sub_index = self.sub_index

    (n_mets, n_path) = Smz.shape

    r_kin_basic = [
        kmax[0] * x[sub_index] / (K[0] + x[sub_index]),
        kmax[1] * x[sub_index] / (K[1] + x[sub_index]),
        kmax[2] * x[sub_index] / (K[2] + x[sub_index]),
        kmax[3] * x[2] / (K[3] + x[2]),
        kmax[4] * x[2] / (K[4] + x[2])]

    rM = r_kin_basic * x[n_mets:]

    rE = ke * r_kin_basic / kmax

    rG = Smz[Biomass_index, :] * rM[:]  # 11: biomass index

    return rM, rE, rG


setattr(Cybernetic_Model_basic, 'rate_def', rate_def)

sol = Cybernetic_Functions.cb_model_simulate(CB_model, tspan, draw=True)

weight_of_method_t_f = 3
weights_of_mets = np.array([1, 5, 3, ]) ** 2

CB_model.weight_of_method_t_f = False
# CB_model.weights_of_mets = weights_of_mets  # /1000
experiment_data_dfs = [experiment_data_df]  # [5],[6],[7],[8]]
para_to_fit = {'kmax': np.array([[0], [1], [2], [3], [4]]),
               'K': np.array([[0], [1], [2], [3], [4]])}

# minimum = Cybernetic_Functions.parameters_fitting([CB_model], experiment_data_dfs, para_to_fit, tspan, draw=True,
#                                                   full_output=False,
#                                                   options={'xtol': 0.01, 'ftol': 0.01, 'maxiter': 200, 'disp': True})
#
# CB_model = Cybernetic_Functions.update_paras_func(minimum.x, CB_model, para_to_fit,
#                                                   retern_model_or_paras='model')[0]

sol = Cybernetic_Functions.cb_model_simulate(CB_model, tspan, draw=True)

np.savetxt('Case2_1_ecoli_core/ecoli_core_sol.csv', sol, delimiter=',')
sol_ = np.loadtxt('Case2_1_ecoli_core/ecoli_core_sol.csv', delimiter=',')

# %% <plot cybernetic model result>ac

# experiment data


fig = plt.figure(figsize=(5.5, 2.5))
ax = fig.add_subplot(111)
# colors = ['blue', 'teal', 'tab:red', 'tab:orange']
color_list = plt.cm.tab10(np.linspace(0, 1, 11))

for index in range(0, CB_model.Smz.shape[0]):
    if index == 1:
        ax1 = ax.twinx()
        ax1.plot(tspan, sol[:, index], color=color_list[index], linewidth=1, label=metabObj[index])

        experiment_p = ax1.plot(experiment_data_df['time'], experiment_data_df[metabObj[index]], '*',
                                color=color_list[index], alpha=0.8,
                                linewidth=1)
    else:
        ax.plot(tspan, sol[:, index], color=color_list[index], linewidth=1, label=metabObj[index])

        experiment_p = ax.plot(experiment_data_df['time'], experiment_data_df[metabObj[index]], '*',
                               color=color_list[index], alpha=0.8,
                               linewidth=1)

ax.set_xlabel('Time (h)', fontsize=10, family='Arial', )
ax.set_ylabel('Concentration (mM)', fontsize=10, family='Arial', )
ax1.set_ylabel('Biomass (g/L)', fontsize=10, family='Arial', )
plt.yticks(fontname="Arial", fontsize=8)
plt.xticks(fontname="Arial", fontsize=8)
L = fig.legend(loc='lower left', bbox_to_anchor=(1, 0.2), ncol=1, fontsize=8, )
plt.setp(L.texts, family='Arial')
fig.tight_layout()
# fig.savefig('Case2_1_ecoli_core/fig4c_core_simulate.pdf', bbox_inches='tight')
# fig.savefig('Case2_1_ecoli_core/fig4c_core_simulate.jpg', bbox_inches='tight', dpi=600)
fig.show()
