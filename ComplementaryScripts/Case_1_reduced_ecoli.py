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
import cobra
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull

os.chdir('../ComplementaryData/')

'''
model:

-METEXT
GLU SUC FOR ACT LAC ETH B CO2 H2

R1 : GLU + PEP => G6P + PYR .
R2 : G6P + ATP => 2 T3P + ADP .
R3 : G6P + 6 NAD => T3P + 6 NADH .
R4 : T3P + NAD + ADP => PEP + NADH + ATP .
R5 : PEP + ADP => PYR + ATP .
R6 : PEP + CO2 + 2 NADH => SUC + 2 NAD .
R7 : PYR + CoA => AcCoA + FOR .
R8 : PYR + NADH => LAC + NAD .
R9 : AcCoA + ADP => ACT + CoA + ATP .
R10 : AcCoA + 2 NADH => ETH + CoA + 2 NAD .
R11 : FOR => CO2 + H2 .
R12 : 6.775 G6P + 82.2 ATP + 4.065 NADH => B + 82.2 ADP + 4.065 NAD .
'''


#%% <EFMs > caculated by matlab emftool. and standardized by step0
print('\n---------- Loading EFMs ... ---------- ')
em_z = np.genfromtxt('Case1_ecoli_reduced/ecoli_reduced_EFMs_standardized.csv', delimiter = ',')
EFMs_all_points = em_z.T
EFMs_all_points = EFMs_all_points[:,1:]   # modes x reas
print('EFMs:', EFMs_all_points.shape)


#%% <FBA modes> from reference and standardized by step0
print('\n---------- Loading FBA modes ... ---------- ')
FBAMs_z = np.genfromtxt('Case1_ecoli_reduced/ecoli_reduced_FBAMs_standardized.csv', delimiter = ',')
FBAMs_all_points = FBAMs_z.T  # modes x rea
FBAMs_all_points = FBAMs_all_points[:, 1:]
print('FBA models:', FBAMs_all_points.shape)


#%% <Our method>
print('\n---------- Caculating Yield space by this method ... ---------- ')
# load a GEM
ecoli_reduced_model = cobra.io.read_sbml_model('Case1_ecoli_reduced/ecoli_reduced_model.xml')
model = ecoli_reduced_model.copy()

production_rea_ids_x = ['R12', 'R9', 'EX_FOR', 'R10', 'R8', 'R6']
production_rea_ids_y = ['R12', 'R9', 'EX_FOR', 'R10', 'R8', 'R6']
carbon_source_rea_id = 'R1'
model.reactions.get_by_id(carbon_source_rea_id).bounds = (0.001, 10)  # set carbon lb as -10
steps = 10
carbon_uptake_direction = 1

# all modes
yield_normalized_df, fluxes_all, hull_index_all = GEM2pathways.get_yield_space_multi(model, production_rea_ids_x,
                                                                                     production_rea_ids_y,
                                                                                     carbon_source_rea_id,
                                                                                     steps=steps,
                                                                                     carbon_uptake_direction=carbon_uptake_direction,
                                                                                     draw=True)

yield_normalized_df.to_csv('Case1_ecoli_reduced/ecoli_reduced_our_yield_normalized_df.csv', index=True, sep=',')
yield_normalized_df_ = pd.read_csv('Case1_ecoli_reduced/ecoli_reduced_our_yield_normalized_df.csv', sep=',',
                                   index_col=0, )
yield_normalized_df_hull = yield_normalized_df_[yield_normalized_df_.columns[hull_index_all]]
yield_normalized_df_hull.to_csv('Case1_ecoli_reduced/ecoli_reduced_our_yield_normalized_df_hull.csv', index=True,
                                sep=',')
yield_normalized_df_hull_ = pd.read_csv('Case1_ecoli_reduced/ecoli_reduced_our_yield_normalized_df_hull.csv', sep=',',
                                        index_col=0, )
our_all_points = yield_normalized_df_hull_.values.T
our_all_points = our_all_points[ abs(our_all_points[:,0]) > 1e-10,:]
our_all_points = our_all_points[:, 1:]  # exclude the carbon colume!!!
print('Our method:', our_all_points.shape)

# %% <MYA >
print('\n---------- Loading experiment data ... ---------- ')
experiment_data_df = pd.read_csv('Case1_ecoli_reduced/ecoli_reduce_experiment_data.txt', delimiter='\t', header=0)
data_cloum_name = ['glc', 'biomass', 'ac', 'for', 'etoh', 'lac', 'succ']
experiment_data_df_trimed = experiment_data_df[data_cloum_name]
experiment_data_df_trimed_values = experiment_data_df_trimed.values[:, :] - experiment_data_df_trimed.values[0, :]
experiment_data_df_trimed_values = experiment_data_df_trimed_values[1:, :]
experiment_data_df_trimed_values = experiment_data_df_trimed_values.T
experiment_data_df_trimed_values = experiment_data_df_trimed_values[:, :] / abs(experiment_data_df_trimed_values[0, :])
experiment_data_df_trimed_values = experiment_data_df_trimed_values.T
experiment_data_df_trimed_values = experiment_data_df_trimed_values[:, 1:]

experiment_datas = [ ]      # TODO experiment data!!!
for i in range(0, experiment_data_df_trimed_values.shape[0]):
    experiment_datas.append(experiment_data_df_trimed_values[i, :])
# %%
qhull_options = 'QJ Qx A0.9999999'  # 'Qt QJ Pp Qw Qx'      #'QG0'mean expect point 0(index) TODO Qx will lose 7 points check it
cutoff_persent = 1  # can't reduce much points because too many dimasions....
# Fixme check the hull_all ??? 169points??? to many !!! why ??? Decimals??
# Note not Decimals problem try other reason
# decimals = 3
# EFMs_all_points_ = np.around(EFMs_all_points,decimals,)
# FBAMs_all_points_ = np.around(FBAMs_all_points,decimals,)
# our_all_points_ = np.around(our_all_points,decimals,)

EFMs_indexes, EFMs_weights, EFMs_estimated_datas, EFMs_in_hulls = ConvexHull_yield.pipeline_mya(EFMs_all_points,
                                                                                                experiment_datas,
                                                                                                qhull_options=qhull_options,
                                                                                                method=1,
                                                                                                cutoff_persent=cutoff_persent)
# # FBAMs not have enouth point to perform hull
# # FBAMs_indexes, FBAMs_weights, FBAMs_estimated_datas, FBAMs_in_hulls = ConvexHull_yield.pipeline_mya(FBAMs_all_points,
# #                                                                                             experiment_datas,
# #                                                                                             qhull_options=qhull_options,
# #                                                                                             method=1,
# #                                                                                             cutoff_persent=cutoff_persent)
#

# experiment_datas = [['', 0.70959994, 0.65218536, 0.75269061, 0.00598727, 0.10675876]]
our_indexes, our_weights, our_estimated_datas, our_in_hulls = ConvexHull_yield.pipeline_mya(our_all_points,
                                                                                            experiment_datas,
                                                                                            qhull_options=qhull_options,
                                                                                            method=1,
                                                                                            cutoff_persent=cutoff_persent)

# EFMs_hull = ConvexHull(EFMs_all_points,qhull_options = qhull_options )
# # FBAMs_hull = ConvexHull(FBAMs_all_points,qhull_options = qhull_options )
our_hull = ConvexHull(our_all_points[:, 1:], qhull_options=qhull_options)
ConvexHull_yield.point_in_hull(experiment_datas[-1][1:], our_hull, tolerance=1e-12)


#%% <plot initial yield>: TODO not divided by carbon, so not yield
cmap_list = ['Pastel1', 'Pastel2', 'Paired', 'Accent',
             'Dark2', 'Set1', 'Set2', 'Set3',
             'tab10', 'tab20', 'tab20b', 'tab20c']

plt.get_cmap('Pastel1')

yield_rea_ids_name = ['EX_ac_e','EX_for_e','EX_etoh_e','EX_lac__D_e','EX_succ_e',]


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

experiment_points = np.array(experiment_datas)

for ax , exmet_reaction in zip(axs,yield_rea_ids_name):
    index = yield_rea_ids_name.index(exmet_reaction)+1
    
    xy_EFMs = EFMs_all_points[:, [0,index]]
    xy_FBAMs = FBAMs_all_points[:, [0, index]]
    xy_our = our_all_points[:, [0,index]]
    xy_exp = experiment_points[-4:-1, [0, index]]

    hull_EFMs = ConvexHull(xy_EFMs, qhull_options='QJ')
    hull_FBAMs = ConvexHull(xy_FBAMs, qhull_options='QJ')
    hull_our = ConvexHull(xy_our, qhull_options='QJ')

    points_our = ax.plot(xy_our[:,0], xy_our[:,1],                 '^',markerfacecolor='none',color = 'tab:blue',
                      alpha = 1,label = 'This_methd',markersize=5)
    for simplex in hull_our.simplices: 
        line_our = ax.plot(xy_our[simplex, 0], xy_our[simplex, 1], '^--',markerfacecolor='none',color = 'tab:blue',
                      alpha = 0.5,label = 'This_methd',markersize=5)

    points_FBAMs = ax.plot(xy_FBAMs[:, 0], xy_FBAMs[:, 1], 'o', markerfacecolor='none', color='black',
                      alpha=1, label='FBA_mode', markersize=10)
    for simplex in hull_FBAMs.simplices:
        line_FBAMs = ax.plot(xy_FBAMs[simplex, 0], xy_FBAMs[simplex, 1], 'o--', markerfacecolor='none', color='black',
                      alpha=0.5, label='FBA_mode', markersize=10)

    points_EFMs = ax.plot(xy_EFMs[:, 0], xy_EFMs[:, 1], 'x', markerfacecolor='none', color='tab:orange',
                          alpha=1, label='EFMs', markersize=10)
    for simplex in hull_EFMs.simplices:
        line_EFMs = ax.plot(xy_EFMs[simplex, 0], xy_EFMs[simplex, 1], 'x-', markerfacecolor='none', color='tab:orange',
                            alpha=0.5, label='EFMs', markersize=10)

    points_exp = ax.plot(xy_exp[:, 0], xy_exp[:, 1], '^', color='tab:red',
                         alpha=1, label='experiment data', markersize=11)

    ax.set_ylabel(yield_rea_ids_name[index-1]+'/Glucose',fontsize = 12)

ax.set_xlabel('Yield Biomass/Glucose',fontsize = 12)
fig.legend((line_EFMs[0], line_FBAMs[0], line_our[0]), ('EFMs', 'FBA modes', 'This study'), bbox_to_anchor=(0.55, 0.25),
           loc='upper left', borderaxespad=0.)
fig.show()
fig.show()
