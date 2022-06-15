#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 6/24/20

"""Case_3_test.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import ConvexHull_yield
import GEM2pathways
import Cybernetic_Functions
import cobra
import My_def
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import statistics
from cobra.flux_analysis import production_envelope
from cobra.flux_analysis import flux_variability_analysis
import pandas as pd

from scipy.spatial import ConvexHull
import os

os.chdir('../ComplementaryData/')

iML1515 = cobra.io.read_sbml_model('Case3_iML1515/iML1515.xml')

# %%

iML1515_2 = iML1515.copy()

iML1515_2.reactions.get_by_id('EX_o2_e').bounds = (-20, 1000.0)

iML1515_2.optimize()

columns = ['EX_glc__D_e', 'BIOMASS_Ec_iML1515_core_75p37M', 'EX_o2_e', 'EX_ac_e', 'EX_for_e', 'EX_etoh_e',
           'EX_lac__D_e', 'EX_succ_e', ]
result = pd.DataFrame(columns=columns)

for i in range(0, 33, 1):
    iML1515_2.reactions.get_by_id('EX_glc__D_e').lower_bound = -i
    solution = iML1515_2.optimize()
    if solution.fluxes['BIOMASS_Ec_iML1515_core_75p37M'] > 0:
        result = result.append(solution.fluxes[columns])
    # print(solution.objective_value)

# plt.rcParams["font.family"] = "Times New Roman"
fig, ax = plt.subplots(1, 1)
x = result['BIOMASS_Ec_iML1515_core_75p37M']
label = {'EX_glc__D_e': 'Glucose', 'EX_o2_e': 'Oxygen',
         'EX_ac_e': 'Acetate', 'EX_for_e': 'Formate', 'EX_etoh_e': 'Ethanol',
         'EX_lac__D_e': 'Lactate', 'EX_succ_e': 'Succinate'}
i = 0
for met_i in columns:
    if met_i != 'BIOMASS_Ec_iML1515_core_75p37M':
        ax.plot(x, abs(result[met_i] + i * 0.01), label=label[met_i], linewidth=2, alpha=1)
    i += 1
ax.set_xlabel('Biomass rate ', fontsize=15)
ax.set_ylabel('Metabolites rates ', fontsize=15)
# ax.set_yticklabels( fontsize=18)
ax.legend(fontsize=12)
if iML1515_2.reactions.get_by_id('EX_o2_e').lower_bound == 0:
    ax.set_title(u'$E.coli$ Anaerobic', fontsize=18)
else:
    ax.set_title('$E.coli$ Aerobic', fontsize=18)
    ax.set_xlim(0, 1.8)
fig.show()

print('\n---------- Loading experiment data ... ---------- ')
experiment_data_df = pd.read_csv('Case1_ecoli_reduced/ecoli_reduce_experiment_data.txt', delimiter='\t', header=0)
data_cloum_name = ['glc', 'biomass', 'ac', 'for', 'etoh', 'lac', 'succ']

experiment_data_df = experiment_data_df[data_cloum_name]
experiment_data_df.columns = ['EX_glc__D_e', 'BIOMASS_Ec_iML1515_core_75p37M', 'EX_ac_e', 'EX_for_e', 'EX_etoh_e',
                              'EX_lac__D_e', 'EX_succ_e', ]

experiment_yield = {}
for met_i in experiment_data_df.columns:
    x = (experiment_data_df['BIOMASS_Ec_iML1515_core_75p37M'].values[-1] -
         experiment_data_df['BIOMASS_Ec_iML1515_core_75p37M'].values[-0]) / \
        (experiment_data_df['EX_glc__D_e'].values[0] - experiment_data_df['EX_glc__D_e'].values[-1])
    y = (experiment_data_df[met_i].values[-1] - experiment_data_df[met_i].values[-0]) / \
        (experiment_data_df['EX_glc__D_e'].values[0] - experiment_data_df['EX_glc__D_e'].values[-1])
    experiment_yield[met_i] = [x, y]

fig, ax = plt.subplots(1, 1)
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

x = abs(result['BIOMASS_Ec_iML1515_core_75p37M'] / result['EX_glc__D_e'])
label = {'EX_glc__D_e': 'Glucose', 'EX_o2_e': 'Oxygen',
         'EX_ac_e': 'Acetate', 'EX_for_e': 'Formate', 'EX_etoh_e': 'Ethanol',
         'EX_lac__D_e': 'Lactate', 'EX_succ_e': 'Succinate'}

# points_exp = ax.plot(experiment_datas[-1][0], experiment_datas[-1][index], '^', color='red',
#                      alpha=1, markersize=15)
# points_exp = ax.plot(experiment_datas[-2][0], experiment_datas[-2][index], '^', color='red',
#                      alpha=1, markersize=15)
# points_exp = ax.plot(experiment_datas[-3][0], experiment_datas[-3][index], '^', color='red',
#                 alpha=1,  markersize=15)
for index in range(len(columns)):
    met_i = columns[index]
    if met_i == 'EX_o2_e':
        continue
    if met_i != 'BIOMASS_Ec_iML1515_core_75p37M':
        if index >= 1:
            index = index - 1
        ax.plot(x, abs(result[met_i] / result['EX_glc__D_e']), 'o', label=label[met_i], alpha=0.8, color=colors[index])

        # points_exp = ax.plot(abs(experiment_yield[met_i][0]), abs(experiment_yield[met_i][1]), '^',color = colors[index],)

ax.set_xlabel('Biomass/ Glucose yield ', fontsize=15)
ax.set_ylabel('Metabolites/ Glucose yield ', fontsize=15)
# ax.set_yticklabels( fontsize=18)
ax.legend(fontsize=12)
if iML1515_2.reactions.get_by_id('EX_o2_e').lower_bound == 0:
    ax.set_title('$E.coli$ Anaerobic', fontsize=18)
else:
    ax.set_title('$E.coli$ Aerobic', fontsize=18)
    ax.set_ylim(-0.050, 1.8)
fig.show()

# %%

iML1515_3 = iML1515.copy()

iML1515_3.reactions.get_by_id('EX_o2_e').bounds = (0, 1000.0)

iML1515_3.optimize()

columns = ['EX_glc__D_e', 'EX_ac_e', 'EX_for_e', 'EX_etoh_e',
           'EX_lac__D_e', 'EX_succ_e', ]
column_max = [i + 'max' for i in columns]
column_min = [i + 'min' for i in columns]

result_max_min = {}
for i in ['biomass'] + column_max + column_min:
    result_max_min[i] = []

iML1515_3.reactions.get_by_id('EX_glc__D_e').lower_bound = -30

iML1515_3.objective = 'BIOMASS_Ec_iML1515_core_75p37M'
iML1515.objective_direction = 'max'

solution = iML1515_3.optimize()

for i in np.arange(0, solution.objective_value, 0.01):
    iML1515_3.reactions.get_by_id('BIOMASS_Ec_iML1515_core_75p37M').bounds = (i, i)
    result_max_min['biomass'].append(i)
    for i in columns:
        iML1515_3.objective = i
        iML1515_3.objective_direction = 'max'
        max_i = iML1515_3.optimize().objective_value
        result_max_min[i + 'max'].append(max_i)

        iML1515_3.objective_direction = 'min'
        min_i = iML1515_3.optimize().objective_value
        result_max_min[i + 'min'].append(min_i)

    # print(solution.objective_value)

for i in range(0, len(columns)):
    if i == 0:
        continue
    met_i = columns[i]
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(result_max_min['biomass'], result_max_min[met_i + 'max'], 'x--', color='tab:blue', alpha=0.8)
    ax.plot(result_max_min['biomass'], result_max_min[met_i + 'min'], 'x--', color='tab:blue', alpha=0.8)

    ax.set_ylabel(label[met_i] + ' rate ')
    ax.set_xlabel('Biomass rate ')
    fig.show()

# %%
iML1515.reactions.get_by_id('EX_o2_e').bounds = (0.0, 1000.0)
# % Constrain the phosphotransferase system
# model = changeRxnBounds(model, 'GLCabcpp', -1000, 'l');
# model = changeRxnBounds(model, 'GLCptspp', -1000, 'l');
# model = changeRxnBounds(model, 'GLCabcpp', 1000, 'u');
# model = changeRxnBounds(model, 'GLCptspp', 1000, 'u');
# model = changeRxnBounds(model, 'GLCt2pp', 0, 'b');
iML1515.reactions.get_by_id('GLCabcpp').bounds = (-1000.0, 1000.0)
iML1515.reactions.get_by_id('GLCptspp').bounds = (-1000.0, 1000.0)
iML1515.reactions.get_by_id('GLCt2pp').bounds = (0.0, 0.0)
# iML1515.reactions.get_by_id('EX_o2_e').bounds = (0.0,1000.0)
solution = iML1515.optimize()
print(solution)

model = iML1515.copy()
production_rea_ids_x = ['BIOMASS_Ec_iML1515_core_75p37M', ]
production_rea_ids_y = ['EX_ac_e', 'EX_for_e', 'EX_etoh_e', 'EX_lac__D_e', 'EX_succ_e', ]
carbon_source_rea_id = 'EX_glc__D_e'

steps = 10
carbon_uptake_direction = -1
model.reactions.get_by_id(carbon_source_rea_id).bounds = (-10, -0.001)

# all modes
# yield_normalized_df, fluxes_all, hull_index_all = GEM2pathways.get_yield_space_multi(model, production_rea_ids_x,
#                                                                                      production_rea_ids_y,
#                                                                                      carbon_source_rea_id,
#                                                                                      steps=steps,
#                                                                                      carbon_uptake_direction=carbon_uptake_direction,
#                                                                                      draw=True)
# yield_normalized_df_hull = yield_normalized_df[yield_normalized_df.columns[hull_index_all]]
#
#
# yield_normalized_df_hull.to_csv('../ComplementaryData/Case3_iML1515/iML151_yield_normalized_df_hull.csv', sep=',',index=True)


yield_normalized_df_hull_ = pd.read_csv('../../ComplementaryData/Case3_iML1515/iML151_yield_normalized_df_hull.csv',
                                        index_col=0)

our_all_points = yield_normalized_df_hull_.values.T
our_all_points = our_all_points[abs(our_all_points[:, 0]) > 1e-10, :]
our_all_points = our_all_points[:, 1:]  # exclude the carbon colume

# NOTE: according to the our_estimated_datas, set the biomacc_mas_coefficient = 5.0 it means 1 in model , 1 *5.0 in experiment data
coefficient = 3

our_all_points[:, 0] = our_all_points[:, 0] * coefficient

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

experiment_datas = []  # TODO experiment data!!!
for i in range(0, experiment_data_df_trimed_values.shape[0]):
    experiment_data = list(experiment_data_df_trimed_values[i, :])
    # experiment_data[0] = ''
    experiment_datas.append(experiment_data)

qhull_options = 'QJ Qx A0.9999999'
cutoff_persent = 1
our_indexes, our_weights, our_estimated_datas, our_in_hulls = ConvexHull_yield.pipeline_mya(our_all_points,
                                                                                            experiment_datas,
                                                                                            qhull_options=qhull_options,
                                                                                            method=1,
                                                                                            cutoff_persent=cutoff_persent)

our_hull = ConvexHull(our_all_points[:, 1:], qhull_options=qhull_options)
ConvexHull_yield.point_in_hull(experiment_datas[-1][1:], our_hull, tolerance=1e-12)

hull_cutoff_index_99 = [0, 2, 3, 9, 10, 11, 16, 17, 18, 19, 20, 23, 26, 30, 33, 34, 35]
hull_active_index = our_indexes[-1][
    -2]  # list(set(our_indexes[-1][-1]) | set(our_indexes[-1][-2]) | set(our_indexes[-1][-3]))

our_all_points_hull_1 = our_all_points[our_indexes[0], :]
our_all_points_hull_99 = our_all_points[hull_cutoff_index_99, :]
our_all_points_hull_act = our_all_points[hull_active_index, :]

# %% <plot initial yield>: TODO not divided by carbon, so not yield

yield_rea_ids_name = ['Acetate', 'Formate', 'Ethanol', 'Lacate', 'Succinate']
markersize_list = [10, 10, 10]
lw_sizelist = [2, 2, 2]


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

    # xy_EFMs = EFMs_all_points[:, [0,index]]
    # xy_EFVs = EFVs_all_points[:, [0,index]]
    # xy_FBAem = FBAem_all_points[:, [0,index]]
    xy_our = our_all_points[:, [0, index]]

    colors_list = ['tab:blue', 'tab:orange', 'tab:red', 'tab:orange']

    # points_EFMs = ax.plot(xy_EFMs[:,0], xy_EFMs[:,1],       '.', markerfacecolor='none', color=colors_list[0],
    #                 alpha = 0.5, markersize=3)
    # points_EFVs = ax.plot(xy_EFVs[:,0], xy_EFVs[:,1],       '.', markerfacecolor='none', color=colors_list[1],
    #                 alpha = 0.5, markersize=3)

    points_our = ax.plot(xy_our[:, 0], xy_our[:, 1], '.', color=colors_list[0],
                         alpha=0.8, markersize=10)

    points_exp = ax.plot(experiment_datas[-1][0], experiment_datas[-1][index], '^', color='red',
                         alpha=1, markersize=15)
    points_exp = ax.plot(experiment_datas[-2][0], experiment_datas[-2][index], '^', color='red',
                         alpha=1, markersize=15)
    points_exp = ax.plot(experiment_datas[-3][0], experiment_datas[-3][index], '^', color='red',
                         alpha=1, markersize=15)
    # points_FBAem = ax.plot(xy_FBAem[:,0], xy_FBAem[:,1],    's', color=colors_list[3],
    #                 alpha=1,  markersize=8)

    draw_hull = True

    if draw_hull:

        xy_our_hull = our_all_points_hull_1[:, [0, index]]
        ax.plot(xy_our_hull[:, 0], xy_our_hull[:, 1], 'x', markerfacecolor='none', color=colors_list[0],
                alpha=0.8, markersize=markersize_list[0])

        hull_our = ConvexHull(xy_our_hull, qhull_options=qhull_options)

        for simplex in hull_our.simplices:
            line_our_1 = ax.plot(xy_our_hull[simplex, 0], xy_our_hull[simplex, 1], 'x-', markerfacecolor='none',
                                 color=colors_list[0],
                                 alpha=0.5, markersize=10, lw=lw_sizelist[0])

        xy_our_hull = our_all_points_hull_99[:, [0, index]]
        ax.plot(xy_our_hull[:, 0], xy_our_hull[:, 1], 'x', markerfacecolor='none', color=colors_list[1],
                alpha=0.8, markersize=markersize_list[1])

        hull_our = ConvexHull(xy_our_hull, qhull_options=qhull_options)

        for simplex in hull_our.simplices:
            line_our_99 = ax.plot(xy_our_hull[simplex, 0], xy_our_hull[simplex, 1], 'x-', markerfacecolor='none',
                                  color=colors_list[1],
                                  alpha=0.5, markersize=10, lw=lw_sizelist[1])

        xy_our_hull = our_all_points_hull_act[:, [0, index]]
        ax.plot(xy_our_hull[:, 0], xy_our_hull[:, 1], 'x', markerfacecolor='none', color=colors_list[2],
                alpha=0.8, markersize=markersize_list[2])

        hull_our = ConvexHull(xy_our_hull, qhull_options=qhull_options)
        for simplex in hull_our.simplices:
            line_our_ac = ax.plot(xy_our_hull[simplex, 0], xy_our_hull[simplex, 1], 'x-', markerfacecolor='none',
                                  color=colors_list[2],
                                  alpha=0.5, markersize=10, lw=lw_sizelist[2])


    else:
        # line_EFMs = points_EFMs
        # line_EFVs = points_EFVs
        # line_FBAem = points_FBAem
        line_our = points_our

    ax.set_ylabel(yield_rea_ids_name[index - 1] + ' / Glc', fontsize=18)

ax.set_xlabel('Yield Biomass/Glucose', fontsize=18)
# fig.legend((line_EFMs[0],line_EFVs[0],line_our[0],line_FBAem[0]),('EFMs','EFVs','This study','FBA modes'),bbox_to_anchor=(0.6, 0.3), loc='upper left', borderaxespad=0.,prop={'size': 18})
fig.legend((points_our[0], points_exp[0], line_our_1[0], line_our_99[0], line_our_ac[0]),
           ('All pathways', 'Experiment data', 'Convex hull', 'Convex hull %99', 'Convex hull with data'),
           bbox_to_anchor=(0.55, 0.35), loc='upper left', borderaxespad=0., prop={'size': 20})
# fig.savefig('../ComplementaryData/Case3_iML1515/4.png',format ='png')
fig.show()
