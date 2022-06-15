#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 3/14/22

"""Case_1_reduced_ecoli_plot.py
:description : script
:param : 
:returns: 
:rtype: 
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull
import ConvexHull_yield

import seaborn as sns
import matplotlib

matplotlib.rc('font', family="Arial")
matplotlib.rcParams["font.family"] = 'Arial'  # 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
os.chdir('../Data/')

# %% < our points>

yield_normalized_df_ = pd.read_csv('Case3_iML1515/iML151_yield_normalized_df.csv', sep=',',
                                   index_col=0, )

yield_normalized_df_hull_ = pd.read_csv('Case3_iML1515/iML151_yield_normalized_df_hull.csv',
                                        index_col=0)

our_all_points = yield_normalized_df_hull_.values.T
our_all_points = our_all_points[abs(our_all_points[:, 0]) > 1e-10, :]
our_all_points = our_all_points[:, 1:]  # exclude the carbon colume

coefficient = 3
our_all_points[:, 0] = our_all_points[:, 0] * coefficient
# %% <experimental data>

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

hull_active_index = [3, 17, 20, 30, 32, 33, 34]
our_all_points_hull_1 = our_all_points[our_indexes[0], :]
our_all_points_hull_99 = our_all_points[hull_cutoff_index_99, :]
our_all_points_hull_act = our_all_points[hull_active_index, :]

# %% <Figure 5a>

yield_rea_ids_name = ['Acetate', 'Ethanol', 'Lactate', 'Succinate', 'Formate', ]  # ,
yield_rea_ids_lables = ['Acetate', 'Ethanol', 'Lactate', 'Succinate', 'Formate', ]
markersize_list = [6, 5, 5]
lw_sizelist = [2, 1, 1]
yield_normalized_df_ = pd.read_csv('Case3_iML1515/iML151_yield_normalized_df.csv', sep=',',
                                   index_col=0, )
our_all_points = yield_normalized_df_.values.T
our_all_points = our_all_points[:, 1:]
our_all_points[:, 0] = our_all_points[:, 0] * coefficient

cmap_list = ['Pastel1', 'Pastel2', 'Paired', 'Accent',
             'Dark2', 'Set1', 'Set2', 'Set3',
             'tab10', 'tab20', 'tab20b', 'tab20c']

plt.get_cmap('Pastel1')
colors_list = ['tab:blue', 'tab:orange', 'tab:red', 'tab:orange']


def trim_axs(axs, N):
    """little helper to massage the axs list to have correct length..."""
    axs = axs.flat
    for ax in axs[N:]:
        ax.remove()
    return axs[:N]


cols = 2
rows = 3
figsize = (7, 4)
fig, axs = plt.subplots(cols, rows, figsize=figsize)
axs = trim_axs(axs, len(yield_rea_ids_name))

experiment_points = np.array(experiment_datas)
yield_rea_ids_lables = ['Acetate', 'Ethanol', 'Lactate', 'Succinate', 'Formate', ]

for ax, exmet_reaction in zip(axs, yield_rea_ids_name):
    index = yield_rea_ids_name.index(exmet_reaction) + 1
    lable = yield_rea_ids_lables[index - 1]
    if lable == 'Ethanol':
        index = 3
    elif lable == 'Acetate':
        index = 1
    elif lable == 'Lactate':
        index = 4
    elif lable == 'Succinate':
        index = 5
    elif lable == 'Formate':
        index = 2

    xy_our = our_all_points[:, [0, index]]
    # points_our = ax.plot(xy_our[:, 0], xy_our[:, 1], '.', color='black', markerfacecolor='none',
    #                      alpha=0.5, markersize=2)
    points_our = ax.plot(xy_our[:, 0], xy_our[:, 1], 'o', markerfacecolor='none', color='black',
                         alpha=0.3, markersize=2)

    points_exp = ax.plot(experiment_datas[-1][0], experiment_datas[-1][index], '*', color='red',
                         alpha=0.8, markersize=5)
    points_exp = ax.plot(experiment_datas[-2][0], experiment_datas[-2][index], '*', color='red',
                         alpha=0.8, markersize=5)
    points_exp = ax.plot(experiment_datas[-3][0], experiment_datas[-3][index], '*', color='red',
                         alpha=0.8, markersize=5)

    xy_our_hull = our_all_points_hull_1[:, [0, index]]

    ax.plot(xy_our_hull[:, 0], xy_our_hull[:, 1], '+', markerfacecolor='none', color=colors_list[0],
            alpha=0.8, markersize=markersize_list[0])

    xy_our_hull = our_all_points[:, [0, index]]
    hull_our = ConvexHull(xy_our_hull, qhull_options=qhull_options)

    for simplex in hull_our.simplices:
        line_our_1 = ax.plot(xy_our_hull[simplex, 0], xy_our_hull[simplex, 1], '+-', markerfacecolor='none',
                             color=colors_list[0],
                             alpha=0.5, markersize=5, lw=lw_sizelist[0])

    xy_our_hull = our_all_points_hull_99[:, [0, index]]
    ax.plot(xy_our_hull[:, 0], xy_our_hull[:, 1], 'x', markerfacecolor='none', color=colors_list[1],
            alpha=0.8, markersize=markersize_list[1])

    hull_our = ConvexHull(xy_our_hull, qhull_options=qhull_options)

    for simplex in hull_our.simplices:
        line_our_99 = ax.plot(xy_our_hull[simplex, 0], xy_our_hull[simplex, 1], 'x-', markerfacecolor='none',
                              color=colors_list[1],
                              alpha=1, markersize=5, lw=lw_sizelist[1])

    xy_our_hull = our_all_points_hull_act[:, [0, index]]
    ax.plot(xy_our_hull[:, 0], xy_our_hull[:, 1], 'x', markerfacecolor='none', color=colors_list[2],
            alpha=0.8, markersize=markersize_list[2])

    hull_our = ConvexHull(xy_our_hull, qhull_options=qhull_options)
    for simplex in hull_our.simplices:
        line_our_ac = ax.plot(xy_our_hull[simplex, 0], xy_our_hull[simplex, 1], 'x-', markerfacecolor='none',
                              color=colors_list[2],
                              alpha=0.5, markersize=5, lw=lw_sizelist[2])

    # # ax.set_xlim((-0.01, 0.1))
    # ax.set_ylim((-0.2, 2.2))
    # ax.set_yticks(np.arange(0, 2.1, 1))
    # ax.set_yticklabels([0, 1, 2, ])
    # if index == 2:
    #     ax.set_ylim((-0.2, 4.4))
    #     ax.set_yticks(np.arange(0, 4.1, 1))
    #     ax.set_yticklabels([0, 1, 2, 3, 4])
    #     # ax.set_yticks(np.arange(0, 1.2, 0.2))
    #     # ax.set_yticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

    plt.yticks(fontname="Arial", fontsize=8)
    plt.xticks(fontname="Arial", fontsize=8)
    ax.tick_params(axis="x", labelsize=8, pad=1)
    ax.tick_params(axis="y", labelsize=8, pad=1),

    ax.set_ylabel(lable + ' Yield', fontsize=8, family='Arial', labelpad=1)
    ax.set_xlabel('Biomass Yield', fontsize=8, family='Arial', )

fig.legend((points_our[0], points_exp[0], line_our_1[0], line_our_99[0], line_our_ac[0]),
           ('All pathways', 'Experiment data', 'Convex hull', 'Convex hull 99%', 'Convex hull with data'),
           loc=4, bbox_to_anchor=(0.92, 0.11), prop={'family': 'Arial', 'size': 8})

fig.tight_layout(pad=1, h_pad=0.5, w_pad=0.5)
fig.savefig('Case3_iML1515/fig5a_yield.pdf', bbox_inches='tight')
fig.savefig('Case3_iML1515/fig5a_yield.jpg', bbox_inches='tight', dpi=600)

fig.show()
plt.close()

# %% <Figure S1>
hull_all_index_dic = {}

qhull_options = 'QJ '

print('ConvexHull cutoff ...')
for cutoff_persent in [1, 0.99, 0.97, 0.95, 0.9, 0.80, 0.70, 0.60, 0.5]:
    if cutoff_persent == 1:
        hull_all = ConvexHull(our_all_points, qhull_options=qhull_options)
        hull_all_index = hull_all.vertices
        print('hull_all_index = ', list(hull_all_index), '\n len:\t', len(list(hull_all_index)))
        hull_all_index_dic[1] = list(hull_all_index)
        hull_all_area = hull_all.area
        hull_all_volume = hull_all.volume

        print('len:\t', len(list(hull_all_index)))
        print('hull area: ', round(hull_all.area, 3), '\t', round(hull_all.area / hull_all_area, 4))
        print('hull volume: ', round(hull_all.volume, 3), '\t', round(hull_all.volume / hull_all_volume, 4))
    else:
        cutoff_v = cutoff_persent * hull_all.volume
        hull_cutoff_index = ConvexHull_yield.get_hull_cutoff(our_all_points, hull_all_index, cutoff_v,
                                                             qhull_options=qhull_options, method=2)
        hull_all_index_dic[cutoff_persent] = hull_cutoff_index
        print('hull_cutoff_index %.2f = ' % cutoff_persent, hull_cutoff_index, '\n len:\t',
              len(list(hull_cutoff_index)))
        temp_points = our_all_points[hull_cutoff_index, :]
        temp_hull = ConvexHull(temp_points, qhull_options=qhull_options)
        print('len:\t', len(list(hull_cutoff_index)))
        print('hull area: ', round(temp_hull.area, 3), '\t', round(temp_hull.area / hull_all_area, 4))
        print('hull volume: ', round(temp_hull.volume, 3), '\t', round(temp_hull.volume / hull_all_volume, 4))

'''# 
hull_all_index =  [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 15, 16, 17, 19, 20, 21, 22, 23, 24, 26, 27, 29, 30, 31, 32, 33, 34, 36, 37, 38, 39, 40, 41, 42, 43, 44, 48, 50, 51, 52, 53, 54, 57, 63, 64, 68, 69, 71, 72, 73, 74, 75, 77, 80, 82, 83, 84, 87, 91, 92, 93, 94, 97, 98, 101, 102, 103, 104, 105, 109, 110, 111, 112, 113, 114, 115, 118, 119, 120, 121, 123, 124, 128, 130, 131, 132, 133, 134, 138, 139, 140, 141, 142, 143, 144, 145, 148, 149, 150, 151, 152, 153, 154, 155, 157, 158, 160, 161, 162, 164, 165, 168, 169, 170, 171, 173, 175, 176, 177, 180, 181, 182, 183, 187, 190, 191, 192, 193, 195, 198, 199, 201, 202, 203, 207, 209] 
 len:	 139

hull_cutoff_index 0.99 =  [0, 130, 3, 4, 7, 8, 9, 201, 75, 12, 77, 142, 16, 150, 87, 152, 155, 93, 158, 190] 
 len:	 20
base point is not enough
hull_cutoff_index 0.97 =  [0, 3, 4, 7, 8, 201, 9, 75, 12, 77, 16, 150, 152, 155, 93, 158] 
 len:	 16
base point is not enough
hull_cutoff_index 0.95 =  [0, 3, 4, 7, 8, 201, 9, 75, 12, 77, 16, 150, 152, 155, 157] 
 len:	 15
base point is not enough
hull_cutoff_index 0.90 =  [0, 3, 4, 7, 8, 201, 9, 12, 77, 16, 150, 54, 152, 155] 
 len:	 14
base point is not enough
hull_cutoff_index 0.80 =  [0, 3, 4, 7, 8, 201, 9, 12, 142, 16, 152, 155] 
 len:	 12
base point is not enough
hull_cutoff_index 0.70 =  [0, 3, 4, 7, 8, 201, 9, 12, 16, 145, 152] 
 len:	 11
base point is not enough
hull_cutoff_index 0.60 =  [0, 3, 4, 7, 8, 201, 9, 12, 13, 16, 152] 
 len:	 11
 hull_cutoff_index 0.50 =  [0, 3, 4, 7, 8, 201, 9, 12, 16, 152] 
 len:	 10
'''
# %%
yield_rea_ids_name = ['Acetate', 'Ethanol', 'Lactate', 'Succinate', 'Formate', ]  # ,
yield_rea_ids_lables = ['Acetate', 'Ethanol', 'Lactate', 'Succinate', 'Formate', ]
markersize_list = [6, 5, 5]
yield_normalized_df_ = pd.read_csv('Case3_iML1515/iML151_yield_normalized_df.csv', sep=',',
                                   index_col=0, )
our_all_points = yield_normalized_df_.values.T
our_all_points = our_all_points[:, 1:]
our_all_points[:, 0] = our_all_points[:, 0] * coefficient

color_list = plt.cm.tab10(np.linspace(0, 1, 11))


def trim_axs(axs, N):
    """little helper to massage the axs list to have correct length..."""
    axs = axs.flat
    for ax in axs[N:]:
        ax.remove()
    return axs[:N]


cols = 3
rows = 2
figsize = (7.5, 7)
fig, axs = plt.subplots(cols, rows, figsize=figsize)
axs = trim_axs(axs, len(yield_rea_ids_name))

experiment_points = np.array(experiment_datas)
yield_rea_ids_lables = ['Acetate', 'Ethanol', 'Lactate', 'Succinate', 'Formate', ]
labels = ['Convex hull', 'Convex hull 99%', 'Convex hull 97%',
          'Convex hull 95%', 'Convex hull 90%', 'Convex hull 80%',
          'Convex hull 70%', 'Convex hull 60%', 'Convex hull 50%']

lw_sizelist = [1, 1, 1, 3, 1, 1, 1, 1, 1]
alphas = [0.5, 0.5, 0.5, 0.3, 0.5, 0.5, 0.5, 0.5, 0.5]
markersize_list = [8, 7, 6, 5, 4, 3, 2, 1, 1]

for ax, exmet_reaction in zip(axs, yield_rea_ids_name):
    index = yield_rea_ids_name.index(exmet_reaction) + 1
    lable = yield_rea_ids_lables[index - 1]
    if lable == 'Ethanol':
        index = 3
    elif lable == 'Acetate':
        index = 1
    elif lable == 'Lactate':
        index = 4
    elif lable == 'Succinate':
        index = 5
    elif lable == 'Formate':
        index = 2
    lines_list = []
    for k, v in hull_all_index_dic.items():

        xy_our = our_all_points[:, [0, index]]
        xy_our = xy_our[v, :]
        index_i = list(hull_all_index_dic.keys()).index(k)
        if k == 1:
            points_our = ax.plot(xy_our[:, 0], xy_our[:, 1], 'o', markerfacecolor='none', color='black',
                                 alpha=0.3, markersize=2, label='All pathways')

        xy_our_hull = ConvexHull(xy_our, qhull_options=qhull_options)
        for simplex in xy_our_hull.simplices:
            line_our = ax.plot(xy_our[simplex, 0], xy_our[simplex, 1], 'o-', markerfacecolor='none',
                               color=color_list[index_i], label=labels[index_i],
                               alpha=alphas[index_i], markersize=markersize_list[index_i], lw=lw_sizelist[index_i])
        lines_list.append(line_our[0])

    plt.yticks(fontname="Arial", fontsize=8)
    plt.xticks(fontname="Arial", fontsize=8)
    ax.tick_params(axis="x", labelsize=8, pad=1)
    ax.tick_params(axis="y", labelsize=8, pad=1),

    ax.set_ylabel(lable + '/Glucose', fontsize=10, family='Arial', labelpad=1)
    ax.set_xlabel('Biomass/Glucose', fontsize=10, family='Arial', )

fig.legend([points_our[0]] + lines_list,
           ['All pathways'] + labels,
           loc=4, bbox_to_anchor=(0.85, 0.1), prop={'family': 'Arial', 'size': 8})

fig.tight_layout(pad=1, h_pad=1, w_pad=1)
fig.savefig('../../Paper_cybernetic_modeling/Figure S3_yield.pdf', bbox_inches='tight')
fig.savefig('../../Paper_cybernetic_modeling/Figure S3_yield.jpg', bbox_inches='tight', dpi=600)
fig.savefig('Case3_iML1515/Figure S3_yield.pdf', bbox_inches='tight')
fig.savefig('Case3_iML1515/Figure S3_yield.jpg', bbox_inches='tight', dpi=600)

fig.show()
plt.close()
# %% <Figure simuation>

sol = np.loadtxt('Case3_iML1515/iML1515_sol.csv', delimiter=',')

metabObj = ['glc', 'biomass', 'ac', 'for', 'etoh', 'lac', 'succ']
labels = ['Glucose', 'Biomass', 'Acetate', 'Formate', 'Ethanol', 'Lactate', 'Succinate']

tStart = 0.0  # DefineTime
tStop = 8.5
tStep = 0.1
tspan = np.linspace(tStart, tStop, int(((tStop - tStart) / tStep) + 1))

fig = plt.figure(figsize=(6, 2.5))
ax = fig.add_subplot(111)
color_list = plt.cm.tab10(np.linspace(0, 1, 11))

for index in range(0, 7):
    if index == 1:
        ax1 = ax.twinx()
        ax1.plot(tspan, sol[:, index], color=color_list[index], linewidth=1, label=labels[index])

        experiment_p = ax1.plot(experiment_data_df['time'], experiment_data_df[metabObj[index]], '.',
                                color=color_list[index], alpha=0.8, markersize=5)

    else:
        ax.plot(tspan, sol[:, index], color=color_list[index], linewidth=1, label=labels[index])

        experiment_p = ax.plot(experiment_data_df['time'], experiment_data_df[metabObj[index]], '.',
                               color=color_list[index], alpha=0.8, markersize=5)

ax.set_xlabel('Time (h)', fontsize=10, family='Arial', )
ax.set_ylabel('Concentration (mM)', fontsize=10, family='Arial', )
ax1.set_ylabel('Biomass (g/L)', fontsize=10, family='Arial', )
ax1.set_ylim(0, 2.1)
ax.set_ylim(0, 28)
plt.yticks(fontname="Arial", fontsize=8)
plt.xticks(fontname="Arial", fontsize=8)
L = fig.legend(loc='lower left', bbox_to_anchor=(1, 0.21), ncol=1, fontsize=8, )
plt.setp(L.texts, family='Arial')
fig.tight_layout()
fig.savefig('Case3_iML1515/fig5b_simulate.pdf', bbox_inches='tight')
fig.savefig('Case3_iML1515/fig5b_simulate.jpg', bbox_inches='tight', dpi=600)

fig.show()
