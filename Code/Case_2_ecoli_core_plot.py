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
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull

import seaborn as sns
import matplotlib

matplotlib.rc('font', family="Arial")
matplotlib.rcParams["font.family"] = 'Arial'  # 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
os.chdir('../Data/')

# %% <method EFMs > can not calculate
print('\n---------- Loading EFMs ... ---------- ')
core_EFMs_z = np.genfromtxt('Case2_1_ecoli_core/e_coli_core_EFMs_standardized.csv', delimiter=',')
EFMs_all_points = core_EFMs_z.T
EFMs_all_points = EFMs_all_points[:, 1:]  # modes x reas
experiment_datas = []
print('EFMs:', EFMs_all_points.shape)

# %% <method EFVs > can not calculate
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

# %% < our points>

yield_normalized_df_ = pd.read_csv('Case2_1_ecoli_core/ecoli_core_our_yield_normalized_df.csv', sep=',',
                                   index_col=0, )

yield_normalized_df_hull_ = pd.read_csv('Case2_1_ecoli_core/ecoli_core_our_yield_normalized_df_hull.csv', sep=',',
                                        index_col=0, )

our_all_points = yield_normalized_df_hull_.values.T
our_all_points = our_all_points[abs(our_all_points[:, 0]) > 1e-10, :]
our_all_points = our_all_points[:, 1:]  # exclude the carbon colume
print('Our method:', our_all_points.shape)

# %% fig 4a
qhull_options = 'QJ Qx A0.9999999'
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
# fig.savefig('Case2_1_ecoli_core/fig4a_ac_biomass.jpg', bbox_inches='tight', dpi=600)

fig.show()
# %%
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
# fig.savefig('Case2_1_ecoli_core/fig4a_all_mets_biomass.jpg', bbox_inches='tight', dpi=600)

# fig.savefig('Case2_1_ecoli_core/fig4a_all_mets_biomass.svg', bbox_inches='tight', transparent=True)

fig.show()

# %% <Figure 4>

yield_rea_ids_name = ['Acetate', 'Ethanol', 'Lactate', 'Succinate', 'Formate', ]  # ,
yield_rea_ids_lables = ['Acetate', 'Ethanol', 'Lactate', 'Succinate', 'Formate', ]

cmap_list = ['Pastel1', 'Pastel2', 'Paired', 'Accent',
             'Dark2', 'Set1', 'Set2', 'Set3',
             'tab10', 'tab20', 'tab20b', 'tab20c']

plt.get_cmap('Pastel1')


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
        qhull_options = 'QJ'
        hull_EFMs = ConvexHull(xy_EFMs, qhull_options=qhull_options)
        hull_EFVs = ConvexHull(xy_EFVs, qhull_options=qhull_options)
        hull_our = ConvexHull(xy_our, qhull_options=qhull_options)

        # print('index_list:', index, lable)
        # print('hull_EFMs area: ', hull_EFMs.area)
        # print('hull_EFVs area: ', hull_EFVs.area)
        # print('hull_all  area: ', hull_our.area)
        #
        # print('hull_EFMs volume: ', hull_EFMs.volume)
        # print('hull_EFVs volume: ', hull_EFVs.volume)
        # print('hull_all  volume: ', hull_our.volume)

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

    ax.set_ylabel(lable + ' Yield', fontsize=8, family='Arial', labelpad=1)
    ax.set_xlabel('Biomass Yield', fontsize=8, family='Arial', )

fig.legend((line_EFMs[0], line_EFVs[0], line_our[0],),
           ('EFMs', 'EFVs', 'Opt-yield-FBA'),
           loc=4, bbox_to_anchor=(0.9, 0.11),
           # loc=2, bbox_to_anchor=(0.15, 0.93),
           prop={'family': 'Arial', 'size': 8})
fig.tight_layout(pad=1, h_pad=0.5, w_pad=0.5)

# fig.savefig('Case2_1_ecoli_core/fig4a_ac_biomass.pdf', bbox_inches='tight')
# fig.savefig('Case2_1_ecoli_core/fig4a+b_all_biomass.svg', bbox_inches='tight')
fig.savefig('Case2_1_ecoli_core/fig4a_ac_biomass.jpg', bbox_inches='tight', dpi=900)
fig.savefig('Case2_1_ecoli_core/fig4a_ac_biomass.png', bbox_inches='tight', dpi=600)

fig.show()

# %% <Figure simuation>
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

sol = np.loadtxt('Case2_1_ecoli_core/ecoli_core_sol.csv', delimiter=',')

metabObj = ['glc', 'biomass', 'ac', 'for', 'etoh', 'lac', 'succ']
labels = ['Glucose', 'Biomass', 'Acetate', 'Formate', 'Ethanol', 'Lactate', 'Succinate']

tStart = 0.0  # DefineTime
tStop = 10
tStep = 0.1
tspan = np.linspace(tStart, tStop, int(((tStop - tStart) / tStep) + 1))

fig = plt.figure(figsize=(6, 2.5))
ax = fig.add_subplot(111)
color_list = plt.cm.tab10(np.linspace(0, 1, 11))

for index in range(0, 3):
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
plt.yticks(fontname="Arial", fontsize=8)
plt.xticks(fontname="Arial", fontsize=8)
L = fig.legend(loc='lower left', bbox_to_anchor=(1, 0.21), ncol=1, prop={'family': 'Arial', 'size': 8} )
plt.setp(L.texts, family='Arial')
fig.tight_layout()
fig.savefig('Case2_1_ecoli_core/fig4c_core_simulate.pdf', bbox_inches='tight')
fig.savefig('Case2_1_ecoli_core/fig4c_core_simulate.jpg', bbox_inches='tight', dpi=600)

fig.show()
# %%
print('title: \t 2-dimensional yield space',
      'hull_EFVs.volume', 'hull_EFMs.volume', 'hull_our.volume',
      'hull_EFVs.area', 'hull_EFMs.area', 'hull_our.area'
      )
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
    xy_EFMs = EFMs_all_points[:, [0, index]]
    xy_EFVs = EFVs_all_points[:, [0, index]]

    our_all_points = yield_normalized_df_.values.T
    our_all_points = our_all_points[:, 1:]
    xy_our = our_all_points[:, [0, index]]

    qhull_options = 'QJ'
    hull_EFMs = ConvexHull(xy_EFMs, qhull_options=qhull_options)
    hull_EFVs = ConvexHull(xy_EFVs, qhull_options=qhull_options)
    hull_our = ConvexHull(xy_our, qhull_options=qhull_options)

    hull_EFVs_volume = hull_EFVs.volume
    hull_EFVs_area = hull_EFVs.area
    volume_list = [str(round(hull_EFVs_volume, 3)), str(round(hull_our.volume, 3)),
                   str(round(hull_our.volume / hull_EFVs_volume, 3)),
                   str(round(hull_EFMs.volume, 3)), str(round(hull_EFMs.volume / hull_EFVs_volume, 3)),
                   ]
    area_list = [str(round(hull_EFVs_area, 3)), str(round(hull_our.area, 3)),
                 str(round(hull_our.area / hull_EFVs_area, 3)),
                 str(round(hull_EFMs.area, 3)), str(round(hull_EFMs.area / hull_EFVs_area, 3)),
                 ]

    print(lable + '\t' + '\t'.join(volume_list) + '\t' + '\t'.join(area_list))

# %%
from itertools import combinations
import pandas as pd

columns = ['Products combinations', 'Combinations index',
           'EFVs volume value',
           'opt-yield-FBA volume value', 'opt-yield-FBA volume percentage', 'EFMs volume value',
           'EFMs volume percentage',
           'EFVs area value',
           'opt-yield-FBA area value', 'opt-yield-FBA area percentage', 'EFMs area value',
           'EFMs area percentage']
mult_yield_space_df = pd.DataFrame(columns=columns, )
products_list = ['Biomass', 'Acetate', 'Ethanol', 'Lactate', 'Succinate', 'Formate', ]

for i in [2, 3, 4, 6]:
    comb = combinations([0, 1, 2, 3, 4, 5], i)

    for index_list in list(comb):
        # index_list = [2, 1, 5]
        print('index_list:', index_list)
        lable_list = ['-'.join([products_list[i] for i in index_list])] + [index_list]

        xy_EFMs = EFMs_all_points[:, index_list]
        xy_EFVs = EFVs_all_points[:, index_list]
        xy_our = our_all_points[:, index_list]

        hull_EFMs = ConvexHull(xy_EFMs, qhull_options=qhull_options)
        hull_EFVs = ConvexHull(xy_EFVs, qhull_options=qhull_options)
        hull_our = ConvexHull(xy_our, qhull_options=qhull_options)

        hull_EFVs_volume = hull_EFVs.volume
        hull_EFVs_area = hull_EFVs.area
        volume_list = [hull_EFVs_volume,
                       hull_our.volume, hull_our.volume / hull_EFVs_volume,
                       hull_EFMs.volume, hull_EFMs.volume / hull_EFVs_volume,
                       ]
        area_list = [hull_EFVs_area,
                     hull_our.area, hull_our.area / hull_EFVs_area,
                     hull_EFMs.area, hull_EFMs.area / hull_EFVs_area,
                     ]

        mult_yield_space_df.loc[len(mult_yield_space_df)] = lable_list + volume_list + area_list

mult_yield_space_df.to_csv('Case2_1_ecoli_core/multidimensional_yield_space_table.tsv', sep='\t')

values_col = ['opt-yield-FBA volume percentage',
              'EFMs volume percentage',
              'opt-yield-FBA area percentage',
              'EFMs area percentage']

line_index = 51
mult_yield_space_df.loc[line_index, 'Products combinations'] = 'mean 2d biomass'
mult_yield_space_df.loc[line_index, values_col] = mult_yield_space_df.loc[0:4, values_col].mean()
line_index += 1
mult_yield_space_df.loc[line_index, 'Products combinations'] = 'std 2d biomass'
mult_yield_space_df.loc[line_index, values_col] = mult_yield_space_df.loc[0:4, values_col].std()
line_index += 1
mult_yield_space_df.loc[line_index, 'Products combinations'] = 'mean 2d'
mult_yield_space_df.loc[line_index, values_col] = mult_yield_space_df.loc[0:14, values_col].mean()
line_index += 1
mult_yield_space_df.loc[line_index, 'Products combinations'] = 'std 2d'
mult_yield_space_df.loc[line_index, values_col] = mult_yield_space_df.loc[0:14, values_col].std()
line_index += 1
mult_yield_space_df.loc[line_index, 'Products combinations'] = 'mean 3d'
mult_yield_space_df.loc[line_index, values_col] = mult_yield_space_df.loc[15:34, values_col].mean()
line_index += 1
mult_yield_space_df.loc[line_index, 'Products combinations'] = 'std 3d'
mult_yield_space_df.loc[line_index, values_col] = mult_yield_space_df.loc[15:34, values_col].std()
line_index += 1
mult_yield_space_df.loc[line_index, 'Products combinations'] = 'mean 4d'
mult_yield_space_df.loc[line_index, values_col] = mult_yield_space_df.loc[35:49, values_col].mean()
line_index += 1
mult_yield_space_df.loc[line_index, 'Products combinations'] = 'std 4d'
mult_yield_space_df.loc[line_index, values_col] = mult_yield_space_df.loc[35:49, values_col].std()
line_index += 1
mult_yield_space_df.loc[line_index, 'Products combinations'] = 'mean all'
mult_yield_space_df.loc[line_index, values_col] = mult_yield_space_df.loc[0:50, values_col].mean()
line_index += 1
mult_yield_space_df.loc[line_index, 'Products combinations'] = 'std all'
mult_yield_space_df.loc[line_index, values_col] = mult_yield_space_df.loc[0:50, values_col].std()

mult_yield_space_df.to_csv('Case2_1_ecoli_core/multidimensional_yield_space_table.tsv', sep='\t')
mult_yield_space_df.to_csv('Case2_1_ecoli_core/multidimensional_yield_space_table.tsv', sep='\t')
mult_yield_space_df.to_csv('../../Paper_cybernetic_modeling/Table_S5_multidimensional_yield_space_table.tsv', sep='\t')
