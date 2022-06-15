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

# %% <EFMs > caculated by matlab emftool. and standardized by step0
print('\n---------- Loading EFMs ... ---------- ')
em_z = np.genfromtxt('Case1_ecoli_reduced/ecoli_reduced_EFMs_standardized.csv', delimiter=',')
EFMs_all_points = em_z.T
EFMs_all_points = EFMs_all_points[:, 1:]  # modes x reas
print('EFMs:', EFMs_all_points.shape)

# %% <FBA modes> from reference and standardized by step0
print('\n---------- Loading FBA modes ... ---------- ')
FBAMs_z = np.genfromtxt('Case1_ecoli_reduced/ecoli_reduced_FBAMs_standardized.csv', delimiter=',')
FBAMs_all_points = FBAMs_z.T  # modes x rea
FBAMs_all_points = FBAMs_all_points[:, 1:]
print('FBA models:', FBAMs_all_points.shape)

# %% < our points>
fluxes_2d = pd.read_csv('Case1_ecoli_reduced/ecoli_reduced_our_fluxes_2d_df.csv', sep=',', index_col=0, )
yield_normalized_df_ = pd.read_csv('Case1_ecoli_reduced/ecoli_reduced_our_yield_normalized_df.csv', sep=',',
                                   index_col=0, )
yield_normalized_df_hull_ = pd.read_csv('Case1_ecoli_reduced/ecoli_reduced_our_yield_normalized_df_hull.csv', sep=',',
                                        index_col=0, )

points_2d = fluxes_2d.values
points_2d = points_2d / points_2d[0, :]
points_2d = points_2d[[8, 11], :]
points_2d = points_2d.T
our_all_points = yield_normalized_df_.values.T
our_all_points = our_all_points[:, 1:]

our_key_points = yield_normalized_df_hull_.values.T
our_key_points = our_key_points[abs(our_key_points[:, 0]) > 1e-10, :]
our_key_points = our_key_points[:, 1:]  # exclude the carbon colume!!!
our_key_points = our_key_points

# %%fig 2a
figsize = (4, 2.5)
fig, ax = plt.subplots(figsize=figsize)
colors = sns.color_palette("Set2")
sns.set_style("ticks", )
xy_our = points_2d

hull_our = ConvexHull(xy_our, qhull_options='QJ Qx A0.9999999')
points_our_all = ax.plot(xy_our[:, 1], xy_our[:, 0], 'o', markerfacecolor='none', color='black',
                         alpha=0.3, label='This_methd', markersize=3)
for simplex in hull_our.simplices:
    points_our_key = ax.plot(xy_our[simplex, 1], xy_our[simplex, 0], '+', color='tab:blue',
                             alpha=0.8, label='This_methd', markersize=8, linewidth=.5, )
    line_our = ax.plot(xy_our[simplex, 1], xy_our[simplex, 0], '--', color='tab:blue',
                       alpha=0.8, label='This_methd', markersize=8, linewidth=.5, )

xy_EFMs = EFMs_all_points[:, [0, 1]]
points_EFMs = ax.plot(xy_EFMs[:, 0], xy_EFMs[:, 1], 'x', markerfacecolor='none', color='tab:orange',
                      alpha=0.8, label='EFMs', markersize=7)

# all points:
# xy_our = our_all_points[:, [0, 1]]
# points_our_all = ax.plot(xy_our[:, 0], xy_our[:, 1], 'o', markerfacecolor='none', color='black',
#                          alpha=0.5, label='This_methd', markersize=3)
ax.set_xlim((0.01, 0.03))
ax.set_ylim((-0.1, 1))
ax.set_ylabel('Yield: Acetate/Glucose', fontsize=10, family='Arial', )
ax.set_xlabel('Yield: Biomass/Glucose', fontsize=10, family='Arial', )

plt.yticks(fontname="Arial", fontsize=8)
plt.xticks(fontname="Arial", fontsize=8)

fig.legend((points_our_all[0], points_our_key[0], points_EFMs[0], line_our[0]),
           ('All points', 'Extreme points', 'EFMs', 'ConvexHull'),
           # bbox_to_anchor=(0.65, -0.12),ncol=2, loc=8,
           loc=4, bbox_to_anchor=(0.93, 0.2),
           # loc=2, bbox_to_anchor=(0.15, 0.93),
           prop={'family': 'Arial', 'size': 7.5})

fig.tight_layout()
fig.savefig('Case1_ecoli_reduced/fig2a_ac_biomass.pdf', bbox_inches='tight')
fig.savefig('Case1_ecoli_reduced/fig2a_ac_biomass.jpg', bbox_inches='tight')

fig.show()

# %% fig 2b
yield_rea_ids_name = ['EX_etoh_e', 'EX_lac__D_e', 'EX_succ_e', 'EX_for_e', ]
yield_rea_ids_lables = ['Ethanol', 'Lactate', 'Succinate', 'Formate', ]


def trim_axs(axs, N):
    """little helper to massage the axs list to have correct length..."""
    axs = axs.flat
    for ax in axs[N:]:
        ax.remove()
    return axs[:N]


cols = 2
rows = 2
figsize = (4, 2.75)
fig, axs = plt.subplots(cols, rows, figsize=figsize)
axs = trim_axs(axs, len(yield_rea_ids_name))

# experiment_points = np.array(experiment_datas)
for ax, exmet_reaction in zip(axs, yield_rea_ids_name):
    index = yield_rea_ids_name.index(exmet_reaction) + 1
    lable = yield_rea_ids_lables[index - 1]
    if lable == 'Ethanol':
        index = 3
    elif lable == 'Lactate':
        index = 4
    elif lable == 'Succinate':
        index = 5
    elif lable == 'Formate':
        index = 2

    xy_EFMs = EFMs_all_points[:, [0, index]]
    xy_our_all = our_all_points[:, [0, index]]
    xy_our_key = our_key_points[:, [0, index]]
    # xy_exp = experiment_points[-4:-1, [0, index]]

    hull_EFMs = ConvexHull(xy_EFMs, qhull_options='QJ')
    hull_our = ConvexHull(xy_our_all, qhull_options='QJ')

    points_our_all = ax.plot(xy_our_all[:, 0], xy_our_all[:, 1], 'o', markerfacecolor='none', color='black',
                             alpha=0.3, label='This_methd', markersize=2)

    points_our_key = ax.plot(xy_our_key[:, 0], xy_our_key[:, 1], '+', color='tab:blue',
                             alpha=0.8, label='This_methd', markersize=5)

    points_EFMs = ax.plot(xy_EFMs[:, 0], xy_EFMs[:, 1], 'x', color='tab:orange',
                          alpha=0.8, label='EFMs', markersize=5)
    # outline
    for simplex in hull_EFMs.simplices:
        line_EFMs = ax.plot(xy_EFMs[simplex, 0], xy_EFMs[simplex, 1], ':', markerfacecolor='none', color='tab:orange',
                            alpha=0.8, label='EFMs', markersize=5, linewidth=.5, )

    for simplex in hull_our.simplices:
        line_our = ax.plot(xy_our_all[simplex, 0], xy_our_all[simplex, 1], '--', color='tab:blue',
                           alpha=0.8, label='This_methd', markersize=3, linewidth=.5, )

    # points_exp = ax.plot(xy_exp[:, 0], xy_exp[:, 1], '^', color='tab:orange',
    #                      alpha=1, label='experiment data', markersize=11)

    ax.set_ylabel(lable + '/Glucose', fontsize=8, family='Arial', labelpad=1)
    ax.set_xlim((0.01, 0.03))
    ax.set_ylim((-0.2, 2))
    if index == 5:
        ax.set_ylim((-0.1, 1))
        ax.set_yticks(np.arange(0, 1.2, 0.2))
        ax.set_yticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

    ax.tick_params(axis="x", labelsize=6, pad=0.5)
    ax.tick_params(axis="y", labelsize=6, pad=0.5),

# fig.legend((points_our_all[0], points_our_key[0], points_EFMs[0], line_our[0]),
#            ('All points', 'Key points', 'EFMs', 'ConvexHull'),
#            bbox_to_anchor=(0.8, -0.22),
#            loc=8, ncol=1, prop={'family': 'Arial', 'size': 8})

fig.tight_layout(pad=1, h_pad=0.5, w_pad=-0.2)
fig.savefig('Case1_ecoli_reduced/fig2b_all_mets_biomass.pdf', bbox_inches='tight', transparent=True)
fig.savefig('Case1_ecoli_reduced/fig2b_all_mets_biomass.svg', bbox_inches='tight', transparent=True)

fig.show()

# %% <Figure 2>
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
    experiment_datas.append(experiment_data_df_trimed_values[i, :])

cmap_list = ['Pastel1', 'Pastel2', 'Paired', 'Accent',
             'Dark2', 'Set1', 'Set2', 'Set3',
             'tab10', 'tab20', 'tab20b', 'tab20c']

plt.get_cmap('Pastel1')

yield_rea_ids_name = ['EX_ac_e', 'EX_for_e', 'EX_etoh_e', 'EX_lac__D_e', 'EX_succ_e', ]


def trim_axs(axs, N):
    """little helper to massage the axs list to have correct length..."""
    axs = axs.flat
    for ax in axs[N:]:
        ax.remove()
    return axs[:N]


cols = 3
rows = 2
figsize = (4, 4.7)
fig, axs = plt.subplots(cols, rows, figsize=figsize)
axs = trim_axs(axs, len(yield_rea_ids_name))

experiment_points = np.array(experiment_datas)
yield_rea_ids_lables = ['Acetate', 'Ethanol', 'Lactate', 'Succinate', 'Formate', ]

for ax, exmet_reaction in zip(axs, yield_rea_ids_name):
    index = yield_rea_ids_name.index(exmet_reaction) + 1
    lable = yield_rea_ids_lables[index - 1]
    if lable == 'Ethanol':
        index = 3
    elif lable == 'Lactate':
        index = 4
    elif lable == 'Succinate':
        index = 5
    elif lable == 'Formate':
        index = 2

    xy_EFMs = EFMs_all_points[:, [0, index]]
    xy_our_all = our_all_points[:, [0, index]]
    xy_our_key = our_key_points[:, [0, index]]
    # xy_exp = experiment_points[-4:-1, [0, index]]

    hull_EFMs = ConvexHull(xy_EFMs, qhull_options='QJ')
    hull_our = ConvexHull(xy_our_all, qhull_options='QJ')

    points_our_all = ax.plot(xy_our_all[:, 0], xy_our_all[:, 1], 'o', markerfacecolor='none', color='black',
                             alpha=0.3, label='This_methd', markersize=2)

    points_our_key = ax.plot(xy_our_key[:, 0], xy_our_key[:, 1], '+', color='tab:blue',
                             alpha=0.8, label='This_methd', markersize=6)
    points_EFMs = ax.plot(xy_EFMs[:, 0], xy_EFMs[:, 1], 'x', color='tab:orange',
                          alpha=0.8, label='EFMs', markersize=5)
    # outline
    for simplex in hull_our.simplices:
        line_our = ax.plot(xy_our_all[simplex, 0], xy_our_all[simplex, 1], '--', color='tab:blue',
                           alpha=0.8, label='This_methd', markersize=3, linewidth=.5, )
    for simplex in hull_EFMs.simplices:
        line_EFMs = ax.plot(xy_EFMs[simplex, 0], xy_EFMs[simplex, 1], ':', markerfacecolor='none', color='tab:orange',
                            alpha=0.8, label='EFMs', markersize=5, linewidth=.5, )

    ax.set_ylabel(yield_rea_ids_name[index - 1] + '/Glucose', fontsize=12)
    ax.set_ylabel(lable + ' Yield', fontsize=8, family='Arial', labelpad=1)
    ax.set_xlabel('Biomass Yield', fontsize=8, family='Arial', )
    ax.set_xlim((0.01, 0.03))
    ax.set_xticks(np.arange(0.01, 0.04, 0.01))
    ax.set_xticklabels([0.01, 0.02, 0.03])
    ax.set_ylim((-0.2, 2))
    ax.set_yticks(np.arange(0, 3, 1))
    ax.set_yticklabels([0, 1, 2])
    if index in {1, 5}:
        ax.set_ylim((-0.1, 1))
        ax.set_yticks(np.arange(0, 1.5, 0.5))
        ax.set_yticklabels([0,'', 1])

    ax.tick_params(axis="x", labelsize=7, pad=0.5)
    ax.tick_params(axis="y", labelsize=7, pad=0.5),

fig.legend((points_our_all[0], points_our_key[0], points_EFMs[0], line_our[0]),
           ('All points', 'Extreme points', 'EFMs', 'ConvexHull'),
           # bbox_to_anchor=(0.65, -0.12),ncol=2, loc=8,
           loc=4, bbox_to_anchor=(0.9, 0.085),
           # loc=2, bbox_to_anchor=(0.15, 0.93),
           prop={'family': 'Arial', 'size': 8})
fig.tight_layout(pad=1, h_pad=0.5, w_pad=-0.2)
fig.savefig('Case1_ecoli_reduced/fig2_all_mets_biomass.pdf', bbox_inches='tight', transparent=True)
fig.savefig('Case1_ecoli_reduced/fig2_all_mets_biomass.svg', bbox_inches='tight', transparent=True)

fig.show()

# %% fig 2a + b

xy_our = points_2d
hull_our = ConvexHull(xy_our, qhull_options='QJ Qx A0.9999999')

colors = sns.color_palette("Set2")
sns.set_style("ticks", )

figsize = (4, 4.7)

fig = plt.figure(figsize=figsize)
gs = gridspec.GridSpec(7, 2)
ax0 = fig.add_subplot(gs[0:3, :])
# ax0.set_aspect('equal', 'box')
points_our_all = ax0.plot(xy_our[:, 1], xy_our[:, 0], 'o', markerfacecolor='none', color='black',
                          alpha=0.3, label='This_methd', markersize=3)
for simplex in hull_our.simplices:
    points_our_key = ax0.plot(xy_our[simplex, 1], xy_our[simplex, 0], '+', color='tab:blue',
                              alpha=0.8, label='This_methd', markersize=8, linewidth=.5, )
    line_our = ax0.plot(xy_our[simplex, 1], xy_our[simplex, 0], '--', color='tab:blue',
                        alpha=0.8, label='This_methd', markersize=8, linewidth=.5, )

xy_EFMs = EFMs_all_points[:, [0, 1]]
points_EFMs = ax0.plot(xy_EFMs[:, 0], xy_EFMs[:, 1], 'x', markerfacecolor='none', color='tab:orange',
                       alpha=0.8, label='EFMs', markersize=7)
ax0.set_xlim((0.01, 0.03))
ax0.set_xlim((0.01, 0.03))
ax0.set_xticks(np.arange(0.01, 0.04, 0.01))
ax0.set_xticklabels([0.01, 0.02, 0.03])
ax0.set_ylim((-0.1, 1))
ax0.set_yticks(np.arange(0, 1.5, 0.5))
# ax0.set_yticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

# ax.tick_params(axis="x", labelsize=6, pad=0.5)
# ax.tick_params(axis="y", labelsize=6, pad=0.5),
ax0.set_ylabel('Yield: Acetate/Glucose', fontsize=8, family='Arial', labelpad=1)
ax0.set_xlabel('Yield: Biomass/Glucose', fontsize=8, family='Arial', labelpad=1)
plt.yticks(fontname="Arial", fontsize=8)
plt.xticks(fontname="Arial", fontsize=8)

yield_rea_ids_name = ['EX_etoh_e', 'EX_lac__D_e', 'EX_succ_e', 'EX_for_e', ]
yield_rea_ids_lables = ['Ethanol', 'Lactate', 'Succinate', 'Formate', ]

for exmet_reaction in yield_rea_ids_name:

    index = yield_rea_ids_name.index(exmet_reaction) + 1
    lable = yield_rea_ids_lables[index - 1]
    if lable == 'Ethanol':
        index = 3
        ax = fig.add_subplot(gs[3:5, 0])
    elif lable == 'Lactate':
        index = 4
        ax = fig.add_subplot(gs[3:5, 1])
    elif lable == 'Succinate':
        index = 5
        ax = fig.add_subplot(gs[5:7, 0])
    elif lable == 'Formate':
        index = 2
        ax = fig.add_subplot(gs[5:7, 1])

    xy_EFMs = EFMs_all_points[:, [0, index]]
    xy_our_all = our_all_points[:, [0, index]]
    xy_our_key = our_key_points[:, [0, index]]
    # xy_exp = experiment_points[-4:-1, [0, index]]

    hull_EFMs = ConvexHull(xy_EFMs, qhull_options='QJ')
    hull_our = ConvexHull(xy_our_all, qhull_options='QJ')

    points_our_all = ax.plot(xy_our_all[:, 0], xy_our_all[:, 1], 'o', markerfacecolor='none', color='black',
                             alpha=0.3, label='This_methd', markersize=2)

    points_our_key = ax.plot(xy_our_key[:, 0], xy_our_key[:, 1], '+', color='tab:blue',
                             alpha=0.8, label='This_methd', markersize=5)

    points_EFMs = ax.plot(xy_EFMs[:, 0], xy_EFMs[:, 1], 'x', color='tab:orange',
                          alpha=0.8, label='EFMs', markersize=5)
    # outline
    for simplex in hull_EFMs.simplices:
        line_EFMs = ax.plot(xy_EFMs[simplex, 0], xy_EFMs[simplex, 1], ':', markerfacecolor='none', color='tab:orange',
                            alpha=0.8, label='EFMs', markersize=5, linewidth=.5, )

    for simplex in hull_our.simplices:
        line_our = ax.plot(xy_our_all[simplex, 0], xy_our_all[simplex, 1], '--', color='tab:blue',
                           alpha=0.8, label='This_methd', markersize=3, linewidth=.5, )

    # points_exp = ax.plot(xy_exp[:, 0], xy_exp[:, 1], '^', color='tab:orange',
    #                      alpha=1, label='experiment data', markersize=11)

    ax.set_ylabel(lable + '/Glucose', fontsize=8, family='Arial', labelpad=1)
    ax.set_xlabel('Biomass/Glucose', fontsize=8, family='Arial', labelpad=1)

    ax.set_xlim((0.01, 0.03))
    ax.set_xticks(np.arange(0.01, 0.04, 0.01))
    ax.set_xticklabels([0.01, 0.02, 0.03])
    ax.set_ylim((-0.2, 2))
    ax.set_yticks(np.arange(0, 3, 1))
    ax.set_yticklabels([0, 1, 2])
    if index in {1, 5}:
        ax.set_ylim((-0.1, 1))
        ax.set_yticks(np.arange(0, 2, 1))
        ax.set_yticklabels([0, 1])

    # ax.tick_params(axis="x", labelsize=6, pad=0.5)
    # ax.tick_params(axis="y", labelsize=6, pad=0.5),

    plt.yticks(fontname="Arial", fontsize=8)
    plt.xticks(fontname="Arial", fontsize=8)
fig.tight_layout(pad=1, h_pad=0.5, w_pad=-0.2)
fig.legend((points_our_all[0], points_our_key[0], points_EFMs[0], line_our[0]),
           ('All points', 'Extreme points', 'EFMs', 'ConvexHull'),
           # bbox_to_anchor=(0.65, -0.12),ncol=2, loc=8,
           loc=4, bbox_to_anchor=(0.94, 0.64),
           # loc=2, bbox_to_anchor=(0.15, 0.93),
           prop={'family': 'Arial', 'size': 7.5})

fig.savefig('Case1_ecoli_reduced/fig2a+b_all_mets_biomass.pdf', bbox_inches='tight', transparent=True)
fig.savefig('Case1_ecoli_reduced/fig2a+b_all_mets_biomass.svg', bbox_inches='tight', transparent=True)
fig.show()

plt.close()

# %% <Figure simuation>
sol = np.loadtxt('Case1_ecoli_reduced/ecoli_reduced_sol.csv', delimiter=',')

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
plt.yticks(fontname="Arial", fontsize=8)
plt.xticks(fontname="Arial", fontsize=8)
ax1.set_ylim(0, 2.1)
ax.set_ylim(0, 28)
L = fig.legend(loc='lower left', bbox_to_anchor=(1, 0.2), ncol=1, fontsize=8, )
plt.setp(L.texts, family='Arial')
fig.tight_layout()
fig.savefig('Case1_ecoli_reduced/fig2c_simulate.pdf', bbox_inches='tight')
fig.savefig('Case1_ecoli_reduced/fig2c_simulate.jpg', bbox_inches='tight')

fig.show()

# %%
xy_EFMs = EFMs_all_points[:, [0, 1, 2, 3, 4, 5]]
xy_our_all = our_all_points[:, [0, 1, 2, 3, 4, 5]]
xy_our_key = our_key_points[:, [0, 1, 2, 3, 4, 5]]

hull_EFMs = ConvexHull(xy_EFMs, qhull_options='QJ')
hull_all = ConvexHull(xy_our_all, qhull_options='QJ')
hull_key = ConvexHull(xy_our_key, qhull_options='QJ')

print('hull_EFMs: ', hull_EFMs.area)
print('hull_all: ', hull_all.area)
print('hull_key: ', hull_key.area)

print('hull_EFMs volume: ', hull_EFMs.volume)
print('hull_all volume: ', hull_all.volume)
print('hull_key volume: ', hull_key.volume)
