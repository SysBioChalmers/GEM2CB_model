#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 3/17/22

"""yeast_three_03_plot.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os
import pickle

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import copy
import seaborn as sns
import matplotlib

matplotlib.rc('font', family="Arial")
matplotlib.rcParams["font.family"] = 'Arial'  # 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']

# %% < def input:>
os.chdir('../../Data/Case_yeast_three_species')

print('\n---------- Loading Experiment Data ... ---------- ')

experiment_data_df_1 = pd.read_excel(r'experiemntdata_to_fit.xlsx', sheet_name='model1')
experiment_data_df_2 = pd.read_excel(r'experiemntdata_to_fit.xlsx', sheet_name='model2')
experiment_data_df_3 = pd.read_excel(r'experiemntdata_to_fit.xlsx', sheet_name='model3')
experiment_data_df_12 = pd.read_excel(r'experiemntdata_to_fit.xlsx', sheet_name='model1_and_model2')
experiment_data_df_23 = pd.read_excel(r'experiemntdata_to_fit.xlsx', sheet_name='model2_and_model3')
# experiment_data_df_13 = pd.read_excel(r'experiemntdata_to_fit.xlsx', sheet_name='model1_and_model3')
# experiment_data_df_123 = pd.read_excel(r'experiemntdata_to_fit.xlsx', sheet_name='model1_and_model2_and_model3')
experiment_data_dfs_copy = [experiment_data_df_1, experiment_data_df_2, experiment_data_df_3,
                            experiment_data_df_12, experiment_data_df_23]
experiment_data_dfs = copy.deepcopy(experiment_data_dfs_copy)

molar_mass = np.array([1040, 46, 180, 150, 180, 180])  # unit conversion : 0.01g/l = 1000*0.01/900 mM(mmol/l)
# for experiment_data_df_i in experiment_data_dfs:
#     experiment_data_df_i.iloc[:, 1:] = experiment_data_df_i.values[:, 1:] * 1000 / molar_mass
#     experiment_data_df_i.iloc[:, 1:] = experiment_data_df_i.values[:, 1:] * molar_mass / 1000

# [sol1, sol2, sol3, sol12, sol23, sol13, sol123] = pickle.load(open('sol_all', 'rb'))
sols = pickle.load(open('sol_all', 'rb'))
results = []
for sol_i in sols:
    result_i = sol_i[:, 0:6] * molar_mass / 1000
    # result_i = sol_i[:, 0:6]
    results.append(result_i)

# time def
tStart = 0.0
tStop = 100
tStep = 1
tspan = np.linspace(tStart, tStop, int(((tStop - tStart) / tStep) + 1))

titles = ['$S. cerevisiae$', '$P. stipitis$', '$K. marxianus$',
          '$S. cerevisiae$ & $P. stipitis$', '$P. stipitis$ & $K. marxianus$', '$S. cerevisiae$ & $K. marxianus$',
          '$S. cerevisiae$ & $P. stipitis$ & $K. marxianus$']
mets_name = list(experiment_data_df_1.columns[1:])  # ['BIOM', 'ETHx', 'GLC', 'XYL', 'MAN', 'GAL']
label_list = ['Biomass', 'Ethanol', 'Glucose', 'Xylose', 'Mannose', 'Galactose']
# colors = sns.color_palette("Set2")
colors = sns.color_palette('muted')
del colors[0]
# colors = sns.color_palette('dark')
sns.set_style("ticks", )
# colors = plt.cm.tab10(np.linspace(0, 1, 12))
for index in range(0, len(results)):
    result = results[index]
    if index in {0, 3}:
        figsize = (8, 2.5)
        fig, axs = plt.subplots(2, 3, figsize=figsize)
    elif index in {6}:
        figsize = (3.85, 2.5)
        fig, axs = plt.subplots(2, 1, figsize=figsize)

    if index < 3:
        ax_0 = axs[0, index]
        ax_1 = axs[1, index]
    elif 2 < index < 6:
        ax_0 = axs[0, index - 3]
        ax_1 = axs[1, index - 3]
    else:
        ax_0 = axs[0]
        ax_1 = axs[1]

    result = results[index]
    for met_i in mets_name:
        ii = mets_name.index(met_i)
        label_i = label_list[ii]
        if met_i == 'BIOM':
            ax_0.plot(tspan, result[:, ii], color=colors[ii], label=label_i, linewidth=1.5)
            ax_0.set_yscale('log')
            # ax_0.plt.yscale()
        elif met_i != 'BIOM':
            ax_1.plot(tspan, result[:, ii], color=colors[ii], label=label_i, linewidth=1.5)

    if index < 5:
        experiment_data_df_i = experiment_data_dfs[index]
        color_i = 0
        for met_i in mets_name:
            ii = mets_name.index(met_i)
            if met_i == 'BIOM':
                ax_0.plot(experiment_data_df_i['t'], experiment_data_df_i[met_i], '.', color=colors[ii], alpha=1,
                          linewidth=1, markersize=5)

                ax_0.set_yscale('log')
                # ax_0.plt.yscale()
            elif met_i != 'BIOM':
                ax_1.plot(tspan, result[:, ii], color=colors[ii], linewidth=1)
                ax_1.plot(experiment_data_df_i['t'], experiment_data_df_i[met_i], '.', color=colors[ii], alpha=1,
                          linewidth=1, markersize=5)

    ax_0.axis(ymin=0.01, ymax=100)
    ax_0.set_yticks([0.01, 1, 100])
    ax_0.set_yticklabels([-2, 0, 2])

    ax_0.set_xticks([0, 50, 100])
    ax_0.set_xticklabels([])
    ax_1.set_xlabel('Time (h)', fontsize=8, family='Arial', )

    if index in {0, 3, 6}:
        ax_0.set_ylabel('Biomass log (g/L)', fontsize=8, family='Arial', )
        ax_1.set_ylabel('Concentration (g/L)', fontsize=8, family='Arial', )
    if index == 6:
        legend0 = ax_0.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, prop={'family': 'Arial', 'size': 8}, )
        legend1 = ax_1.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, prop={'family': 'Arial', 'size': 8}, )
    ax_0.set_title(titles[index], fontsize=8, family='Arial', )
    plt.yticks(fontname="Arial", fontsize=8)
    plt.xticks(fontname="Arial", fontsize=8)

    ax_0.tick_params(axis="x", labelsize=8, pad=-1)
    ax_0.tick_params(axis="y", labelsize=8, pad=1),
    # ax_1.tick_params(axis="x", labelsize=8, pad=4)
    ax_1.tick_params(axis="y", labelsize=8, pad=1),
    if index in {2, 5, 6}:
        fig.tight_layout(pad=1, h_pad=1, w_pad=0.5)
        fig.savefig('fig6a_%s_2.pdf' % titles[index], bbox_inches='tight')
        fig.savefig('fig6a_%s_2.jpg' % titles[index], bbox_inches='tight', dpi=600.)
        fig.show()

# for index in range(0, len(results)):
#     figsize = (2.7, 2.5)
#     if index == 6:
#         figsize = (3.65, 2.5)
#     fig, axs = plt.subplots(2, 1, figsize=figsize)
#     ax_0 = axs[0]
#     ax_1 = axs[1]
#     colors = plt.cm.tab10(np.linspace(0, 1, 11))
#
#     result = results[index]
#     color_i = 0
#     for met_i in mets_name:
#         label_i = label_list[mets_name.index(met_i)]
#         if met_i == 'BIOM':
#             ax_0.plot(tspan, result[:, mets_name.index(met_i)], color=colors[color_i], label=label_i, linewidth=1)
#             ax_0.set_yscale('log')
#             # ax_0.plt.yscale()
#         elif met_i != 'BIOM':
#             ax_1.plot(tspan, result[:, mets_name.index(met_i)], color=colors[color_i], label=label_i, linewidth=1)
#         color_i += 1
#
#     if index < 5:
#         experiment_data_df_i = experiment_data_dfs[index]
#         color_i = 0
#         for met_i in mets_name:
#             # label_i = label_list[mets_name.index(met_i)]
#             if met_i == 'BIOM':
#                 ax_0.plot(experiment_data_df_i['t'], experiment_data_df_i[met_i], '.', color=colors[color_i], alpha=0.8,
#                           linewidth=1, markersize=5)
#
#                 ax_0.set_yscale('log')
#                 # ax_0.plt.yscale()
#             elif met_i != 'BIOM':
#                 ax_1.plot(tspan, result[:, mets_name.index(met_i)], color=colors[color_i], linewidth=1)
#                 ax_1.plot(experiment_data_df_i['t'], experiment_data_df_i[met_i], '.', color=colors[color_i], alpha=0.8,
#                           linewidth=1, markersize=5)
#             color_i += 1
#
#     ax_0.axis(ymin=0.01, ymax=100)
#     ax_0.set_yticks([0.01, 1, 100])
#     ax_0.set_yticklabels([-2, 0, 2])
#
#     ax_0.set_xticks([0, 50, 100])
#     ax_0.set_xticklabels([])
#
#     ax_0.set_ylabel('Biomass log (g/L)', fontsize=8, family='Arial', )
#     ax_1.set_xlabel('Time (h)', fontsize=8, family='Arial', )
#     ax_1.set_ylabel('Concentration (g/L)', fontsize=8, family='Arial', )
#     if index == 6:
#         legend0 = ax_0.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, prop={'family': 'Arial', 'size': 8}, )
#         legend1 = ax_1.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, prop={'family': 'Arial', 'size': 8}, )
#     ax_0.set_title(titles[index], fontsize=8, family='Arial', )
#     plt.yticks(fontname="Arial", fontsize=8)
#     plt.xticks(fontname="Arial", fontsize=8)
#
#     ax_0.tick_params(axis="x", labelsize=8, pad=-1)
#     ax_0.tick_params(axis="y", labelsize=8, pad=1),
#     # ax_1.tick_params(axis="x", labelsize=8, pad=4)
#     ax_1.tick_params(axis="y", labelsize=8, pad=1),
#     fig.tight_layout(pad=1, h_pad=1, w_pad=0.5)
#     fig.savefig('fig6a_%s.pdf' % titles[index], bbox_inches='tight')
#     fig.savefig('fig6a_%s.jpg' % titles[index], bbox_inches='tight', dpi=600.)
#     fig.show()
