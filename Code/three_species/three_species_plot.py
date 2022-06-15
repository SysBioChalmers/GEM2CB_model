#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 7/29/20

"""draw_results.py
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
import matplotlib
import seaborn as sns

matplotlib.rc('font', family="Arial")
matplotlib.rcParams["font.family"] = 'Arial'  # 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']

os.chdir('../../Data/three_species/')
# %% <experiment data> experiment data ti get initial mets
print('\n---------- Loading Experiment Data ... ---------- ')
file_name = 'experiment_data/experiment_data.xlsx'
RI_experimet_data = pd.read_excel(file_name, sheet_name='RI_14', index_col=0, header=0, )
FP_experimet_data = pd.read_excel(file_name, sheet_name='FP_4', index_col=0, header=0, )
BH_experimet_data = pd.read_excel(file_name, sheet_name='BH_14', index_col=0, header=0, )

RI_FP_experimet_data = pd.read_excel(file_name, sheet_name='RI_FP_8', index_col=0, header=0, )
RI_BH_experimet_data = pd.read_excel(file_name, sheet_name='RI_BH_4', index_col=0, header=0, )
FP_BH_experimet_data = pd.read_excel(file_name, sheet_name='FP_BH_1', index_col=0, header=0, )

RI_BH_experimet_data_no_ac = pd.read_excel(file_name, sheet_name='RI_BH_7', index_col=0, header=0, )
FP_BH_experimet_data_no_ac = pd.read_excel(file_name, sheet_name='FP_BH_3', index_col=0, header=0, )

RI_FP_BH_10_experimet_data = pd.read_excel(file_name, sheet_name='RI_FP_BH_10', index_col=0, header=0, )
RI_FP_BH_12_experimet_data = pd.read_excel(file_name, sheet_name='RI_FP_BH_12', index_col=0, header=0, )

# %% <simulation result>
tStart = 0.0  # DefineTime
tStop = 50
tStep = 0.5
tspan = np.linspace(tStart, tStop, int(((tStop - tStart) / tStep) + 1))

name_list = ['mono', 'bi_ac', 'tri_1', 'tri_2', 'mono_bi', 'mono_bi_tri_1', 'mono_bi_tri_2']

simulation_dic = {}
for name in name_list:
    simulation_file_name = name
    with open('result_no_constraint/' + simulation_file_name + '.pickle', 'rb') as file:
        simulation = pickle.load(file)
        # simulation = [sol1, sol2, sol3, sol12, sol13, sol23, sol13_no_ac, sol23_no_ac, sol123, sol123_2]
    # print(simulation)
    simulation_dic[name] = simulation

experiment_data_dfs = [RI_experimet_data, FP_experimet_data, BH_experimet_data] + \
                      [RI_FP_experimet_data, RI_BH_experimet_data, FP_BH_experimet_data] + \
                      [RI_FP_BH_10_experimet_data, RI_FP_BH_12_experimet_data, ]
biomass_ids = [[1], [2], [3],
               [1, 2], [1, 3], [2, 3],
               [1, 2, 3], [1, 2, 3]]
titles = ['$R.i$', '$F.p$', '$B.h$',
          '$R.i$ & $F.p$', '$R.i$ & $B.h$', '$F.p$ & $B.h$',
          '$R.i$ & $F.p$ & $B.h$', '$R.i$ & $F.p$ & $B.h$']
# %%
color_list = plt.cm.tab10(np.linspace(0, 1, 12))
# color_list = sns.color_palette('muted')
# del color_list[0]
# colors = sns.color_palette('dark')
sns.set_style("ticks", )
mets_name = ['fru', 'R.i', 'F.p', 'B.h', 'ac', 'for', 'but']
labels = ['Fructose', '$R.i$', '$F.p$', '$B.h$', 'Acetate', 'Formate', 'Butyrate']
for k, v in simulation_dic.items():
    # %%
    # k = 'mono_bi'
    # v = simulation_dic[k]

    simulation = v

    for i in range(0, 8):
        simulation_i = simulation[i]
        experiment_data_df = experiment_data_dfs[i]
        if i in {0, 3}:
            figsize = (8, 2.5)
            fig, axs = plt.subplots(2, 3, figsize=figsize)
        elif i in {6}:
            figsize = (6, 2.5)
            fig, axs = plt.subplots(2, 2, figsize=figsize)

        if i < 3:
            ax_0 = axs[0, i]
            ax_1 = axs[1, i]
        elif 2 < i < 6:
            ax_0 = axs[0, i - 3]
            ax_1 = axs[1, i - 3]
        else:
            simulation_i = simulation[i + 2]
            ax_0 = axs[0, i - 6]
            ax_1 = axs[1, i - 6]

        for index in range(0, len(mets_name)):
            if index in biomass_ids[i]:
                ax_0.plot(tspan, simulation_i[:, index], color=color_list[index + 1], linewidth=1,
                          label=labels[index])
                experiment_p = ax_0.plot(experiment_data_df['time'], experiment_data_df[mets_name[index]], '.',
                                         color=color_list[index + 1],
                                         markersize=5)

            if index not in [1, 2, 3]:
                ax_1.plot(tspan, simulation_i[:, index], color=color_list[index + 1], linewidth=1,
                          label=labels[index])
                experiment_p = ax_1.plot(experiment_data_df['time'], experiment_data_df[mets_name[index]], '.',
                                         color=color_list[index + 1],
                                         markersize=5)
        ax_0.set_ylim((0, 150))
        ax_1.set_ylim((0, 100))
        ax_0.tick_params(axis="x", labelsize=8, pad=1)
        ax_0.tick_params(axis="y", labelsize=8, pad=1),
        ax_1.tick_params(axis="x", labelsize=8, pad=1)
        ax_1.tick_params(axis="y", labelsize=8, pad=1),
        ax_1.set_xlabel('Time (h)', fontsize=8, family='Arial', )
        ax_0.set_title(titles[i], fontsize=8, family='Arial', )

        if i in {0, 3, 6}:
            ax_0.set_ylabel('Counts ($10^8/$mL)', fontsize=8, family='Arial', )
            ax_1.set_ylabel('Concentration (mM)', fontsize=8, family='Arial', )
        if i in {7}:
            ax_0.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, prop={'family': 'Arial', 'size': 8}, )
            ax_1.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, prop={'family': 'Arial', 'size': 8}, )

        if i in {2, 5, 7}:
            fig.tight_layout(pad=1, h_pad=0.5, w_pad=0.5)
            fig.savefig('result_no_constraint/fig7_' + k + titles[i] + '_.pdf', bbox_inches='tight')
            fig.savefig('result_no_constraint/fig7_' + k + titles[i] + '_.jpg', bbox_inches='tight')
            fig.show()
