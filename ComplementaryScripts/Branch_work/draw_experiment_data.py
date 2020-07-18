#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 7/18/20

"""draw_experiment_data.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os

import pandas as pd
import matplotlib.pyplot as plt

os.chdir('../../ComplementaryData/three_species/experiment_data/')

coloums = ['time', 'fru', 'fru_STDV',
           'R.i', 'R.i_STDV', 'F.p', 'F.p_STDV', 'B.h', 'B.h_STDV',
           'ac', 'ac_STDV', 'for', 'for_STDV', 'but', 'but_STDV', 'lac', 'lac_STDV']

xls = pd.ExcelFile('experiment_data_trimmed.xlsx')
# color_list = plt.cm.tab10(np.linspace(0, 1, 12))
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

for sheet_name in xls.sheet_names:
    df = xls.parse(sheet_name, index_col=0, header=0)

    figsize = (8, 6)
    # fig = plt.figure(figsize=(6, 2.5))
    fig, axs = plt.subplots(3, 1, figsize=figsize)
    # ax_1 = fig.add_subplot(111)
    ax_0 = axs[0]
    ax_1 = axs[2]
    ax_2 = axs[1]

    mets_name = ['fru',
                 'R.i', 'F.p', 'B.h',
                 'ac', 'for', 'but', 'lac', ]

    tspan = df['time']
    color_i = 0
    for met_i in mets_name:
        if met_i in ['R.i', 'F.p', 'B.h']:
            ax_0.errorbar(tspan, df[met_i], fmt='o--', yerr=df[met_i + '_STDV'], color=colors[color_i], label=met_i)
            ax_2.errorbar(tspan, df[met_i] * 1e8, fmt='o--', yerr=df[met_i + '_STDV'] * 1e8, color=colors[color_i],
                          label=met_i)
            ax_2.set_yscale('log')
        else:
            ax_1.errorbar(tspan, df[met_i], fmt='o--', yerr=df[met_i + '_STDV'], color=colors[color_i], label=met_i)
        color_i += 1

    ax_0.axis(ymin=0.01, ymax=150)
    ax_2.axis(ymin=1e4, ymax=1e12)
    ax_1.axis(ymin=0.01, ymax=80)
    ax_0.set_ylabel('Biomass (g/L)', fontsize=12)
    ax_1.set_xlabel('Time (h)', fontsize=12)
    ax_1.set_ylabel('Concentration (g/L)', fontsize=12)

    ax_0.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)
    ax_1.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)
    ax_0.set_title(sheet_name)
    # plt.title(titles[index])
    plt.show()
    fig.savefig(sheet_name + '_trimmed.png')
