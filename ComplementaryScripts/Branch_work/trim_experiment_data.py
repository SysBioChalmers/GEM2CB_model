#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 7/18/20

"""trim_experiment_data.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

os.chdir('../../ComplementaryData/three_species/experiment_data/')
# %%
trimme_data = False
if trimme_data:

    coloums = ['time', 'fru', 'fru_STDV',
               'R.i', 'R.i_STDV', 'F.p', 'F.p_STDV', 'B.h', 'B.h_STDV',
               'ac', 'ac_STDV', 'for', 'for_STDV', 'but', 'but_STDV', 'lac', 'lac_STDV']
    xls = pd.ExcelFile('experiment_data.xlsx')
    df_list = {}
    for sheet_name in xls.sheet_names:
        df = xls.parse(sheet_name, index_col=0, header=0)
        for biom in ['R.i', 'F.p', 'B.h', ]:
            a = df[biom].values
            new_a = []
            for i in range(0, len(a)):
                if a[i] < 300:
                    new_a.append(a[i])
                else:
                    new_a.append(new_a[-1])
                new_a[-1] = max(new_a)
            df[biom] = np.array(new_a)
        df_list[sheet_name] = df

    with pd.ExcelWriter('experiment_data_trimmed.xlsx') as writer:
        for k, v in df_list.items():
            v.to_excel(writer, sheet_name=k)
# %%
yiels_view = True

if yiels_view:
    data_cloum_name = ['fru', 'R.i', 'F.p', 'B.h', 'ac', 'for', 'but', 'lac']
    xls_2 = pd.ExcelFile('experiment_data_trimmed.xlsx')

    yield_list = {}
    for sheet_name in xls_2.sheet_names:
        df = xls_2.parse(sheet_name, index_col=0, header=0)

        experiment_data_df_trimed = df[data_cloum_name]
        experiment_data_df_trimed_values = experiment_data_df_trimed.values[:, :] - experiment_data_df_trimed.values[0,
                                                                                    :]  # - time0 value
        experiment_data_df_trimed_values = experiment_data_df_trimed_values[1:, :]  # remove time 0
        experiment_data_df_trimed_values = experiment_data_df_trimed_values.T
        experiment_data_df_trimed_values = experiment_data_df_trimed_values[:, :] / abs(
            experiment_data_df_trimed_values[0, :])  # / carbon uptake
        experiment_data_df_trimed_values = experiment_data_df_trimed_values.T
        experiment_data_df_trimed_values = experiment_data_df_trimed_values[:, 1:]  # remove carbon (all are -1 )

        final_yield_row = experiment_data_df_trimed_values[-4:, :]
        final_yield_row = final_yield_row.mean(axis=0)
        yield_list[sheet_name] = final_yield_row
    df = pd.DataFrame.from_dict(yield_list, orient='index',columns=['R.i', 'F.p', 'B.h', 'ac', 'for', 'but', 'lac'])
    df['biom'] = (6-df['ac']*2-df['for']*1-df['but']*4-df['lac']*3)/6
    df.to_csv('yield_summary.csv')
    for i in df.columns:
        plt.scatter(x = df.index,y = df[i])
        plt.title(i)
        plt.xticks(rotation=90)
        # df[i].plot(kind='scatter',x = df.index ,title = i)
        plt.show()
        # plt.savefig(i+'.png')
        plt.clf()

