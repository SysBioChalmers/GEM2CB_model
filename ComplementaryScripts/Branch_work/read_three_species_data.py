#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 7/17/20

"""read_three_species_data.py
:description : script
:param : 
:returns: 
:rtype: 
"""
import os

import pandas as pd

os.chdir('../../ComplementaryData/three_species/initial_data/Fermenter/')

coloums = ['time', 'fru', 'fru_STDV',
           'R.i', 'R.i_STDV', 'F.p', 'F.p_STDV', 'B.h', 'B.h_STDV',
           'ac', 'ac_STDV', 'for', 'for_STDV', 'but', 'but_STDV', 'lac', 'lac_STDV']

co_list = ['RI_FP_8', 'RI_FP_9', 'RI_BH_4', 'RI_BH_6', 'RI_BH_5', 'RI_BH_7', 'FP_BH_1', 'FP_BH_2', 'FP_BH_3']
tri_list = ['RI_FP_BH_10', 'RI_FP_BH_11', 'RI_FP_BH_12', 'RI_FP_BH_13', 'RI_FP_BH_14', 'RI_FP_BH_15']
mono_list = ['RI_8', 'RI_14', 'RI_15', 'RI_16', 'FP_4', 'FP_14', 'FP_15', 'BH_14', 'BH_15', 'BH_16']

FP_sheet_set = set(['Determination cell counts FP', 'Determination cell count F. p'])
BH_sheet_set = set(['Determination cell counts BH', 'Determination cell count BH'])
mono_sheet_set = set(['Determination cell count GOOD', 'Determination cell count'
                         , 'Determination cell counts RI', 'Determination cell counts'])
df_list = {}

for file_name in mono_list + co_list + tri_list:
    print(file_name)
    meta_df = pd.read_excel(file_name + '.xlsx', sheet_name='Metabolites', header=None, usecols='D:U')
    # print(meta_df[[4]].iloc[[2,20]])
    # print(meta_df[[4]].iloc[[23,41]])
    range_index_1 = range(3, 20)
    range_index_2 = range(24, 41)
    if file_name in ['RI_8']:
        range_index_1 = range(3, 19)
        range_index_2 = range(23, 39)
    if file_name in ['FP_4', 'FP_14', 'FP_15', 'BH_15', 'BH_16']:
        range_index_1 = range(3, 21)
        range_index_2 = range(25, 43)

    df_i = pd.DataFrame(columns=coloums)
    df_i[['time', 'fru', 'fru_STDV', 'ac', 'ac_STDV', 'for', 'for_STDV', 'lac', 'lac_STDV']] = \
        meta_df[[4, 5, 6, 13, 14, 17, 18, 9, 10]].iloc[range_index_1]
    df_i[['but', 'but_STDV']] = meta_df[[9, 10]].iloc[range_index_2].values

    xl = pd.ExcelFile(file_name + '.xlsx', )

    if 'RI' in file_name and file_name not in mono_list:
        sheet_name = 'Determination cell counts RI'
        cell_df_RI = pd.read_excel(file_name + '.xlsx', sheet_name=sheet_name, header=None,
                                   usecols='P:Q')
        df_i[['R.i', 'R.i_STDV', ]] = cell_df_RI.iloc[range_index_1].values / 1e8

    if 'FP' in file_name and file_name not in mono_list:
        sheet_name = (FP_sheet_set & set(xl.sheet_names)).pop()
        cell_df_FP = pd.read_excel(file_name + '.xlsx', sheet_name=sheet_name, header=None,
                                   usecols='P:Q')
        df_i[['F.p', 'F.p_STDV', ]] = cell_df_FP.iloc[range_index_1].values / 1e8
    if 'BH' in file_name and file_name not in mono_list:
        sheet_name = (BH_sheet_set & set(xl.sheet_names)).pop()
        cell_df_BH = pd.read_excel(file_name + '.xlsx', sheet_name=sheet_name, header=None,
                                   usecols='P:Q')
        df_i[['B.h', 'B.h_STDV', ]] = cell_df_BH.iloc[range_index_1].values / 1e8

    if file_name in mono_list:
        sheet_name = (mono_sheet_set & set(xl.sheet_names)).pop()
        cell_df = pd.read_excel(file_name + '.xlsx', sheet_name=sheet_name, header=None,
                                usecols='P:Q')
        if 'RI' in file_name:
            df_i[['R.i', 'R.i_STDV', ]] = cell_df.iloc[range_index_1].values / 1e8
        elif 'FP' in file_name:
            df_i[['F.p', 'F.p_STDV', ]] = cell_df.iloc[range_index_1].values / 1e8
        else:
            df_i[['B.h', 'B.h_STDV', ]] = cell_df.iloc[range_index_1].values / 1e8

    df_i = df_i.fillna(0)

    df_i.to_csv('../../experiment_data/' + file_name + '_.csv', sep='\t')
    df_list[file_name] = df_i

with pd.ExcelWriter('../../experiment_data/experiment_data.xlsx') as writer:
    for k, v in df_list.items():
        v.to_excel(writer, sheet_name=k)
