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

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

os.chdir('../../ComplementaryData/three_species/experiment_data/')

coloums = ['time', 'fru', 'fru_STDV',
           'R.i', 'R.i_STDV', 'F.p', 'F.p_STDV', 'B.h', 'B.h_STDV',
           'ac', 'ac_STDV', 'for', 'for_STDV', 'but', 'but_STDV', 'lac', 'lac_STDV']
df_list = {}

xls = pd.ExcelFile('experiment_data.xlsx')
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
