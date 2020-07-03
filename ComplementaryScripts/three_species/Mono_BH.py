#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 5/5/20

"""Branch_work.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os

import numpy as np
import pandas as pd

os.chdir('../../ComplementaryData/three_species/')


def yield_view(data_df, carbon_s='fru', produces=['ac', 'for', 'but', 'lac']):
    basic_s = data_df[carbon_s].values[-1] - data_df[carbon_s].values[0]
    yiedd_list = [-1]
    for i in produces:
        yield_i = (data_df[i].values[-1] - data_df[i].values[0]) / (-basic_s)
        yiedd_list.append(yield_i)

    return yiedd_list


RI_experimet_data = pd.read_excel('experiment_data.xlsx', sheet_name='RI_14', header=0, usecols='A:I')
FP_experimet_data = pd.read_excel('experiment_data.xlsx', sheet_name='FP_4', header=0, usecols='A:I')
BH_experimet_data = pd.read_excel('experiment_data.xlsx', sheet_name='BH_14', header=0, usecols='A:I')

RI_BH_experimet_data = pd.read_excel('experiment_data.xlsx', sheet_name='RI_BH_4', header=0, usecols='A:I')
RI_FP_experimet_data = pd.read_excel('experiment_data.xlsx', sheet_name='RI_FP_8', header=0, usecols='A:I')
FP_BH_experimet_data = pd.read_excel('experiment_data.xlsx', sheet_name='FP_BH_1', header=0, usecols='A:I')

experimet_data = BH_experimet_data
cell_number = np.ptp(experimet_data['B.h']) / (1e8)
# 20.7
convert_rate = 0.877 * 15.5 / 10 / 20  # cell_number = 20.7

experimet_data['B.h'] = experimet_data['B.h'] * convert_rate / 1e8

metabObj = ['fru', 'R.i', 'F.p', 'B.h', 'ac', 'for', 'but', 'lac']

yiedd_list = yield_view(experimet_data)
Smz = np.array([-1, 0, 0.0, 0.0887, 1.2, -2.48, 0, 0.38]).reshape(8, 1)

print('CB modeing')
tStart = 0.0  # DefineTime
tStop = 50
tStep = 0.5
tspan = np.linspace(tStart, tStop, (tStop - tStart) / tStep)

(n_mets, n_path) = Smz.shape
experiment_data_df = experimet_data[metabObj]
initial_mets = experiment_data_df.values[0, :]
# initial_mets[0] = initial_mets[0]-37
initial_mets[0] = initial_mets[0]
initial_enzyme = np.array([0.9] * n_path)
initial_x0 = np.concatenate((initial_mets, initial_enzyme))

# Enzyme Rate Parameters: alpha,beta,ke : de/dt =  alpha + rE(ke) * u - (beta + mu) * e
alpha = np.array([0.04] * n_path)
beta = np.array([0.05] * n_path)
ke = np.array([0.5] * n_path)  # or 0.5

kmax = np.array([
    7,
])

K = np.array([
    10,
])

# carbon number for each pathways
n_carbon = np.array([6, ])
Biomass_index = 3
sub_index = 0
import Cybernetic_Functions

CB_model = Cybernetic_Functions.Cybernetic_Model('CB model')
CB_model.Smz = Smz
CB_model.x0 = initial_x0
CB_model.kmax = kmax
CB_model.K = K
CB_model.ke = ke
CB_model.alpha = alpha
CB_model.beta = beta
CB_model.n_carbon = n_carbon
CB_model['sub_index'] = sub_index
CB_model['Biomass_index'] = Biomass_index


def rate_def(x, CB_model):
    Smz = CB_model.Smz
    kmax = CB_model.kmax
    K = CB_model.K
    ke = CB_model.ke
    alpha = CB_model.alpha
    beta = CB_model.beta
    Biomass_index = CB_model.Biomass_index
    sub_index = CB_model['sub_index']

    (n_mets, n_path) = Smz.shape

    r_kin_basic = [
        (x[sub_index] - 37) / (K[0] + x[sub_index] - 37),
    ]

    rM = r_kin_basic * x[n_mets:] * kmax

    rE = ke * r_kin_basic

    rG = Smz[Biomass_index, :] * rM[:]  # 11: biomass index

    return rM, rE, rG


def cybernetic_var_def(rM, CB_model):
    # print('def cybernetic_var')
    (n_mets, n_path) = CB_model.Smz.shape
    # cybernetic_var = abs(Smz[sub_index, :] * rM * n_carbon)
    cybernetic_var = rM * CB_model.n_carbon
    # cybernetic_var[cybernetic_var==0] = 1
    cybernetic_var[cybernetic_var < 0] = 0
    if sum(cybernetic_var) > 0:
        u = cybernetic_var / sum(cybernetic_var)
        v = cybernetic_var / np.max(abs(cybernetic_var))
    else:
        u = v = np.zeros(n_path)

    if CB_model.name in ['CB model for Rint', 'CB model for Ecoli reduced matrix ',
                         'CB model for Ecoli iML1515 matrix ']:
        u[-1] = 1.0
        v[-1] = 1.0
    return u, v


sol = Cybernetic_Functions.cb_model_simulate(CB_model, tspan, draw=True)
np.savetxt('sol_BH.txt', sol, delimiter=',')
# # %%
# CB_model['metas_names'] = ['fru','biomass', 'ac', 'for', 'but','lac' ]
# para_to_fit = {'kmax': [0, 1, 2, 3, 4,5,6], 'K': [0, 1, 2, 3, 4,5,6]}
#
# # experiment data
# experiment_data_to_fit = experiment_data_df[['time']+metabObj].iloc[0:-1]
# experiment_data_to_fit['biomass'] = experiment_data_to_fit['biomass'] * 1e-9
# # minimum = Cybernetic_Functions.parameters_fitting(CB_model, experiment_data_to_fit, para_to_fit, tspan, draw=True)
#
# # %% <plot cybernetic model result>


experiment_data_df = experimet_data
experiment_data_df['B.h'] = experiment_data_df['B.h'] / convert_rate
simulation = sol
simulation[:, 3] = simulation[:, 3] / convert_rate
# simulation[:, 0] = simulation[:, 0] + 37
biomassindexs = [3]

import matplotlib.pyplot as plt

figsize = (8, 4)
# fig = plt.figure(figsize=(6, 2.5))
fig, axs = plt.subplots(2, 1, figsize=figsize)
# ax_1 = fig.add_subplot(111)
ax_1 = axs[1]
ax_0 = axs[0]
colors = ['blue', 'teal', 'tab:red', 'tab:orange']
color_list = plt.cm.tab10(np.linspace(0, 1, 12))

metabObj = ['fru', 'R.i', 'F.p', 'B.h', 'ac', 'for', 'but', 'lac']

for index in range(0, CB_model.Smz.shape[0]):
    if index in biomassindexs:
        ax_0.plot(tspan, simulation[:, index], color=color_list[index + 1], linewidth=1, label=metabObj[index])
        experiment_p = ax_0.plot(experiment_data_df['time'], experiment_data_df[metabObj[index]], 'o',
                                 color=color_list[index + 1],
                                 linewidth=0.5)

    if index not in [1, 2, 3]:
        ax_1.plot(tspan, simulation[:, index], color=color_list[index + 1], linewidth=1, label=metabObj[index])
        experiment_p = ax_1.plot(experiment_data_df['time'], experiment_data_df[metabObj[index]], 'o',
                                 color=color_list[index + 1],
                                 linewidth=0.5)
ax_0.axis(ymin=-5, ymax=100)
ax_0.set_ylabel('Counts (1e8/ml)', fontsize=12)
ax_1.set_xlabel('Time (h)', fontsize=12)
ax_1.set_ylabel('Concentration (mM)', fontsize=12)

ax_0.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)
ax_1.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)

# fig.savefig('Branch_work.png')
fig.show()
