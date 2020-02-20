#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2/19/20

"""bi_Ri_Bh_ac.py
:description : script
:param : 
:returns: 
:rtype: 
"""
import ConvexHull_yield
import GEM2pathways
import Cybernetic_Functions
import cobra
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import gemstool
from scipy.spatial import ConvexHull
import os

os.chdir('../../ComplementaryData/three_species/')

Smz_Ri = np.loadtxt('Smz_Rint.txt', delimiter=',')  # 6 x 7
Smz_Bh = np.loadtxt('Smz_Bhyd.txt', delimiter=',')  # 5 x 6
Ri_metas = ['fru', 'R.i', 'ac', 'for', 'but', 'lac']
Bh_metas = ['fru', 'B.h', 'ac', 'for', 'lac']

bio_met = np.array([0, 0, 0, 0, 0, 0, 0])
Smz_Ri = np.insert(Smz_Ri, 2, values=bio_met, axis=0)
Smz_Ri[4, :] = Smz_Ri[4, :]
Ri_metas = ['fru', 'R.i', 'B.h', 'ac', 'for', 'but', 'lac']

but_met = np.array([0, 0, 0, 0, 0, 0])
Smz_Bh = np.insert(Smz_Bh, 4, values=but_met, axis=0)
bio_met = np.array([0, 0, 0, 0, 0, 0])
Smz_Bh = np.insert(Smz_Bh, 1, values=but_met, axis=0)
Bh_metas = ['fru', 'R.i', 'B.h', 'ac', 'for', 'but', 'lac']

metas_name = Bh_metas
Smz = np.concatenate((Smz_Ri, Smz_Bh), axis=1)
np.savetxt('Smz_bi_Ri_Bh.txt', Smz, delimiter=',')
# metabolites and pathways number
(n_mets, n_path) = Smz.shape
# matrix Z: reactions x pathways; and Smz metabolites x pathways , S metabolites x reactions : Smz = Sm @ Z

np.set_printoptions(suppress=True)

print('CB modeing')
tStart = 0.0  # DefineTime
tStop = 30
tStep = 0.5
tspan = np.linspace(tStart, tStop, (tStop - tStart) / tStep)

experiment_data_df = pd.read_csv('bi_Ri_Bh_experiment_data_ac.txt', delimiter='\t', header=0)
experiment_data_df['R.i'] = experiment_data_df['R.i'] * 5e-3
experiment_data_df['B.h'] = experiment_data_df['B.h'] * 5e-3

# initial metabolites at tome 0, t0: initial_mets.shape = (n_mets,)
initial_mets = experiment_data_df[metas_name].values[0, :]
# initial_mets = [50,	6.929513*5e-3, 1,	50,	0.754884547]
initial_enzyme = np.array([0.9] * n_path)

# initial data x0 :initial_x0.shape  = (n_mets + n_path,)
initial_x0 = np.concatenate((initial_mets, initial_enzyme))

# Enzyme Rate Parameters: alpha,beta,ke : de/dt =  alpha + rE(ke) * u - (beta + mu) * e
alpha = np.array([0.04] * n_path)
beta = np.array([0.05] * n_path)
ke = np.array([0.5] * n_path)  # or 0.5

# Metabolites rate Parameters kmax , Ki : dm/dt =  Smz @ V @ rM(kmax,K) * c
# Ri + Bh
kmax = np.array([
    16,
    14,
    15,
    30,
    40,
    40,
    15,

    10,
    11,
    10,
    10.4,
    10,
    25,
])

# # K : n_path
K = np.array([
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
])

# carbon number for each pathways
n_carbon = np.array([6, 6, 6, 6, 6, 6, 2] + [6, 6, 6, 6, 6, 1])
Biomass_index = 1
sub_index = 0

biRiBh_cb = Cybernetic_Functions.Cybernetic_Model('CB model for bi R.int and B.hyd')
biRiBh_cb.Smz = Smz
biRiBh_cb.x0 = initial_x0
biRiBh_cb.kmax = kmax
biRiBh_cb.K = K
biRiBh_cb.ke = ke
biRiBh_cb.alpha = alpha
biRiBh_cb.beta = beta
biRiBh_cb.n_carbon = n_carbon
biRiBh_cb['sub_index'] = sub_index
biRiBh_cb['Biomass_index'] = Biomass_index

CB_model = biRiBh_cb


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
        x[sub_index] / (K[0] + x[sub_index]),
        x[sub_index] / (K[1] + x[sub_index]),
        x[sub_index] / (K[2] + x[sub_index]),
        x[sub_index] / (K[3] + x[sub_index]),
        x[sub_index] / (K[4] + x[sub_index]),
        x[sub_index] / (K[5] + x[sub_index]),
        x[sub_index] / (K[6] + x[sub_index]),

        x[sub_index + 1] / (K[7] + x[sub_index + 1]),
        x[sub_index + 1] / (K[8] + x[sub_index + 1]),
        x[sub_index + 1] / (K[9] + x[sub_index + 1]),
        x[sub_index + 1] / (K[10] + x[sub_index + 1]),
        x[sub_index + 1] / (K[11] + x[sub_index + 1]),
        x[sub_index + 1] / (K[12] + x[sub_index + 1]),

    ]

    rM = r_kin_basic * x[n_mets:] * kmax

    rE = ke * r_kin_basic

    rG = Smz[Biomass_index, :] * rM[:]  # + Smz[Biomass_index+1, :] * rM[:] # 11: biomass index

    return rM, rE, rG


def cybernetic_var_def(rM, CB_model):
    # print('def cybernetic_var')
    (n_mets, n_path) = CB_model.Smz.shape
    # cybernetic_var = abs(Smz[sub_index, :] * rM * n_carbon)
    cybernetic_var = rM * CB_model.n_carbon
    cybernetic_var[cybernetic_var < 0] = 0
    cybernetic_var_1 = cybernetic_var[0:7]
    cybernetic_var_2 = cybernetic_var[7:]
    u = v = np.array([1] * n_path)

    if sum(cybernetic_var_1) > 0:
        u[0:7] = cybernetic_var_1 / sum(cybernetic_var_1)
        v[0:7] = cybernetic_var_1 / np.max(abs(cybernetic_var_1))
    else:
        u[0:7] = v[0:7] = np.zeros(n_path)

    if sum(cybernetic_var_2) > 0:
        u[7:] = cybernetic_var_2 / sum(cybernetic_var_2)
        v[7:] = cybernetic_var_2 / np.max(abs(cybernetic_var_2))
    else:
        u[7:] = v[7:] = np.zeros(n_path)

    if True:
        u[6] = 1.0
        v[6] = 1.0
        u[-1] = 1.0
        v[-1] = 1.0
    return u, v


sol = Cybernetic_Functions.cb_model_simulate(CB_model, tspan, draw=True)
metabObj = metas_name

# %% < Fiting >
CB_model['metas_names'] = metas_name
para_to_fit = {'kmax': [0, 1, 2, 3, 4, 5], 'K': []}

# experiment data
experiment_data_df = pd.read_csv('bi_Ri_Bh_experiment_data_ac.txt', delimiter='\t', header=0)
experiment_data_df['R.i'] = experiment_data_df['R.i'] * 5e-3
experiment_data_df['B.h'] = experiment_data_df['B.h'] * 5e-3

experiment_data_to_fit = experiment_data_df[['time'] + metabObj].iloc[0:-1]
# minimum = Cybernetic_Functions.parameters_fitting(CB_model, experiment_data_to_fit, para_to_fit, tspan, draw=True)


# %% <plot cybernetic model result>
# experiment data

fig = plt.figure(figsize=(6, 2.5))
ax = fig.add_subplot(111)
colors = ['blue', 'teal', 'tab:red', 'tab:orange']
color_list = plt.cm.tab10(np.linspace(0, 1, 12))
experiment_data_df = pd.read_csv('bi_Ri_Bh_experiment_data_ac.txt', delimiter='\t', header=0)
experiment_data_df['R.i'] = experiment_data_df['R.i'] * 5e-3
experiment_data_df['B.h'] = experiment_data_df['B.h'] * 5e-3
experiment_data_df = experiment_data_df[0:-1]
for index in range(0, CB_model.Smz.shape[0]):
    ax.plot(tspan, sol[:, index], color=color_list[index + 1], linewidth=2, label=metabObj[index])
    experiment_p = ax.plot(experiment_data_df['time'], experiment_data_df[metabObj[index]], 'o--',
                           color=color_list[index + 1],
                           linewidth=1)

ax.set_xlabel('Time (h)', fontsize=12)
ax.set_ylabel('Concentration (mM)', fontsize=12)
# box = ax.get_position()
# ax.set_position([box.x0, box.y0, box.width , box.height* 0.8])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)

# fig.savefig('test.png')
fig.show()
