#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-09-12

"""Step1_iML1515_test.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import Cybernetic_Functions

os.chdir('../Data/')

print('\n---------- Loading experiment data ... ---------- ')
experiment_data_df = pd.read_csv('Case1_ecoli_reduced/ecoli_reduce_experiment_data.txt', delimiter='\t', header=0)

# %% <Cybernetic model simulations>:

tStart = 0.0  # DefineTime
tStop = 8.5
tStep = 0.1
# tspan = np.linspace(tStart, tStop, (tStop - tStart) / tStep)
tspan = np.linspace(tStart, tStop, int(((tStop - tStart) / tStep) + 1))

# matrix Z: reactions x pathways; and Smz metabolites x pathways , S metabolites x reactions : Smz = Sm @ Z

# np.savetxt('Case3_iML1515/iML1515_Smz.txt', Smz, delimiter=',')
Smz = np.loadtxt('Case3_iML1515/iML1515_Smz.txt', delimiter=',')
(n_mets, n_path) = Smz.shape
# experiment data
metabObj = ['glc', 'biomass', 'ac', 'for', 'etoh', 'lac', 'succ']

# initial metabolites at tome 0, t0: initial_mets.shape = (n_mets,)
initial_mets = experiment_data_df[metabObj].values[0, :]

# initial:Enzyme: initial_enzyme.shape = (n_path,)
initial_enzyme = np.array([0.9] * n_path)

# initial data x0 :initial_x0.shape  = (n_mets + n_path,)
initial_x0 = np.concatenate((initial_mets, initial_enzyme))

# Enzyme Rate Parameters: alpha,beta,ke : de/dt =  alpha + rE(ke) * u - (beta + mu) * e
alpha = np.array([0.04] * n_path)
beta = np.array([0.05] * n_path)
ke = np.array([0.620342] * n_path)  # or 0.5

# Metabolites rate Parameters kmax , Ki : dm/dt =  Smz @ V @ rM(kmax,K) * c
# kmax : n_path
kmax = np.array([
    2.8334e-02,
    1.3296e-01,
    3.3288e+00,
    6.7559e+00,
    1.7234e+01,
    3.6402e-01,
    6.9528e+00,
    2.9645e+02
])

# K : n_path
K = np.array([
    2.2462e-04,
    1.8866e+00,
    1.3772e-01,
    2.3646e+01,
    1.0928e+01,
    1.7817e+00,
    5.6565e-04,
    9.0691e+01,
])

# carbon number for each pathways
n_carbon = np.array([6, 6, 6, 6, 6, 6, 6, 0])
Biomass_index = 1
sub_index = 0


class Cybernetic_Model_basic(Cybernetic_Functions.Cybernetic_Model):
    pass


ecoli_iML1515_our_cb = Cybernetic_Model_basic('CB model for Ecoli iML1515 matrix ')

# ecoli_iML1515_our_cb = Cybernetic_Functions.Cybernetic_Model('CB model for Ecoli iML1515 matrix ')
ecoli_iML1515_our_cb.Smz = Smz
ecoli_iML1515_our_cb.x0 = initial_x0
ecoli_iML1515_our_cb.kmax = kmax
ecoli_iML1515_our_cb.mets_name = ['glc', 'biomass', 'ac', 'for', 'etoh', 'lac', 'succ']
ecoli_iML1515_our_cb.K = K
ecoli_iML1515_our_cb.ke = ke
ecoli_iML1515_our_cb.alpha = alpha
ecoli_iML1515_our_cb.beta = beta
ecoli_iML1515_our_cb.n_carbon = n_carbon
ecoli_iML1515_our_cb['sub_index'] = sub_index
ecoli_iML1515_our_cb['Biomass_index'] = Biomass_index
ecoli_iML1515_our_cb.experiment_data_df = experiment_data_df

CB_model = ecoli_iML1515_our_cb


def rate_def(self, x):
    name = self.name
    Smz = self.Smz
    kmax = self.kmax
    K = self.K
    ke = self.ke
    Biomass_index = self.Biomass_index
    sub_index = self.sub_index

    (n_mets, n_path) = Smz.shape

    r_kin_basic = [
        x[sub_index] / (K[0] + x[sub_index]),
        x[sub_index] / (K[1] + x[sub_index]),
        x[sub_index] / (K[2] + x[sub_index]),
        x[sub_index] / (K[3] + x[sub_index]),
        x[sub_index] / (K[4] + x[sub_index]),
        x[sub_index] / (K[5] + x[sub_index]),
        x[sub_index] / (K[6] + x[sub_index]),
        x[3] / (K[7] + x[3])]

    rM = r_kin_basic * x[n_mets:] * kmax

    rE = ke * r_kin_basic

    rG = Smz[Biomass_index, :] * rM[:]  # 11: biomass index

    return rM, rE, rG


def cybernetic_var_def(self, rM):
    # print('def cybernetic_var')
    (n_mets, n_path) = self.Smz.shape
    # cybernetic_var = abs(Smz[sub_index, :] * rM * n_carbon)
    cybernetic_var = rM * self.n_carbon
    # cybernetic_var[cybernetic_var==0] = 1
    cybernetic_var[cybernetic_var < 0] = 0
    if sum(cybernetic_var) > 0:
        u = cybernetic_var / sum(cybernetic_var)
        v = cybernetic_var / np.max(abs(cybernetic_var))
    else:
        u = v = np.zeros(n_path)

    if CB_model.name in ['CB model for Ecoli iML1515 matrix ']:
        # print(123)
        u[-1] = 1.0
        v[-1] = 1.0
    return u, v


setattr(Cybernetic_Model_basic, 'rate_def', rate_def)
setattr(Cybernetic_Model_basic, 'cybernetic_var_def', cybernetic_var_def)

sol = Cybernetic_Functions.cb_model_simulate(CB_model, tspan, draw=True)
# np.savetxt('Case3_iML1515/iML1515_sol.csv', sol, delimiter=',')
# sol_ = np.loadtxt('Case3_iML1515/iML1515_sol.csv', delimiter=',')

# %% <fitting>
minimums = []
for i in [0, 1, 2, 3, 4, 5, 6, 7]:
    print(i)
    CB_model.weights_of_time = [0, 1]
    CB_model.weights_of_mets = np.array([1, 2, 1, 1, 1, 1, 1, ])
    CB_model.weight_of_method_t_f = i
    # cb_model_1.weights_of_mets = weights_of_mets/1000
    cb_models = [CB_model, ]  # [1], [2], [3],[4],[5], [6], [7],
    experiment_data_dfs = [experiment_data_df, ]
    para_to_fit = {'kmax': np.array([[0], [1], [2], [3], [4], [5], [6], [7]]),
                   'K': np.array([[0], [1], [2], [3], [4], [5], [6], [7]])}

    minimum = Cybernetic_Functions.parameters_fitting(cb_models, experiment_data_dfs,
                                                      para_to_fit, tspan, draw=False,
                                                      full_output=False,
                                                      options={'xtol': 0.5, 'ftol': 0.1, 'maxiter': 100,
                                                               'disp': True})

    minimums.append(minimum)

for minimum in minimums:
    CB_models = Cybernetic_Functions.update_paras_func(minimum.x, [CB_model], para_to_fit,
                                                       retern_model_or_paras='model')
    sol3 = Cybernetic_Functions.cb_model_simulate(CB_models[0], tspan, draw=False)


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
            ax1.plot(tspan, sol3[:, index], color=color_list[index], linewidth=1, label=labels[index])

            experiment_p = ax1.plot(experiment_data_df['time'], experiment_data_df[metabObj[index]], '.',
                                    color=color_list[index], alpha=0.8, markersize=5)

        else:
            ax.plot(tspan, sol3[:, index], color=color_list[index], linewidth=1, label=labels[index])

            experiment_p = ax.plot(experiment_data_df['time'], experiment_data_df[metabObj[index]], '.',
                                   color=color_list[index], alpha=0.8, markersize=5)

    ax.set_xlabel('Time (h)', fontsize=10, family='Arial', )
    ax.set_ylabel('Concentration (mM)', fontsize=10, family='Arial', )
    ax1.set_ylabel('Biomass (g/L)', fontsize=10, family='Arial', )
    ax1.set_ylim(0,2)
    plt.yticks(fontname="Arial", fontsize=8)
    plt.xticks(fontname="Arial", fontsize=8)
    L = fig.legend(loc='lower left', bbox_to_anchor=(1, 0.21), ncol=1, fontsize=8, )
    plt.setp(L.texts, family='Arial')
    fig.tight_layout()


    fig.show()



# %% <plot cybernetic model result>

# experiment data
fig = plt.figure()
ax = fig.add_subplot(111)
colors = ['blue', 'teal', 'tab:red', 'tab:orange']
color_list = plt.cm.tab10(np.linspace(0, 1, 12))

for index in range(0, CB_model.Smz.shape[0]):
    ax.plot(tspan, sol[:, index], color=color_list[index + 1], linewidth=2, label=metabObj[index])
    experiment_p = ax.plot(experiment_data_df['time'], experiment_data_df[metabObj[index]], 'o--',
                           color=color_list[index + 1],
                           linewidth=1)

ax.set_xlabel('Time (h)', fontsize=12)
ax.set_ylabel('Concentration (mM)', fontsize=12)
fig.legend(loc='center left', bbox_to_anchor=(0.1, 0.57), ncol=2)

fig.show()



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
ax1.set_ylim(0,2)
plt.yticks(fontname="Arial", fontsize=8)
plt.xticks(fontname="Arial", fontsize=8)
L = fig.legend(loc='lower left', bbox_to_anchor=(1, 0.21), ncol=1, fontsize=8, )
plt.setp(L.texts, family='Arial')
fig.tight_layout()


fig.show()
