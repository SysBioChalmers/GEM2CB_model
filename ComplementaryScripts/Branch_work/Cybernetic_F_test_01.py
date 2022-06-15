#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 6/12/20

"""Cybernetic_F_test_01.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import copy
import os

import matplotlib
import numpy as np
import pandas as pd

import ComplementaryScripts.Cybernetic_Functions as Cybernetic_Functions

matplotlib.use("macosx")

# %% < def input:>
os.chdir('../../ComplementaryData/')

tStart = 0.0
tStop = 15
tStep = 0.1
tspan = np.linspace(tStart, tStop, int(((tStop - tStart) / tStep) + 1))

# matrix Z: reactions x pathways; and Smz metabolites x pathways , S metabolites x reactions : Smz = Sm @ Z
print('\n---------- Loading FBA modes ... ---------- ')
# Z = np.genfromtxt('Case1_ecoli_reduced/FBA_em_reduced.csv', delimiter=',')
# Z = Z.T
# Smz = Z
# Smz[6, :] = Smz[6, :] - Smz[10, :]
# Smz = Smz[[0, 5, 6, 7, 8, 9, 11], :]
# Smz[0, :] = -Smz[0, :]
Smz = np.array([[-1., -1., -1., -1., -1., 0.],
                [0.01123911, 0.02796746, 0.02796746, 0.01960329, 0.02294896, 0.],
                [0., 0.8673643, 0.8673643, 0.01980717, 0.89116456, 0.],
                [1., 1.62104086, 0., 0.01980717, 1.15071905, -1.],
                [1., 0.75367656, 0.75367656, 0., 0.25955448, 0.],
                [0., 0., 0., 1.70458824, 0., 0.],
                [0.57722474, 0., 0., 0., 0.53832256, 0.]])

# metabolites and pathways number
(n_mets, n_path) = Smz.shape

# experiment data
print('\n---------- Loading Experiment Data ... ---------- ')
experiment_data_df = pd.read_csv('Case1_ecoli_reduced/ecoli_reduce_experiment_data.txt', delimiter='\t',
                                 header=0)
# metabObj = ['glc', 'succ', 'for', 'lac', 'ac', 'etoh', 'biomass', ]
metabObj = ['glc', 'biomass', 'ac', 'for', 'etoh', 'lac', 'succ', ]

# initial metabolites at tome 0, t0: initial_mets.shape = (n_mets,)
initial_mets = experiment_data_df[metabObj].values[0, :]
# initial_mets = initial_mets[[0, 6, 4, 2, 5, 3, 1]]

# initial:Enzyme: initial_enzyme.shape = (n_path,)
initial_enzyme = np.array([0.9, 0.9, 0.9, 0.9, 0.9, 1])

# initial data x0 :initial_x0.shape  = (n_mets + n_path,)
initial_x0 = np.concatenate((initial_mets, initial_enzyme))

# Enzyme Rate Parameters: alpha,beta,ke : de/dt =  alpha + rE(ke) * u - (beta + mu) * e
alpha = np.array([0.04] * n_path)
beta = np.array([0.05] * n_path)
ke = np.array([0.620342] * n_path)  # or 0.5

# Metabolites rate Parameters kmax , Ki : dm/dt =  Smz @ V @ rM(kmax,K) * c
# TODO the equations not the same , should defined by user

# kmax : n_path

kmax = np.array([
    10,
    35.17,
    4.1210e-09,
    2.31,
    4.12e-7,
    21.08
])

# K : n_path
K = np.array([
    2.5235,
    24.39,
    33.50,
    1,
    1,
    1.042,
])

# carbon number for each pathways
n_carbon = np.array([6, 6, 6, 6, 6, 0])

# Construct the model:
ecoli_reduced_cb = Cybernetic_Functions.Cybernetic_Model('CB model for Ecoli reduced matrix ')
ecoli_reduced_cb.Smz = Smz
ecoli_reduced_cb.x0 = initial_x0
ecoli_reduced_cb.kmax = kmax
ecoli_reduced_cb.K = K
ecoli_reduced_cb.ke = ke
ecoli_reduced_cb.alpha = alpha
ecoli_reduced_cb.beta = beta
ecoli_reduced_cb.n_carbon = n_carbon
ecoli_reduced_cb.sub_index = 0
ecoli_reduced_cb.Biomass_index = 1
CB_model = copy.deepcopy(ecoli_reduced_cb)
CB_model.mets_name = metabObj


# CB_model = copy.deepcopy(ecoli_reduced_cb)
# CB_model.kmax[[1,2,3]] = [11,22,33]
# CB_model.K[[1,2]] = [77,55]

# if num_kmax > 0:
#     index = para_to_fit['kmax']
#     CB_model.kmax[index] = [11,22,33]
# num_K = len(para_to_fit['K'])
# if num_K > 0:
#     index = para_to_fit['K']
#     CB_model.K[index] = [44,55]

def rate_def(self, x):
    # print('def rate')
    Smz = self.Smz
    kmax = self.kmax
    K = self.K
    ke = self.ke
    alpha = self.alpha
    beta = self.beta
    Biomass_index = self.Biomass_index
    sub_index = self.sub_index

    (n_mets, n_path) = Smz.shape

    r_kin_basic = [
        kmax[0] * x[sub_index] / (K[0] + x[sub_index]),
        kmax[1] * x[sub_index] / (K[1] + x[sub_index]),
        kmax[2] * x[sub_index] / (K[2] + x[sub_index]),
        kmax[3] * x[sub_index] / (K[3] + x[sub_index]),
        kmax[4] * x[sub_index] / (K[4] + x[sub_index]),
        kmax[5] * (x[3] ** 2) / ((K[5] ** 2) + (x[3] ** 2)), ]

    rM = r_kin_basic * x[n_mets:]

    rE = ke * r_kin_basic / kmax

    rG = Smz[Biomass_index, :] * rM[:]  # 11: biomass index

    return rM, rE, rG


setattr(Cybernetic_Functions.Cybernetic_Model, 'rate_def', rate_def)


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

    if CB_model.name in ['CB model for Ecoli reduced matrix ', 'CB model for Ecoli iML1515 matrix ']:
        u[-1] = 1.0
        v[-1] = 1.0
    return u, v


setattr(Cybernetic_Functions.Cybernetic_Model, 'cybernetic_var_def', cybernetic_var_def)
# CB_model.cybernetic_var_def = cybernetic_var_def

CB_model.experiment_data_df = experiment_data_df
CB_model.weights_of_mets = 1
sol = Cybernetic_Functions.cb_model_simulate(CB_model, tspan, draw=True)

# %% Branch_work mutile models
# para_to_fit = {'kmax': np.array([[0, np.nan],
#                                  [1, np.nan],
#                                  [4, 4],
#                                  [5, 5]]),
#
#                'K': np.array([[1, 1],
#                               [2, 2],
#                               [4, 3],
#                               [5, 5]])}
#
# CB_models = [CB_model, CB_model]
# experiment_data_dfs = [experiment_data_df, experiment_data_df]
#
# minimum = Cybernetic_Functions.parameters_fitting(CB_models, experiment_data_dfs, para_to_fit, tspan, draw=True,
#                                                   full_output=False,
#                                                   options={'xtol': 0.01, 'ftol': 0.01, 'maxiter': 50,
#                                                            'disp': True})
# CB_models = Cybernetic_Functions.update_paras_func(minimum.x, CB_models, para_to_fit,
#                                                    retern_model_or_paras='model')
# sol1 = Cybernetic_Functions.cb_model_simulate(CB_models[0], tspan, draw=True)
# sol2 = Cybernetic_Functions.cb_model_simulate(CB_models[1], tspan, draw=True)

# %% Branch_work weight_of_method_t_f

para_to_fit = {'kmax': np.array([[0], [1], [4], [5]]), 'K': np.array([[1], [2], [4], [5]])}

CB_models = [CB_model]
experiment_data_dfs = [experiment_data_df]

minimums = []
for i in range(0, 7):
    print(i)
    CB_model.weight_of_method_t_f = i
    minimum = Cybernetic_Functions.parameters_fitting(CB_models, experiment_data_dfs, para_to_fit, tspan, draw=True,
                                                      full_output=False,
                                                      options={'xtol': 1, 'ftol': 0.1, 'maxiter': 50,
                                                               'disp': True})
    minimums.append(minimum)

for minimum in minimums:
    CB_models = Cybernetic_Functions.update_paras_func(minimum.x, [CB_model], para_to_fit,
                                                       retern_model_or_paras='model')
    sol3 = Cybernetic_Functions.cb_model_simulate(CB_models[0], tspan, draw=True)
# result： weight_of_method_t_f ： 2 > 0 > 1 > others
