#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 6/12/20

"""Cybernetic_F_test_01.py
:description : script
:param :
:returns:
:rtype:
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import ComplementaryScripts.Cybernetic_Functions as Cybernetic_Functions

# %% < def input:>
os.chdir('../../ComplementaryData/test')

# time def
tStart = 0.0
tStop = 100
tStep = 1
tspan = np.linspace(tStart, tStop, int(((tStop - tStart) / tStep) + 1))

# matrix SmZ: reactions x pathways; and Smz metabolites x pathways , S metabolites x reactions : Smz = Sm @ Z
print('\n---------- Loading FBA modes ... ---------- ')

Smz1 = np.loadtxt('yest_smz1.txt')[0:6, :]
Smz2 = np.loadtxt('yest_smz2.txt')[0:6, :]
Smz3 = np.loadtxt('yest_smz3.txt')[0:6, :]

# metabolites and pathways number
(n_mets1, n_path1) = Smz1.shape
(n_mets2, n_path2) = Smz2.shape
(n_mets3, n_path3) = Smz3.shape

# experiment data ti get initial mets
print('\n---------- Loading Experiment Data ... ---------- ')

# experiment_data_df = pd.read_csv('Case1_ecoli_reduced/ecoli_reduce_experiment_data.txt', delimiter='\t',
#                                  header=0)


experiment_data_df_1 = pd.read_excel(r'experiemntdata_to_fit.xlsx', sheet_name='model1')

experiment_data_df_2 = pd.read_excel(r'experiemntdata_to_fit.xlsx', sheet_name='model2')

experiment_data_df_3 = pd.read_excel(r'experiemntdata_to_fit.xlsx', sheet_name='model3')

experiment_data_df_12 = pd.read_excel(r'experiemntdata_to_fit.xlsx', sheet_name='model1_and_model2')

experiment_data_df_23 = pd.read_excel(r'experiemntdata_to_fit.xlsx', sheet_name='model2_and_model3')

# experiment_data_df_13 = pd.read_excel(r'experiemntdata_to_fit.xlsx', sheet_name='model1_and_model3')
#
# experiment_data_df_123 = pd.read_excel(r'experiemntdata_to_fit.xlsx', sheet_name='model1_and_model2_and_model3')


mets_name = list(experiment_data_df_1.columns[1:])  # ['GLC', 'BIOM', 'ETHx', 'XYL', 'MAN', 'GAL' ]

# initial metabolites at tome 0, t0: initial_mets.shape = (n_mets,)
molar_mass = np.array([44.03 * 26.04, 46, 180, 150, 180, 180])  # unit conversion : 0.01g/l = 1000*0.01/900 mM(mmol/l)

experiment_data_df_1.iloc[:, 1:] = experiment_data_df_1.values[:, 1:] * 1000 / molar_mass

experiment_data_df_2.iloc[:, 1:] = experiment_data_df_2.values[:, 1:] * 1000 / molar_mass

experiment_data_df_3.iloc[:, 1:] = experiment_data_df_3.values[:, 1:] * 1000 / molar_mass

experiment_data_df_12.iloc[:, 1:] = experiment_data_df_12.values[:, 1:] * 1000 / molar_mass

experiment_data_df_23.iloc[:, 1:] = experiment_data_df_23.values[:, 1:] * 1000 / molar_mass

initial_mets1 = experiment_data_df_1[mets_name].values[0, :]

initial_mets2 = experiment_data_df_2[mets_name].values[0, :]

initial_mets3 = experiment_data_df_3[mets_name].values[0, :]

initial_mets12 = experiment_data_df_12[mets_name].values[0, :]

initial_mets23 = experiment_data_df_23[mets_name].values[0, :]

# initial_mets13 = experiment_data_df_13[mets_name].values[0, :] * 1000 / molar_mass
#
# initial_mets123 = experiment_data_df_123[mets_name].values[0, :] * 1000 / molar_mass


# Metabolites rate Parameters kmax , Ki : dm/dt =  Smz @ V @ rM(kmax,K) * c
# TODO the equations not the same , should defined by user

# kmax : n_path

kmax1 = np.array([
    27.44,
    11.72,
    39.35,
    22.27,
    14.65,
    1.208,
    0.78,
    1.033,
])

kmax2 = np.array([
    13.94,
    94.39,
    27.64,
    106.9,
    20.68,
    18.73,
    199.9,
    15.65,
    17.41,
    89.13,
])

kmax3 = np.array([
    20.35,
    105.4,
    14.28,
    16.83,
    22.07,
    275.6,
    8.094,
    34.29,
])

# K : n_path
K1 = np.array([
    3.4,
    3.4,
    3.4,
    71,
    71,
    3.2,
    3.2,
    3.2,

    200
])

K2 = np.array([
    4.2,
    4.2,
    300,
    300,
    71,
    71,
    71,
    50,
    50,
    50,

    200
])

K3 = np.array([
    4.2,
    4.2,
    300,
    300,
    71,
    71,
    50,
    50,

    200
])


# Construct the model:
class Cybernetic_Model_basic(Cybernetic_Functions.Cybernetic_Model):
    pass


cb_model_1 = Cybernetic_Model_basic('CB model 1')
cb_model_1.Smz = Smz1
cb_model_1.mets_name = mets_name
cb_model_1.initial_mets = initial_mets1
cb_model_1.initial_enzymes = np.array([0.09] * n_path1)
cb_model_1.kmax = kmax1
cb_model_1.K = K1
cb_model_1.ke = np.array([1] * n_path1)
cb_model_1.alpha = np.array([0.1] * n_path1)
cb_model_1.beta = np.array([0.2] * n_path1)
cb_model_1.n_carbon = np.array([6] * n_path1)
cb_model_1.sub_index = np.array([2] * n_path1)
cb_model_1.Biomass_index = np.array([0] * n_path1)  # np.array([0] * n_path1)
cb_model_1.mets_name = mets_name
cb_model_1.experiment_data_df = experiment_data_df_1

cb_model_2 = Cybernetic_Model_basic('CB model 2')
cb_model_2.Smz = Smz2
cb_model_2.mets_name = mets_name
cb_model_2.initial_mets = initial_mets2
cb_model_2.initial_enzymes = np.array([0.09] * n_path2)
cb_model_2.kmax = kmax2
cb_model_2.K = K2
cb_model_2.ke = np.array([1] * n_path2)
cb_model_2.alpha = np.array([0.1] * n_path2)
cb_model_2.beta = np.array([0.2] * n_path2)
cb_model_2.n_carbon = np.array([6, 6, 5, 5, 6, 6, 6, 6, 6, 6])
cb_model_2.sub_index = np.array([2] * n_path2)
cb_model_2.Biomass_index = np.array([0] * n_path2)  # np.array([0] * n_path2)
cb_model_2.mets_name = mets_name
cb_model_2.experiment_data_df = experiment_data_df_2

cb_model_3 = Cybernetic_Model_basic('CB model 3')
cb_model_3.Smz = Smz3
cb_model_3.mets_name = mets_name
cb_model_3.initial_mets = initial_mets3
cb_model_3.initial_enzymes = np.array([0.09] * n_path3)
cb_model_3.kmax = kmax3
cb_model_3.K = K3
cb_model_3.ke = np.array([1] * n_path3)
cb_model_3.alpha = np.array([0.1] * n_path3)
cb_model_3.beta = np.array([0.2] * n_path3)
cb_model_3.n_carbon = np.array([6, 6, 5, 5, 6, 6, 6, 6])
cb_model_3.sub_index = np.array([2] * n_path3)
cb_model_3.Biomass_index = np.array([0] * n_path3)  # np.array([0] * n_path2)
cb_model_3.mets_name = mets_name
cb_model_3.experiment_data_df = experiment_data_df_3

cb_model_12 = cb_model_1.expend_model(cb_model_2, name='CB model 12')
cb_model_12.initial_mets = initial_mets12
cb_model_12.experiment_data_df = experiment_data_df_12

cb_model_23 = cb_model_2.expend_model(cb_model_3, name='CB model 23')
cb_model_23.initial_mets = initial_mets23
cb_model_23.experiment_data_df = experiment_data_df_23

cb_model_13 = cb_model_1.expend_model(cb_model_3, name='CB model 13')
# cb_model_13.initial_mets = initial_mets13

cb_model_123 = cb_model_12.expend_model(cb_model_23, name='CB model 123')


# cb_model_123.initial_mets = initial_mets123


# CB_model = copy.deepcopy(ecoli_reduced_cb)

def rate_def(self, x):
    # print('def rate')
    name = self.name
    Smz = self.Smz
    kmax = self.kmax
    K = self.K
    ke = self.ke
    Biomass_index = self.Biomass_index
    sub_index = self.sub_index

    (n_mets, n_path) = Smz.shape
    c_glc = x[2]
    c_xyl = x[3]
    c_man = x[4]
    c_gal = x[5]
    c_eth = x[1]
    KinhIS = K[-1]

    if name == 'CB model 1':

        sub_arrary = np.array([c_glc, c_glc, c_glc,
                               c_man, c_man,
                               c_gal, c_gal, c_gal])

        # r_kin_basic = [
        #     kmax[0] * c_glc / (K[0] + c_glc) / (1 + c_eth / KinhIS),
        #     kmax[1] * c_glc / (K[1] + c_glc) / (1 + c_eth / KinhIS),
        #     kmax[2] * c_glc / (K[2] + c_glc) / (1 + c_eth / KinhIS),
        #
        #     kmax[3] * c_man / (K[3] + c_man) / (1 + c_eth / KinhIS),
        #     kmax[4] * c_man / (K[4] + c_man) / (1 + c_eth / KinhIS),
        #
        #     kmax[5] * c_gal / (K[5] + c_gal) / (1 + c_eth / KinhIS),
        #     kmax[6] * c_gal / (K[6] + c_gal) / (1 + c_eth / KinhIS),
        #     kmax[7] * c_gal / (K[7] + c_gal) / (1 + c_eth / KinhIS),
        # ]

    elif name == 'CB model 2':

        sub_arrary = np.array([c_glc, c_glc,
                               c_xyl, c_xyl,
                               c_man, c_man, c_man,
                               c_gal, c_gal, c_gal])

    elif name == 'CB model 3':

        sub_arrary = np.array([c_glc, c_glc,
                               c_xyl, c_xyl,
                               c_man, c_man,
                               c_gal, c_gal, ])

    elif name == 'CB model 12':
        sub_arrary = np.array([c_glc, c_glc, c_glc,
                               c_man, c_man,
                               c_gal, c_gal, c_gal,

                               c_glc, c_glc,
                               c_xyl, c_xyl,
                               c_man, c_man, c_man,
                               c_gal, c_gal, c_gal
                               ])

    elif name == 'CB model 23':
        sub_arrary = np.array([
            c_glc, c_glc,
            c_xyl, c_xyl,
            c_man, c_man, c_man,
            c_gal, c_gal, c_gal,

            c_glc, c_glc,
            c_xyl, c_xyl,
            c_man, c_man,
            c_gal, c_gal,
        ])


    elif name == 'CB model 13':
        sub_arrary = np.array([c_glc, c_glc, c_glc,
                               c_man, c_man,
                               c_gal, c_gal, c_gal,

                               c_glc, c_glc,
                               c_xyl, c_xyl,
                               c_man, c_man,
                               c_gal, c_gal,
                               ])

    elif name == 'CB model 123':
        sub_arrary = np.array([c_glc, c_glc, c_glc,
                               c_man, c_man,
                               c_gal, c_gal, c_gal,

                               c_glc, c_glc,
                               c_xyl, c_xyl,
                               c_man, c_man, c_man,
                               c_gal, c_gal, c_gal,

                               c_glc, c_glc,
                               c_xyl, c_xyl,
                               c_man, c_man,
                               c_gal, c_gal,
                               ])
    else:
        print('Error name')

    r_kin_basic = kmax * sub_arrary / (K[0:-1] + sub_arrary) / (1 + c_eth / KinhIS)

    rM = r_kin_basic * x[n_mets:]

    rE = ke * r_kin_basic / kmax

    if type(Biomass_index) != int and len(Biomass_index) == n_path:
        rG = rM[:] * np.array([Smz[Biomass_index[i], i] for i in np.arange(0, n_path)])
    else:
        rG = rM[:] * Smz[Biomass_index, :]  # biomass index

    return rM, rE, rG


setattr(Cybernetic_Model_basic, 'rate_def', rate_def)

# def cybernetic_var_def(self,rM,i_index):
#     # print('def cybernetic_var')
#     (n_mets, n_path) = self.Smz.shape
#     # cybernetic_var = abs(Smz[sub_index, :] * rM * n_carbon)
#     cybernetic_var = rM * self.n_carbon
#     # cybernetic_var[cybernetic_var==0] = 1
#     cybernetic_var[cybernetic_var < 0] = 0
#     if sum(cybernetic_var) > 0:
#         u = cybernetic_var / sum(cybernetic_var)
#         v = cybernetic_var / np.max(abs(cybernetic_var))
#     else:
#         u = v = np.zeros(n_path)
#
#     if self.name in ['CB model for Ecoli reduced matrix ', 'CB model for Ecoli iML1515 matrix ']:
#         u[-1] = 1.0
#         v[-1] = 1.0
#     return u, v
# setattr(Cybernetic_Model_basic, 'cybernetic_var_def', cybernetic_var_def)


# CB_model.cybernetic_var_def = cybernetic_var_def

sol1 = Cybernetic_Functions.cb_model_simulate(cb_model_1, tspan, draw=True)
sol2 = Cybernetic_Functions.cb_model_simulate(cb_model_2, tspan, draw=True)
sol3 = Cybernetic_Functions.cb_model_simulate(cb_model_3, tspan, draw=True)

sol12 = Cybernetic_Functions.cb_model_simulate(cb_model_12, tspan, draw=True)
sol23 = Cybernetic_Functions.cb_model_simulate(cb_model_23, tspan, draw=True)

# %%

para_to_fit = {'kmax': np.array([[0, np.nan, 0],
                                 [1, np.nan, 1],
                                 [2, np.nan, 2],
                                 [3, np.nan, 3],
                                 [4, np.nan, 4],
                                 [5, np.nan, 5],
                                 [6, np.nan, 6],
                                 [7, np.nan, 7],

                                 [np.nan, 0, 0 + 7],
                                 [np.nan, 1, 1 + 7],
                                 [np.nan, 2, 2 + 7],
                                 [np.nan, 3, 3 + 7],
                                 [np.nan, 4, 4 + 7],
                                 [np.nan, 5, 5 + 7],
                                 [np.nan, 6, 6 + 7],
                                 [np.nan, 7, 7 + 7],
                                 ]),

               'K': np.array([
                   [0, np.nan, 0],
                   [1, np.nan, 1],
                   [2, np.nan, 2],
                   [3, np.nan, 3],
                   [4, np.nan, 4],
                   [5, np.nan, 5],
                   [6, np.nan, 6],
                   [7, np.nan, 7],

                   [np.nan, 0, 0 + 7],
                   [np.nan, 1, 1 + 7],
                   [np.nan, 2, 2 + 7],
                   [np.nan, 3, 3 + 7],
                   [np.nan, 4, 4 + 7],
                   [np.nan, 5, 5 + 7],
                   [np.nan, 6, 6 + 7],
                   [np.nan, 7, 7 + 7], ])}

cb_models = [cb_model_1, cb_model_2, cb_model_12]
experiment_data_dfs = [experiment_data_df_1, experiment_data_df_2, experiment_data_df_12]

# minimum = Cybernetic_Functions.parameters_fitting(cb_models, experiment_data_dfs,
#                                                   para_to_fit, tspan, draw=True,
#                                                   full_output=False,
#                                                   options={'xtol': 0.01, 'ftol': 0.01, 'maxiter': 50, 'disp': True})
#
# CB_models = Cybernetic_Functions.update_paras_func(minimum.x, cb_models, para_to_fit,
#                                                    retern_model_or_paras='model')
# sol = Cybernetic_Functions.cb_model_simulate(CB_models[0], tspan, draw=True)
# sol = Cybernetic_Functions.cb_model_simulate(CB_models[1], tspan, draw=True)

# %%

sols = [sol1, sol2, sol3, sol12]
experiment_data_df_s = [experiment_data_df_1, experiment_data_df_2, experiment_data_df_3,  # ,experiment_data_df_23
                        experiment_data_df_12, ]

results = []
for sol_i in sols:
    result_i = sol_i[:, 0:6] * molar_mass / 1000
    results.append(result_i)

for experiment_data_df_i in experiment_data_df_s:
    experiment_data_df_i.iloc[:, 1:] = experiment_data_df_i.values[:, 1:] * molar_mass / 1000

# %%
fig = plt.figure(figsize=(5.5, 10))
axs = [fig.add_subplot(len(results), 1, i + 1) for i in range(len(results))]
results
# fig,axs = plt.subplots(len(results), 1,)

for i in range(len(results)):
    axs[i].plot(tspan, results[i], )
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    color_i = 0
    for met_i in mets_name:
        axs[i].plot(experiment_data_df_s[i]['t'], experiment_data_df_s[i][met_i], '*', color=colors[color_i])
        color_i += 1
    # axs[i].set_xlim(-1,90)

axs[0].legend(mets_name)

plt.show()

# %%
# CB_model['metas_names'] = ['glc', 'succ', 'for', 'lac', 'ac', 'etoh', 'biomass', ]
# para_to_fit = {'kmax': [0, 1, 2, 3, 4, 5], 'K': [0, 1, 2, 3, 4, 5]}
# # minimum = Cybernetic_Functions.parameters_fitting(CB_model, experiment_data_df, para_to_fit, tspan, draw=True)
#
# # matplotlib.use("TkAgg")
# # CB_model = Cybernetic_Functions.update_paras_func(minimum[0], CB_model, para_to_fit, retern_model_or_paras='model')
# sol = Cybernetic_Functions.cb_model_simulate(CB_model, tspan, draw=True)
