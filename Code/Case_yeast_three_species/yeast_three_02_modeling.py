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
import pickle

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import Code.Cybernetic_Functions as Cybernetic_Functions

# %% < def input:>
os.chdir('../../Data/Case_yeast_three_species')

# time def
tStart = 0.0
tStop = 100
tStep = 1
tspan = np.linspace(tStart, tStop, int(((tStop - tStart) / tStep) + 1))

# matrix SmZ: reactions x pathways; and Smz metabolites x pathways , S metabolites x reactions : Smz = Sm @ Z
print('\n---------- Loading FBA modes ... ---------- ')
Smz1_ = np.loadtxt('yest_smz1.txt')[0:6, :]
Smz2_ = np.loadtxt('yest_smz2.txt')[0:6, :]
Smz3_ = np.loadtxt('yest_smz3.txt')[0:6, :]

model_smz = pickle.load(open('model_smz', 'rb'))

Smz1 = model_smz[0]
Smz2 = model_smz[1]
Smz3 = model_smz[2]

Smz1[[0, 1, 2, 3, 4, 5], :] = Smz1[[4, 5, 0, 1, 2, 3], :]
Smz2[[0, 1, 2, 3, 4, 5], :] = Smz2[[4, 5, 0, 1, 2, 3], :]
Smz3[[0, 1, 2, 3, 4, 5], :] = Smz3[[4, 5, 0, 1, 2, 3], :]

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


mets_name = list(experiment_data_df_1.columns[1:])  # ['BIOM', 'ETHx', 'GLC', 'XYL', 'MAN', 'GAL']

# initial metabolites at tome 0, t0: initial_mets.shape = (n_mets,)
molar_mass = np.array([1040, 46, 180, 150, 180, 180])  # unit conversion : 0.01g/l = 1000*0.01/900 mM(mmol/l)

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

kmax1 = np.array(
    [110.91280957986572, 1.1854850622932065, 13.682176672912881, 2.6167266723859317, 3.932561053042331,
     11.018377153500872, 0.39686934689377384, 1.4919147628099374, 0.606451588748099, ])
kmax2 = np.array([23.23092935425619, 0.17319751654601012, 1.2866188498323905, 20.113238231799215, 11.191576816882012,
                  7.029611542618358, 32.885878446402515, 15.059416224769848, 4.4997263093736635, 19.63780820651631, ])
kmax3 = np.array([17.869681531987467, 12.630842188798638, 23.948765514322844, 2.333676728126207, 5.6808897607972675,
                  0.3965147830533877, 0.5685312862222127, 42.021225141220526, 0.5747711617822211, 8.730434600267824])
# np.around(Smz1,4)
# K : n_path
K1 = np.array([
    3.2,
    3.2,
    3.2,
    71,
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
    4.2,
    300,
    300,
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
    300,
    71,
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
cb_model_1.experiment_data_df = experiment_data_df_1.copy()

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
cb_model_2.n_carbon = np.array([6, 6, 6, 5, 5, 6, 6, 6, 6, 6])
cb_model_2.sub_index = np.array([2] * n_path2)
cb_model_2.Biomass_index = np.array([0] * n_path2)  # np.array([0] * n_path2)
cb_model_2.mets_name = mets_name
cb_model_2.experiment_data_df = experiment_data_df_2.copy()

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
cb_model_3.n_carbon = np.array([6, 6, 5, 5, 5, 6, 6, 6, 6, 6])
cb_model_3.sub_index = np.array([2] * n_path3)
cb_model_3.Biomass_index = np.array([0] * n_path3)  # np.array([0] * n_path2)
cb_model_3.mets_name = mets_name
cb_model_3.experiment_data_df = experiment_data_df_3.copy()

cb_model_12 = cb_model_1.expend_model(cb_model_2, name='CB model 12')
cb_model_12.initial_mets = initial_mets12
cb_model_12.experiment_data_df = experiment_data_df_12.copy()
cb_model_12.K = np.concatenate((cb_model_1.K[0:-1], cb_model_2.K), )
cb_model_23 = cb_model_2.expend_model(cb_model_3, name='CB model 23')
cb_model_23.initial_mets = initial_mets23
cb_model_23.experiment_data_df = experiment_data_df_23.copy()
cb_model_23.K = np.concatenate((cb_model_2.K[0:-1], cb_model_2.K), )

cb_model_13 = cb_model_1.expend_model(cb_model_3, name='CB model 13')
cb_model_13.K = np.concatenate((cb_model_1.K[0:-1], cb_model_3.K), )

# cb_model_13.initial_mets = initial_mets13

cb_model_123 = cb_model_12.expend_model(cb_model_3, name='CB model 123')
cb_model_123.K = np.concatenate((cb_model_12.K[0:-1], cb_model_3.K), )


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

    sub_arrary_1 = np.array([c_glc, c_glc, c_glc,
                             c_man, c_man, c_man,
                             c_gal, c_gal, c_gal])

    sub_arrary_2 = np.array([c_glc, c_glc, c_glc,
                             c_xyl, c_xyl,
                             c_man, c_man,
                             c_gal, c_gal, c_gal])

    sub_arrary_3 = np.array([c_glc, c_glc,
                             c_xyl, c_xyl, c_xyl,
                             c_man, c_man, c_man,
                             c_gal, c_gal, ])

    if name == 'CB model 1':

        sub_arrary = sub_arrary_1

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
        sub_arrary = sub_arrary_2

    elif name == 'CB model 3':
        sub_arrary = sub_arrary_3


    elif name == 'CB model 12':
        sub_arrary = np.append(sub_arrary_1, sub_arrary_2)

    elif name == 'CB model 23':
        sub_arrary = np.append(sub_arrary_2, sub_arrary_3)

    elif name == 'CB model 13':
        sub_arrary = np.append(sub_arrary_1, sub_arrary_3)


    elif name == 'CB model 123':
        sub_arrary = np.append(sub_arrary_1, sub_arrary_2)
        sub_arrary = np.append(sub_arrary, sub_arrary_3)

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

# def cybernetic_var_def(self, rM, i_index):
#     name = self.name
#     # print('def cybernetic_var')
#     (n_mets, n_path) = self.Smz.shape
#     # cybernetic_var = abs(Smz[sub_index, :] * rM * n_carbon)
#     cybernetic_var = rM * self.n_carbon
#     # cybernetic_var[cybernetic_var==0] = 1
#     cybernetic_var[cybernetic_var < 0] = 0
#     u = v = np.ones((1, n_path))
#
#     if name in ['CB model 12', 'CB model 13']:
#         fam = [9, 10]
#     elif name in ['CB model 23']:
#         fam = [10, 10]
#     elif name in ['CB model 123']:
#         fam = [9, 10, 10]
#     else:
#         fam = [-1]
#     for index_i in range(0, len(fam)):
#         if index_i == 0:
#             index_strat = 0
#         else:
#             index_strat = fam[index_i - 1]
#         index_end = fam[index_i]
#
#         cybernetic_var_i = cybernetic_var[index_strat:index_end]
#         # u_i = u[index_strat:index_end]
#         # v_i = v[index_strat:index_end]
#
#         if sum(cybernetic_var_i) > 0:
#             u_i = cybernetic_var_i / sum(cybernetic_var_i)
#             v_i = cybernetic_var_i / np.max(abs(cybernetic_var_i))
#         else:
#             u_i = v_i = np.zeros(len(cybernetic_var_i))
#
#         u[index_strat:index_end] = u_i
#         v[index_strat:index_end] = v_i
#
#     return u, v
#
#
# setattr(Cybernetic_Model_basic, 'cybernetic_var_def', cybernetic_var_def)

# CB_model.cybernetic_var_def = cybernetic_var_def
# %% initial model results:
sol1 = Cybernetic_Functions.cb_model_simulate(cb_model_1, tspan, draw=True)
sol2 = Cybernetic_Functions.cb_model_simulate(cb_model_2, tspan, draw=True)
sol3 = Cybernetic_Functions.cb_model_simulate(cb_model_3, tspan, draw=True)

sol12 = Cybernetic_Functions.cb_model_simulate(cb_model_12, tspan, draw=True)
sol23 = Cybernetic_Functions.cb_model_simulate(cb_model_23, tspan, draw=True)
sol13 = Cybernetic_Functions.cb_model_simulate(cb_model_13, tspan, draw=True)
sol123 = Cybernetic_Functions.cb_model_simulate(cb_model_123, tspan, draw=True)

file = open('sol_all', 'wb')
pickle.dump([sol1, sol2, sol3, sol12, sol23, sol13, sol123], file)
file.close()

weights_of_mets = molar_mass / 1000

# %% model1

if False:
    cb_model_1.weight_of_method_t_f = 3
    cb_model_1.weights_of_mets = np.array([5, 0.1, 2, 0, 2, 2])  # /1000
    cb_models = [cb_model_1, ]
    experiment_data_dfs = [experiment_data_df_1, ]  # [5],[6],[7],[8]]
    para_to_fit = {'kmax': np.array([[0], [1], [2], [3], [4], [5], [6], [7], [8]]),
                   'K': np.array([])}

    minimum_1 = Cybernetic_Functions.parameters_fitting(cb_models, experiment_data_dfs,
                                                        para_to_fit, tspan, draw=True,
                                                        full_output=False,
                                                        options={'xtol': 1, 'ftol': 0.01, 'maxiter': 100,
                                                                 'disp': True})

    CB_models = Cybernetic_Functions.update_paras_func(minimum_1.x, cb_models, para_to_fit,
                                                       retern_model_or_paras='model')
    sol1 = Cybernetic_Functions.cb_model_simulate(CB_models[0], tspan, draw=True)

    cb_model_2.weight_of_method_t_f = 3
    cb_models = [cb_model_2, ]
    experiment_data_dfs = [experiment_data_df_2, ]
    para_to_fit = {'kmax': np.array([[0], [1], [2], [3], [4], [5], [6], [7], [8], [9]]),
                   'K': np.array([])}

    minimum_2 = Cybernetic_Functions.parameters_fitting(cb_models, experiment_data_dfs,
                                                        para_to_fit, tspan, draw=True,
                                                        full_output=False,
                                                        options={'xtol': 0.5, 'ftol': 0.5, 'maxiter': 50,
                                                                 'disp': True})

    CB_models = Cybernetic_Functions.update_paras_func(minimum_2.x, cb_models, para_to_fit,
                                                       retern_model_or_paras='model')
    sol2 = Cybernetic_Functions.cb_model_simulate(CB_models[0], tspan, draw=True)

    cb_model_3.weight_of_method_t_f = 3
    cb_models = [cb_model_3, ]
    experiment_data_dfs = [experiment_data_df_3, ]
    para_to_fit = {'kmax': np.array([[0], [1], [2], [3], [4], [5], [6], [7], [8], [9]]),
                   'K': np.array([])}

    minimum_3 = Cybernetic_Functions.parameters_fitting(cb_models, experiment_data_dfs,
                                                        para_to_fit, tspan, draw=True,
                                                        full_output=False,
                                                        options={'xtol': 0.5, 'ftol': 0.1, 'maxiter': 50,
                                                                 'disp': True})

    CB_models = Cybernetic_Functions.update_paras_func(minimum_3.x, cb_models, para_to_fit,
                                                       retern_model_or_paras='model')
    sol3 = Cybernetic_Functions.cb_model_simulate(CB_models[0], tspan, draw=True)

# %%  model1 and model2

if False:
    weight_of_method_t_f = 3
    # cb_model_1.weights_of_mets = np.array([5 , 0.1, 2 , 2 , 2 , 2 ])

    cb_models = [cb_model_1, cb_model_2, cb_model_12]
    for model_i in cb_models:
        model_i.weight_of_method_t_f = weight_of_method_t_f
        model_i.weights_of_mets = np.array([100, 0.1, 1, 1, 2, 2])

    experiment_data_dfs = [experiment_data_df_1, experiment_data_df_2, experiment_data_df_12]
    para_to_fit = {'kmax': np.array([
        [0, np.nan, 0],
        [1, np.nan, 1],
        [2, np.nan, 2],
        [3, np.nan, 3],
        [4, np.nan, 4],
        [5, np.nan, 5],
        [6, np.nan, 6],
        [7, np.nan, 7],
        [8, np.nan, 8],

        [np.nan, 0, 0 + n_path1],
        [np.nan, 1, 1 + n_path1],
        [np.nan, 2, 2 + n_path1],
        [np.nan, 3, 3 + n_path1],
        [np.nan, 4, 4 + n_path1],
        [np.nan, 5, 5 + n_path1],
        [np.nan, 6, 6 + n_path1],
        [np.nan, 7, 7 + n_path1],
        [np.nan, 8, 8 + n_path1],
        [np.nan, 9, 9 + n_path1],
    ]),

        'K': np.array([])}

    minimum = Cybernetic_Functions.parameters_fitting(cb_models, experiment_data_dfs,
                                                      para_to_fit, tspan, draw=True,
                                                      full_output=False,
                                                      options={'xtol': 0.5, 'ftol': 0.01, 'maxiter': 50,
                                                               'disp': True})

    CB_models = Cybernetic_Functions.update_paras_func(minimum.x, cb_models, para_to_fit,
                                                       retern_model_or_paras='model')
    sol2 = Cybernetic_Functions.cb_model_simulate(CB_models[0], tspan, draw=True)
    sol3 = Cybernetic_Functions.cb_model_simulate(CB_models[1], tspan, draw=True)
    sol23 = Cybernetic_Functions.cb_model_simulate(CB_models[2], tspan, draw=True)
    print(stop)

# %% model2 and model3

if False:
    weight_of_method_t_f = 3
    cb_models = [cb_model_2, cb_model_3, cb_model_23]
    for model_i in cb_models:
        model_i.weight_of_method_t_f = weight_of_method_t_f
        model_i.weights_of_mets = np.array([100, 0.1, 1, 1, 2, 2])

    experiment_data_dfs = [experiment_data_df_2, experiment_data_df_3, experiment_data_df_23]
    para_to_fit = {'kmax': np.array([
        [0, np.nan, 0],
        [1, np.nan, 1],
        [2, np.nan, 2],
        [3, np.nan, 3],
        [4, np.nan, 4],
        [5, np.nan, 5],
        [6, np.nan, 6],
        [7, np.nan, 7],
        [8, np.nan, 8],
        [9, np.nan, 9],

        [np.nan, 0, 0 + n_path2],
        [np.nan, 1, 1 + n_path2],
        [np.nan, 2, 2 + n_path2],
        [np.nan, 3, 3 + n_path2],
        [np.nan, 4, 4 + n_path2],
        [np.nan, 5, 5 + n_path2],
        [np.nan, 6, 6 + n_path2],
        [np.nan, 7, 7 + n_path2],
        [np.nan, 8, 8 + n_path2],
        [np.nan, 9, 9 + n_path2],
    ]),

        'K': np.array([])}

    minimum = Cybernetic_Functions.parameters_fitting(cb_models, experiment_data_dfs,
                                                      para_to_fit, tspan, draw=True,
                                                      full_output=False,
                                                      options={'xtol': 0.5, 'ftol': 0.01, 'maxiter': 50, 'disp': True})

    CB_models = Cybernetic_Functions.update_paras_func(minimum.x, cb_models, para_to_fit,
                                                       retern_model_or_paras='model')
    sol2 = Cybernetic_Functions.cb_model_simulate(CB_models[0], tspan, draw=True)
    sol3 = Cybernetic_Functions.cb_model_simulate(CB_models[1], tspan, draw=True)
    sol23 = Cybernetic_Functions.cb_model_simulate(CB_models[2], tspan, draw=True)

# %% model1 and model2 and model 3

if False:
    weight_of_method_t_f = 3
    weights_of_mets = molar_mass / 1000
    cb_models = [cb_model_1, cb_model_2, cb_model_3, cb_model_12, cb_model_23]

    for model_i in cb_models:
        model_i.weight_of_method_t_f = weight_of_method_t_f
        model_i.weights_of_mets = np.array([100, 0.1, 1, 1, 2, 2])
    experiment_data_dfs = [experiment_data_df_1, experiment_data_df_2, experiment_data_df_3,
                           experiment_data_df_12, experiment_data_df_23]

    para_to_fit = {'kmax': np.array([[0, np.nan, np.nan, 0, np.nan],
                                     [1, np.nan, np.nan, 1, np.nan],
                                     [2, np.nan, np.nan, 2, np.nan],
                                     [3, np.nan, np.nan, 3, np.nan],
                                     [4, np.nan, np.nan, 4, np.nan],
                                     [5, np.nan, np.nan, 5, np.nan],
                                     [6, np.nan, np.nan, 6, np.nan],
                                     [7, np.nan, np.nan, 7, np.nan],
                                     [8, np.nan, np.nan, 8, np.nan],

                                     [np.nan, 0, np.nan, 0 + n_path1, 0],
                                     [np.nan, 1, np.nan, 1 + n_path1, 1],
                                     [np.nan, 2, np.nan, 2 + n_path1, 2],
                                     [np.nan, 3, np.nan, 3 + n_path1, 3],
                                     [np.nan, 4, np.nan, 4 + n_path1, 4],
                                     [np.nan, 5, np.nan, 5 + n_path1, 5],
                                     [np.nan, 6, np.nan, 6 + n_path1, 6],
                                     [np.nan, 7, np.nan, 7 + n_path1, 7],
                                     [np.nan, 8, np.nan, 8 + n_path1, 8],
                                     [np.nan, 9, np.nan, 9 + n_path1, 9],

                                     [np.nan, np.nan, 0, np.nan, 0 + n_path2, ],
                                     [np.nan, np.nan, 1, np.nan, 1 + n_path2, ],
                                     [np.nan, np.nan, 2, np.nan, 2 + n_path2, ],
                                     [np.nan, np.nan, 3, np.nan, 3 + n_path2, ],
                                     [np.nan, np.nan, 4, np.nan, 4 + n_path2, ],
                                     [np.nan, np.nan, 5, np.nan, 5 + n_path2, ],
                                     [np.nan, np.nan, 6, np.nan, 6 + n_path2, ],
                                     [np.nan, np.nan, 7, np.nan, 7 + n_path2, ],
                                     [np.nan, np.nan, 8, np.nan, 8 + n_path2, ],
                                     [np.nan, np.nan, 9, np.nan, 9 + n_path2, ],
                                     ]),
                   'K': []}

    minimum = Cybernetic_Functions.parameters_fitting(cb_models, experiment_data_dfs,
                                                      para_to_fit, tspan, draw=True,
                                                      full_output=False,
                                                      options={'xtol': 5, 'ftol': 0.01, 'maxiter': 50,
                                                               'disp': True})

    CB_models = Cybernetic_Functions.update_paras_func(minimum.x, cb_models, para_to_fit,
                                                       retern_model_or_paras='model')
    sol1 = Cybernetic_Functions.cb_model_simulate(CB_models[0], tspan, draw=True)
    sol2 = Cybernetic_Functions.cb_model_simulate(CB_models[1], tspan, draw=True)
    sol3 = Cybernetic_Functions.cb_model_simulate(CB_models[2], tspan, draw=True)

    sol12 = Cybernetic_Functions.cb_model_simulate(CB_models[3], tspan, draw=True)
    sol23 = Cybernetic_Functions.cb_model_simulate(CB_models[4], tspan, draw=True)
# %% draw
sols = [sol1, sol2, sol3, sol12, sol23, sol123]
experiment_data_dfs_copy = [experiment_data_df_1, experiment_data_df_2, experiment_data_df_3,
                            experiment_data_df_12, experiment_data_df_23]
import copy

experiment_data_dfs = copy.deepcopy(experiment_data_dfs_copy)

# experiment_data_dfs = experiment_data_dfs_copy.copy()
# results = []
# for sol_i in sols:
#     result_i = sol_i[:, 0:6] * molar_mass / 1000
#     results.append(result_i)
#
# for experiment_data_df_i in experiment_data_dfs:
#     experiment_data_df_i.iloc[:, 1:] = experiment_data_df_i.values[:, 1:] * molar_mass / 1000
#
# fig = plt.figure(figsize=(5.5, 10))
# axs = [fig.add_subplot(len(results), 1, i + 1) for i in range(len(results))]
# results
# # fig,axs = plt.subplots(len(results), 1,)
#
# for i in range(len(results)):
#     axs[i].plot(tspan, results[i], )
#     prop_cycle = plt.rcParams['axes.prop_cycle']
#     colors = prop_cycle.by_key()['color']
#     color_i = 0
#     for met_i in mets_name:
#         axs[i].plot(experiment_data_dfs[i]['t'], experiment_data_dfs[i][met_i], '*', color=colors[color_i])
#         color_i += 1
#     # axs[i].set_xlim(-1,90)
#
# axs[0].legend(mets_name)
#
# plt.show()

# %%
# experiment_data_dfs = copy.deepcopy(experiment_data_dfs_copy)
# results = []
# for sol_i in sols:
#     result_i = sol_i[:, 0:6] * molar_mass / 1000
#     results.append(result_i)
#
# for experiment_data_df_i in experiment_data_dfs:
#     experiment_data_df_i.iloc[:, 1:] = experiment_data_df_i.values[:, 1:] * molar_mass / 1000
#
# fig = plt.figure(figsize=(5.5, 10))
# axs = [fig.add_subplot(len(results), 1, i + 1) for i in range(len(results))]
# results
# # fig,axs = plt.subplots(len(results), 1,)
#
# for i in range(len(results)):
#     axs[i].plot(tspan, results[i], )
#     prop_cycle = plt.rcParams['axes.prop_cycle']
#     colors = prop_cycle.by_key()['color']
#     color_i = 0
#     for met_i in mets_name:
#         axs[i].plot(experiment_data_dfs[i]['t'], experiment_data_dfs[i][met_i], '*', color=colors[color_i])
#         color_i += 1
#     # axs[i].set_xlim(-1,90)
#
# axs[0].legend(mets_name)
#
# plt.show()

# %%
experiment_data_dfs = copy.deepcopy(experiment_data_dfs_copy)

# colors = ['blue', 'teal', 'tab:red', 'tab:orange']
results = []
for sol_i in sols:
    result_i = sol_i[:, 0:6] * molar_mass / 1000
    results.append(result_i)

for experiment_data_df_i in experiment_data_dfs:
    experiment_data_df_i.iloc[:, 1:] = experiment_data_df_i.values[:, 1:] * molar_mass / 1000

titles = ['$S. cerevisiae$', '$P. stipitis$', '$K. marxianus$',
          '$S. cerevisiae$ & $P. stipitis$', '$P. stipitis$ & $K. marxianus$',
          '$S. cerevisiae$ & $P. stipitis$ & $K. marxianus$']

for index in range(0, len(results)):
    figsize = (2.7, 2.5)
    if index == 5:
        figsize = (3.42, 2.5)
    fig, axs = plt.subplots(2, 1, figsize=figsize)
    ax_0 = axs[0]
    ax_1 = axs[1]
    # color_list = plt.cm.tab10(np.linspace(0, 1, 12))
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']

    result = results[index]
    color_i = 0
    for met_i in mets_name:
        if met_i == 'BIOM':
            ax_0.plot(tspan, result[:, mets_name.index(met_i)], color=colors[color_i], label=met_i, linewidth=1)
            ax_0.set_yscale('log')
            # ax_0.plt.yscale()
        elif met_i != 'BIOM':
            ax_1.plot(tspan, result[:, mets_name.index(met_i)], color=colors[color_i], label=met_i, linewidth=1)
        color_i += 1

    if index < 5:
        experiment_data_df_i = experiment_data_dfs[index]
        color_i = 0
        for met_i in mets_name:
            if met_i == 'BIOM':
                ax_0.plot(experiment_data_df_i['t'], experiment_data_df_i[met_i], '.', color=colors[color_i], alpha=0.8,
                          linewidth=1, markersize=5)

                ax_0.set_yscale('log')
                # ax_0.plt.yscale()
            elif met_i != 'BIOM':
                ax_1.plot(tspan, result[:, mets_name.index(met_i)], color=colors[color_i], linewidth=1)
                ax_1.plot(experiment_data_df_i['t'], experiment_data_df_i[met_i], '.', color=colors[color_i], alpha=0.8,
                          linewidth=1, markersize=5)
            color_i += 1

    ax_0.axis(ymin=0.01, ymax=100)
    ax_0.set_yticks([0.01, 1, 100])
    ax_0.set_yticklabels([-2, 0, 2])

    ax_0.set_xticks([0, 50, 100])
    ax_0.set_xticklabels([])

    ax_0.set_ylabel('Biomass log (g/L)', fontsize=8, family='Arial', )
    ax_1.set_xlabel('Time (h)', fontsize=8, family='Arial', )
    ax_1.set_ylabel('Concentration (g/L)', fontsize=8, family='Arial', )
    if index == 5:
        legend0 = ax_0.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, prop={'family': 'Arial', 'size': 8}, )
        legend1 = ax_1.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, prop={'family': 'Arial', 'size': 8}, )
    ax_0.set_title(titles[index], fontsize=8, family='Arial', )
    plt.yticks(fontname="Arial", fontsize=8)
    plt.xticks(fontname="Arial", fontsize=8)

    ax_0.tick_params(axis="x", labelsize=8, pad=-1)
    ax_0.tick_params(axis="y", labelsize=8, pad=1),
    # ax_1.tick_params(axis="x", labelsize=8, pad=4)
    ax_1.tick_params(axis="y", labelsize=8, pad=1),
    fig.tight_layout(pad=1, h_pad=1, w_pad=0.5)
    fig.savefig('fig6a_%s.pdf' % titles[index], bbox_inches='tight')
    fig.savefig('fig6a_%s.jpg' % titles[index], bbox_inches='tight', dpi=600.)
    fig.show()

legend0
# CB_model['metas_names'] = ['glc', 'succ', 'for', 'lac', 'ac', 'etoh', 'biomass', ]
# para_to_fit = {'kmax': [0, 1, 2, 3, 4, 5], 'K': [0, 1, 2, 3, 4, 5]}
# # minimum = Cybernetic_Functions.parameters_fitting(CB_model, experiment_data_df, para_to_fit, tspan, draw=True)
#
# # matplotlib.use("TkAgg")
# # CB_model = Cybernetic_Functions.update_paras_func(minimum[0], CB_model, para_to_fit, retern_model_or_paras='model')
# sol = Cybernetic_Functions.cb_model_simulate(CB_model, tspan, draw=True)
