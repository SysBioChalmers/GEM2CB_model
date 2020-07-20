#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 5/5/20

""".py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os

import numpy as np
import pandas as pd

import ComplementaryScripts.Cybernetic_Functions as Cybernetic_Functions

os.chdir('../../ComplementaryData/three_species/')
'''
Mono-cultures
=============

RI: RI_8, RI_14, RI_15, RI_16
FP: FP_4, FP_14, FP_15 
BH: BH_14, BH_15, BH_16

Co-cultures
===========

RI-FP: RI_FP_8, RI_FP_9
RI-BH with initial acetate: RI_BH_4, RI_BH_6
RI-BH without initial acetate: RI_BH_5, RI_BH_7 
FP-BH with initial acetate: FP_BH_1, FP_BH_2
FP-BH without initial acetate: FP_BH_3

Tri-cultures
============

RI-FP-BH: RI_FP_BH_10, RI_FP_BH_11, RI_FP_BH_12, RI_FP_BH_13, RI_FP_BH_14, RI_FP_BH_15
'''
# %% <experiment data> experiment data ti get initial mets
print('\n---------- Loading Experiment Data ... ---------- ')
file_name = 'experiment_data/experiment_data_trimmed.xlsx'
RI_experimet_data = pd.read_excel(file_name, sheet_name='RI_14', index_col=0, header=0, )
FP_experimet_data = pd.read_excel(file_name, sheet_name='FP_4', index_col=0, header=0, )
BH_experimet_data = pd.read_excel(file_name, sheet_name='BH_14', index_col=0, header=0, )

RI_BH_experimet_data = pd.read_excel(file_name, sheet_name='RI_BH_4', index_col=0, header=0, )
RI_FP_experimet_data = pd.read_excel(file_name, sheet_name='RI_FP_8', index_col=0, header=0, )
FP_BH_experimet_data = pd.read_excel(file_name, sheet_name='FP_BH_1', index_col=0, header=0, )

# data_list = [RI_experimet_data,FP_experimet_data,BH_experimet_data,RI_BH_experimet_data,RI_FP_experimet_data,FP_BH_experimet_data]
# for data in data_list:
#     data[['R.i', 'F.p', 'B.h']] = data[['R.i', 'F.p', 'B.h']] / 1e8

# mets_name = list(RI_experimet_data.columns[1:-1])
mets_name = ['fru', 'R.i', 'F.p', 'B.h', 'ac', 'for', 'but']

# %% <Cybernetic model>
print('\n---------- Loading modes ... ---------- ')
Smz_Ri = np.loadtxt('R_int_Smz.csv', delimiter=',')
Smz_Fp = np.loadtxt('F_pra_Smz.csv', delimiter=',')
Smz_Bh = np.loadtxt('B_hyd_Smz.csv', delimiter=',')
Smz_Bh[4, 0:-1] = Smz_Bh[4, 0:-1] / 3

print('\n---------- CB modeling ... ---------- ')
tStart = 0.0  # DefineTime
tStop = 50
tStep = 0.5
tspan = np.linspace(tStart, tStop, int(((tStop - tStart) / tStep) + 1))


# Construct the model:
class Cybernetic_Model_basic(Cybernetic_Functions.Cybernetic_Model):
    pass


# Enzyme Rate Parameters: alpha,beta,ke : de/dt =  alpha + rE(ke) * u - (beta + mu) * e

kmax_Ri = np.array([
    5.695011316129802, 1.2362283799908131, 1.667282755978127, 6.819342968503566,
])

K_Ri = np.array([
    200, 200, 200, 100
])
# [3.4524982181822796, 0.7592036532956223, 0.7353427943754969, 85.993126622221, 78.07181958062853, 38.969324388864244]
kmax_Fp = np.array([
    1.4492041469463008, 1.279073770351991, 1.7301980743535568, 2.209722853192763, 11.825606699093186
])

K_Fp = np.array([
    200, 200, 200, 200, 100
])

kmax_Bh = np.array([
    1.0938534718220758, 1.08030996596831, 1.45922493291296, 165.42609677251724
])

K_Bh = np.array([
    200, 200, 200, 150
])

### model structure
cb_model_Ri = Cybernetic_Model_basic('CB model 1')
cb_model_Ri.Smz = Smz_Ri
(n_mets, n_path) = cb_model_Ri.Smz.shape
cb_model_Ri.experiment_data_df = RI_experimet_data
cb_model_Ri.initial_mets = RI_experimet_data[mets_name].values[0, :]
cb_model_Ri.mets_name = mets_name
cb_model_Ri.kmax = kmax_Ri
cb_model_Ri.K = K_Ri
cb_model_Ri.initial_enzymes = np.array([0.9] * n_path)
cb_model_Ri.ke = np.array([0.620342] * n_path)
cb_model_Ri.alpha = np.array([0.04] * n_path)
cb_model_Ri.beta = np.array([0.05] * n_path)
cb_model_Ri.n_carbon = np.array([6, 6, 6, 2])
cb_model_Ri.sub_index = np.array([0] * n_path)
cb_model_Ri.Biomass_index = np.array([1] * n_path)

cb_model_Fp = Cybernetic_Model_basic('CB model 2')
cb_model_Fp.Smz = Smz_Fp
(n_mets, n_path) = cb_model_Fp.Smz.shape
cb_model_Fp.experiment_data_df = FP_experimet_data
cb_model_Fp.initial_mets = FP_experimet_data[mets_name].values[0, :]
cb_model_Fp.mets_name = mets_name
cb_model_Fp.kmax = kmax_Fp
cb_model_Fp.K = K_Fp
cb_model_Fp.initial_enzymes = np.array([0.9] * n_path)
cb_model_Fp.ke = np.array([0.620342] * n_path)
cb_model_Fp.alpha = np.array([0.04] * n_path)
cb_model_Fp.beta = np.array([0.05] * n_path)
cb_model_Fp.n_carbon = np.array([6, 6, 6, 6, 2])
cb_model_Fp.sub_index = np.array([0] * n_path)
cb_model_Fp.Biomass_index = np.array([2] * n_path)

cb_model_Bh = Cybernetic_Model_basic('CB model 3')
cb_model_Bh.Smz = Smz_Bh
(n_mets, n_path) = cb_model_Bh.Smz.shape
cb_model_Bh.experiment_data_df = BH_experimet_data
cb_model_Bh.initial_mets = BH_experimet_data[mets_name].values[0, :]
cb_model_Bh.mets_name = mets_name
cb_model_Bh.kmax = kmax_Bh
cb_model_Bh.K = K_Bh
cb_model_Bh.initial_enzymes = np.array([0.9] * n_path)
cb_model_Bh.ke = np.array([0.620342] * n_path)
cb_model_Bh.alpha = np.array([0.04] * n_path)
cb_model_Bh.beta = np.array([0.05] * n_path)
cb_model_Bh.n_carbon = np.array([6, 6, 6, 2])
cb_model_Bh.sub_index = np.array([0] * n_path)
cb_model_Bh.Biomass_index = np.array([3] * n_path)

cb_model_1 = cb_model_Ri
cb_model_2 = cb_model_Fp
cb_model_3 = cb_model_Bh

cb_model_12 = cb_model_Ri.expend_model(cb_model_Fp, name='CB model 12')
# cb_model_12.resources_group = [4, 5]
cb_model_12.experiment_data_df = RI_FP_experimet_data.copy()
cb_model_12.initial_mets = cb_model_12.experiment_data_df[mets_name].values[0, :]

cb_model_23 = cb_model_Fp.expend_model(cb_model_Bh, name='CB model 23')
# cb_model_23.resources_group = [5, 4]
cb_model_23.experiment_data_df = FP_BH_experimet_data.copy()
cb_model_23.initial_mets = cb_model_23.experiment_data_df[mets_name].values[0, :]

cb_model_13 = cb_model_Ri.expend_model(cb_model_Bh, name='CB model 13')
# cb_model_13.resources_group = [4, 4]
cb_model_13.experiment_data_df = RI_BH_experimet_data.copy()
cb_model_13.initial_mets = cb_model_13.experiment_data_df[mets_name].values[0, :]

cb_model_123 = cb_model_12.expend_model(cb_model_Bh, name='CB model 123')


# cb_model_123.resources_group = [4, 5, 4]


# cb_model_123.initial_mets = initial_mets123

### functions
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
    c_fru = x[0]
    c_Ri = x[1]
    c_Fp = x[2]
    c_Bh = x[3]
    c_ac = x[4]
    c_for = x[5]
    c_but = x[6]
    # Michaelis–Menten kinetics basic :r_kin_basic
    # r_kin_basic = [
    #     kmax[0] * c_fru / (K[0] + c_fru) ,
    #     kmax[1] * c_fru / (K[1] + c_fru) ,
    #     kmax[2] * c_fru / (K[2] + c_fru) ,
    #     kmax[3] * (c_fru / (K[3] + c_fru)) * (c_ac / (K[3] + c_ac)) ,]
    r_kin_basic_numerator_1 = np.array([c_fru, c_fru, c_fru, c_ac * c_fru])

    r_kin_basic_numerator_2 = np.array([c_fru, c_fru, c_fru, c_fru, c_ac * c_fru])

    r_kin_basic_numerator_3 = np.array([c_fru, c_fru, c_fru, c_for])

    if name == 'CB model 1':
        r_kin_basic_numerator = r_kin_basic_numerator_1
        r_kin_basic_denominator = (K + r_kin_basic_numerator)
        r_kin_basic_denominator[-1] = (K[-1] + c_fru) * (K[-1] + c_ac)

    elif name == 'CB model 2':
        r_kin_basic_numerator = r_kin_basic_numerator_2
        r_kin_basic_denominator = (K + r_kin_basic_numerator)
        r_kin_basic_denominator[-1] = (K[-1] + c_fru) * (K[-1] + c_ac)

    elif name == 'CB model 3':
        r_kin_basic_numerator = r_kin_basic_numerator_3
        r_kin_basic_denominator = (K + r_kin_basic_numerator)
        r_kin_basic_denominator[-1] = (K[-1] + c_fru) * (K[-1] + c_for)

    elif name == 'CB model 12':
        r_kin_basic_numerator = np.append(r_kin_basic_numerator_1, r_kin_basic_numerator_2)
        r_kin_basic_denominator = (K + r_kin_basic_numerator)
        r_kin_basic_denominator[3] = (K[3] + c_fru) * (K[3] + c_ac)
        r_kin_basic_denominator[-1] = (K[-1] + c_fru) * (K[-1] + c_ac)

    elif name == 'CB model 23':
        r_kin_basic_numerator = np.append(r_kin_basic_numerator_2, r_kin_basic_numerator_3)
        r_kin_basic_denominator = (K + r_kin_basic_numerator)
        r_kin_basic_denominator[4] = (K[4] + c_fru) * (K[4] + c_ac)
        r_kin_basic_denominator[-1] = (K[-1] + c_fru) * (K[-1] + c_for)

    elif name == 'CB model 13':
        r_kin_basic_numerator = np.append(r_kin_basic_numerator_1, r_kin_basic_numerator_3)
        r_kin_basic_denominator = (K + r_kin_basic_numerator)
        r_kin_basic_denominator[3] = (K[3] + c_fru) * (K[3] + c_ac)
        r_kin_basic_denominator[-1] = (K[-1] + c_fru) * (K[-1] + c_for)

    elif name == 'CB model 123':
        r_kin_basic_numerator = np.append(r_kin_basic_numerator_1, r_kin_basic_numerator_2)
        r_kin_basic_numerator = np.append(r_kin_basic_numerator, r_kin_basic_numerator_3)
        r_kin_basic_denominator = (K + r_kin_basic_numerator)
        r_kin_basic_denominator[3] = (K[3] + c_fru) * (K[3] + c_ac)
        r_kin_basic_denominator[8] = (K[8] + c_fru) * (K[8] + c_ac)
        r_kin_basic_denominator[-1] = (K[-1] + c_fru) * (K[-1] + c_for)

    else:
        print('Error name')

    r_kin_basic = kmax * r_kin_basic_numerator / r_kin_basic_denominator

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
#     resources_group = self.resources_group
#
#     if name in ['CB model 12']:
#         resources_group = [4, 5]
#     elif name in ['CB model 23']:
#         resources_group = [4, 3]
#     elif name in ['CB model 23']:
#         resources_group = [5, 3]
#     elif name in ['CB model 123']:
#         resources_group = [4, 5, 3]
#     else:
#         resources_group = [-1]
#     for index_i in range(0, len(resources_group)):
#         if index_i == 0:
#             index_strat = 0
#         else:
#             index_strat = resources_group[index_i - 1]
#         index_end = resources_group[index_i]
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
# setattr(Cybernetic_Model_basic, 'cybernetic_var_def', cybernetic_var_def)

sol1 = Cybernetic_Functions.cb_model_simulate(cb_model_1, tspan, draw=True)
sol2 = Cybernetic_Functions.cb_model_simulate(cb_model_2, tspan, draw=True)
sol3 = Cybernetic_Functions.cb_model_simulate(cb_model_3, tspan, draw=True)

sol12 = Cybernetic_Functions.cb_model_simulate(cb_model_12, tspan, draw=True)
sol23 = Cybernetic_Functions.cb_model_simulate(cb_model_23, tspan, draw=True)
# sol123 = Cybernetic_Functions.cb_model_simulate(cb_model_123, tspan, draw=True)

# %% <fitting>
single_simualte = False
bi_simualte = True
weight_of_method_t_f = 0
weights_of_mets = np.array([1, 5, 5, 5, 1, 2, 1]) ** 2
weights_of_model = np.array([2, 1, 1, 3, 3, 3]) ** 2
tri_siumlate = False
n_path1 = 4
n_path2 = 5
if single_simualte:
    cb_model_1.weight_of_method_t_f = weight_of_method_t_f
    cb_model_1.weights_of_mets = weights_of_mets  # /1000
    cb_models = [cb_model_1, ]
    experiment_data_dfs = [cb_model_1.experiment_data_df, ]  # [5],[6],[7],[8]]
    para_to_fit = {'kmax': np.array([[0], [1], [2], [3]]),
                   'K': np.array([])}

    # minimum_1 = Cybernetic_Functions.parameters_fitting(cb_models, experiment_data_dfs,
    #                                                     para_to_fit, tspan, draw=True,
    #                                                     full_output=False,
    #                                                     options={'xtol': 3, 'ftol': 0.01, 'maxiter': 100,
    #                                                              'disp': True})
    #
    # CB_models = Cybernetic_Functions.update_paras_func(minimum_1.x, cb_models, para_to_fit,
    #                                                    retern_model_or_paras='model')
    # sol1 = Cybernetic_Functions.cb_model_simulate(CB_models[0], tspan, draw=True)

    cb_model_2.weight_of_method_t_f = weight_of_method_t_f
    cb_model_2.weights_of_mets = weights_of_mets  # /1000
    cb_models = [cb_model_2, ]
    experiment_data_dfs = [cb_model_2.experiment_data_df.loc[0:14], ]  # [5],[6],[7],[8]]
    para_to_fit = {'kmax': np.array([[0], [1], [2], [3], [4]]),
                   'K': np.array([])}

    # minimum_2 = Cybernetic_Functions.parameters_fitting(cb_models, experiment_data_dfs,
    #                                                     para_to_fit, tspan, draw=True,
    #                                                     full_output=False,
    #                                                     options={'xtol': 3, 'ftol': 0.01, 'maxiter': 100,
    #                                                              'disp': True})
    #
    # CB_models = Cybernetic_Functions.update_paras_func(minimum_2.x, cb_models, para_to_fit,
    #                                                    retern_model_or_paras='model')
    # sol2 = Cybernetic_Functions.cb_model_simulate(CB_models[0], tspan, draw=True)

    cb_model_3.weight_of_method_t_f = weight_of_method_t_f
    # cb_model_3.weights_of_mets = weights_of_mets    # /1000
    cb_models = [cb_model_3, ]
    experiment_data_dfs = [cb_model_3.experiment_data_df.loc[0:15], ]  # [5],[6],[7],[8]]
    para_to_fit = {'kmax': np.array([[0], [1], [2], [3], ]),
                   'K': np.array([])}

    minimum_3 = Cybernetic_Functions.parameters_fitting(cb_models, experiment_data_dfs,
                                                        para_to_fit, tspan, draw=True,
                                                        full_output=False,
                                                        options={'xtol': 3, 'ftol': 0.01, 'maxiter': 200,
                                                                 'disp': True})

    CB_models = Cybernetic_Functions.update_paras_func(minimum_3.x, cb_models, para_to_fit,
                                                       retern_model_or_paras='model')
    sol3 = Cybernetic_Functions.cb_model_simulate(CB_models[0], tspan, draw=True)

if bi_simualte:
    cb_models = [cb_model_1, cb_model_2, cb_model_3,
                 cb_model_12, cb_model_13, cb_model_23]
    for index_i in range(len(cb_models)):
        model_i = cb_models[index_i]
        model_i.weight_of_method_t_f = weight_of_method_t_f
        model_i.weights_of_mets = weights_of_mets * weights_of_model[index_i]
    # cb_models[-1].weights_of_mets = np.array([4,4,500,500,4,16,4])

    experiment_data_dfs = [cb_model_1.experiment_data_df, cb_model_2.experiment_data_df.loc[0:14],
                           cb_model_3.experiment_data_df.loc[0:13],
                           cb_model_12.experiment_data_df, cb_model_13.experiment_data_df,
                           cb_model_23.experiment_data_df]

    para_to_fit = {'kmax': np.array([
        [0, np.nan, np.nan, 0, 0, np.nan, ],
        [1, np.nan, np.nan, 1, 1, np.nan, ],
        [2, np.nan, np.nan, 2, 2, np.nan, ],
        [3, np.nan, np.nan, 3, 3, np.nan, ],
        [np.nan, 0, np.nan, n_path1 + 0, np.nan, 0, ],
        [np.nan, 1, np.nan, n_path1 + 1, np.nan, 1, ],
        [np.nan, 2, np.nan, n_path1 + 2, np.nan, 2, ],
        [np.nan, 3, np.nan, n_path1 + 3, np.nan, 3, ],
        [np.nan, 4, np.nan, n_path1 + 4, np.nan, 4, ],
        [np.nan, np.nan, 0, np.nan, n_path1 + 0, n_path2 + 0, ],
        [np.nan, np.nan, 1, np.nan, n_path1 + 1, n_path2 + 1, ],
        [np.nan, np.nan, 2, np.nan, n_path1 + 2, n_path2 + 2, ],
        [np.nan, np.nan, 3, np.nan, n_path1 + 3, n_path2 + 3, ],
    ]),

        'K': np.array([
            [0, np.nan, np.nan, 0, 0, np.nan, ],
            [1, np.nan, np.nan, 1, 1, np.nan, ],
            [2, np.nan, np.nan, 2, 2, np.nan, ],
            [3, np.nan, np.nan, 3, 3, np.nan, ],
            [np.nan, 0, np.nan, n_path1 + 0, np.nan, 0, ],
            [np.nan, 1, np.nan, n_path1 + 1, np.nan, 1, ],
            [np.nan, 2, np.nan, n_path1 + 2, np.nan, 2, ],
            [np.nan, 3, np.nan, n_path1 + 3, np.nan, 3, ],
            [np.nan, 4, np.nan, n_path1 + 4, np.nan, 4, ],
            [np.nan, np.nan, 0, np.nan, n_path1 + 0, n_path2 + 0, ],
            [np.nan, np.nan, 1, np.nan, n_path1 + 1, n_path2 + 1, ],
            [np.nan, np.nan, 2, np.nan, n_path1 + 2, n_path2 + 2, ],
            [np.nan, np.nan, 3, np.nan, n_path1 + 3, n_path2 + 3, ],
        ])}

    minimum = Cybernetic_Functions.parameters_fitting(cb_models, experiment_data_dfs,
                                                      para_to_fit, tspan, draw=True,
                                                      full_output=False,
                                                      options={'xtol': 3, 'ftol': 0.01, 'maxiter': 300,
                                                               'disp': True})

    CB_models = Cybernetic_Functions.update_paras_func(minimum.x, cb_models, para_to_fit,
                                                       retern_model_or_paras='model')
    sol1 = Cybernetic_Functions.cb_model_simulate(CB_models[0], tspan, draw=True)
    sol2 = Cybernetic_Functions.cb_model_simulate(CB_models[1], tspan, draw=True)
    sol3 = Cybernetic_Functions.cb_model_simulate(CB_models[2], tspan, draw=True)
    sol12 = Cybernetic_Functions.cb_model_simulate(CB_models[3], tspan, draw=True)
    sol13 = Cybernetic_Functions.cb_model_simulate(CB_models[4], tspan, draw=True)
    sol14 = Cybernetic_Functions.cb_model_simulate(CB_models[5], tspan, draw=True)
if tri_siumlate:
    pass
    # [0,	np.nan,	np.nan,	0,	0,	np.nan,	0,	]
    # [1,	np.nan,	np.nan,	1,	1,	np.nan,	1,	]
    # [2,	np.nan,	np.nan,	2,	2,	np.nan,	2,	]
    # [3,	np.nan,	np.nan,	3,	3,	np.nan,	3,	]
    # [np.nan,	0,	np.nan,	n_path1+0,	np.nan,	0,	4,	]
    # [np.nan,	1,	np.nan,	n_path1+1,	np.nan,	1,	5,	]
    # [np.nan,	2,	np.nan,	n_path1+2,	np.nan,	2,	6,	]
    # [np.nan,	3,	np.nan,	n_path1+3,	np.nan,	3,	7,	]
    # [np.nan,	4,	np.nan,	n_path1+4,	np.nan,	4,	8,	]
    # [np.nan,	np.nan,	0,	np.nan,	n_path1+0,	n_path2+0,	9,	]
    # [np.nan,	np.nan,	1,	np.nan,	n_path1+1,	n_path2+1,	10,	]
    # [np.nan,	np.nan,	2,	np.nan,	n_path1+2,	n_path2+2,	11,	]
    # [np.nan,	np.nan,	3,	np.nan,	n_path1+3,	n_path2+3,	12,	]

# # %% <plot cybernetic model result>


# np.savetxt('sol_Rint.txt', sol1, delimiter=',')
