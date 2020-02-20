#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2/17/20

"""mono_Byd_archive.py
:description : script
:param :
:returns:
:rtype:
"""

# !/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2020-02-17

"""mono_Ritint.py
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

Bhyd = cobra.io.load_json_model('Bhyd_growth.json')

# %%
Bhyd.reactions.get_by_id('EX_fru_e').bounds = (-10, 1000)
Bhyd.objective = 'biomass'
print(Bhyd.optimize())

Bhyd.reactions.get_by_id('EX_glc__D_e').bounds = (0, 1000)
Bhyd.reactions.get_by_id('EX_fru_e').bounds = (0, 1000)
Bhyd.reactions.get_by_id('EX_ac_e').bounds = (0, 1000)
Bhyd.reactions.get_by_id('EX_for_e').bounds = (0, 1000)
Bhyd.reactions.get_by_id('EX_for_e').bounds = (0, 1000)
Bhyd.reactions.get_by_id('EX_4abut_e').bounds = (0, 1000)
model = Bhyd.copy()

# %% <all modes for fru>
model = Bhyd.copy()
model.reactions.get_by_id('EX_4abut_e').bounds = (0, 0)
model.reactions.get_by_id('EX_ac_e').bounds = (0, 1000)
model.reactions.get_by_id('EX_for_e').bounds = (0, 0)
production_rea_ids_x = ['biomass']
production_rea_ids_y = ['EX_ac_e', 'EX_for_e', 'EX_4abut_e', 'EX_lac__L_e', ]
carbon_source_rea_id = 'EX_fru_e'  # 'EX_ac_e',

steps = 10
carbon_uptake_direction = -1
model.reactions.get_by_id(carbon_source_rea_id).bounds = (-10, -0.1)

# yield_normalized_df_fru, fluxes_all, hull_index_all = GEM2pathways.get_yield_space_multi(model, production_rea_ids_x,
#                                                                                      production_rea_ids_y,
#                                                                                      carbon_source_rea_id,
#                                                                                      steps=steps,
#                                                                                      carbon_uptake_direction=carbon_uptake_direction,
#                                                                                      draw=True)
# yield_normalized_df_hull_fru = yield_normalized_df_fru[yield_normalized_df_fru.columns[hull_index_all]]
#
# yield_normalized_df_hull_fru.to_csv('Bhyd_yield_normalized_df_hull_fru.csv', sep=',',index=True)


# %%
yield_normalized_df_hull_fru_ = pd.read_csv('Bhyd_yield_normalized_df_hull_fru.csv',
                                            index_col=0)

yield_normalized_df_hull_ = yield_normalized_df_hull_fru_.loc[
    ['EX_fru_e', 'biomass', 'EX_ac_e', 'EX_for_e', 'EX_4abut_e', 'EX_lac__L_e', ]]

our_all_points = yield_normalized_df_hull_.values.T

print('Our method:', our_all_points.shape)

# %%<MYA all>
our_all_points_all = yield_normalized_df_hull_.values.T
our_all_points = yield_normalized_df_hull_.loc[['EX_fru_e', 'EX_ac_e', 'EX_4abut_e', 'EX_lac__L_e', ]].values.T
experiment_data_df = pd.read_csv('Bhyd_experiment_data.txt', delimiter='\t', header=0)

data_cloum_name = ['fru', 'ac', 'but', 'lac', ]  # biomass
experiment_data_df_trimed = experiment_data_df[data_cloum_name]
experiment_data_df_trimed_values = experiment_data_df_trimed.values[:, :] - experiment_data_df_trimed.values[0, :]
experiment_data_df_trimed_values = experiment_data_df_trimed_values[1:, :]
experiment_data_df_trimed_values = experiment_data_df_trimed_values.T
experiment_data_df_trimed_values = experiment_data_df_trimed_values[:, :] / abs(experiment_data_df_trimed_values[0, :])
experiment_data_df_trimed_values = experiment_data_df_trimed_values.T  # TODO
# experiment_data_df_trimed_values = experiment_data_df_trimed_values[:, 1:]
experiment_datas = []  # TODO experiment data!!!
for i in range(0, experiment_data_df_trimed_values.shape[0]):
    experiment_data = list(experiment_data_df_trimed_values[i, :])
    # experiment_data[0] = ''
    experiment_datas.append(experiment_data)

qhull_options = 'QJ Qx A0.9999999'
cutoff_persent = 0.90
# our_indexes, our_weights, our_estimated_datas, our_in_hulls = ConvexHull_yield.pipeline_mya(our_all_points,
#                                                                                             experiment_datas,
#                                                                                             qhull_options=qhull_options,
#                                                                                             method=1,
#                                                                                             cutoff_persent=cutoff_persent)
# # hull_active_index = our_indexes[-1][-1]  # list(set(our_indexes[-1][-1]) | set(our_indexes[-1][-2]) | set(our_indexes[-1][-3]))


hull_active_index = [0, 1, 12, 21, 23]

# %% <Cybernetic model simulations>:

print('CB modeing')
tStart = 0.0  # DefineTime
tStop = 30
tStep = 0.5
tspan = np.linspace(tStart, tStop, (tStop - tStart) / tStep)

# matrix Z: reactions x pathways; and Smz metabolites x pathways , S metabolites x reactions : Smz = Sm @ Z
final_index = hull_active_index
Smz = yield_normalized_df_hull_.values[:, final_index]
# Smz[1, :] = Smz[1, :]
Smz = Smz[[0, 1, 2, 3, 5]]
path_for = np.array([-0.0, 0, 0, -1, 0])  # Note !!! this pathway is nessary  to simulate the for experimentdata
Smz = np.insert(Smz, Smz.shape[1], values=path_for, axis=1)
# Smz[3,:] = Smz[3,:]*0.1 # TODO :for not correct
# Smz[4,:] = Smz[4,:]*6 # TODO :but not correct
# Smz[5,:] = Smz[5,:]*0.1
# Smz = np.column_stack((Smz,path_for))

# metabolites and pathways number
(n_mets, n_path) = Smz.shape
np.savetxt('Smz_Bhyd.txt', Smz, delimiter=',')

# initial metabolites at tome 0, t0: initial_mets.shape = (n_mets,)
initial_mets = experiment_data_df[['fru', 'biomass', 'ac', 'for', 'lac']].values[0, :]
initial_mets = [50, 6.929513 * 5e-3, 1, 50, 0.754884547]
# initial_mets =np.array([50.  , 50.  ,  0.  ,  0.64])
# initial:Enzyme: initial_enzyme.shape = (n_path,)
initial_enzyme = np.array([0.9] * n_path)

# initial data x0 :initial_x0.shape  = (n_mets + n_path,)
initial_x0 = np.concatenate((initial_mets, initial_enzyme))

# Enzyme Rate Parameters: alpha,beta,ke : de/dt =  alpha + rE(ke) * u - (beta + mu) * e
alpha = np.array([0.04] * n_path)
beta = np.array([0.05] * n_path)
ke = np.array([0.5] * n_path)  # or 0.5

# Metabolites rate Parameters kmax , Ki : dm/dt =  Smz @ V @ rM(kmax,K) * c

kmax = np.array([
    4,
    1.1,
    1,
    10.4,
    8.9,
    25,
])
# K = np.array( [10 , 1.4939661  , 1.51974243  ,1.09385034  ,1.43940618 ,11.54933097
#                 ] )
# kmax = np.array( [1.14627263 , 1.17942602 , 0.19772118  ,0.6890593   ,0.79972604 ,11.55897386
#                 ] )

#
# # K : n_path
K = np.array([
    1,
    1,
    1,
    1,
    1,
    1,
])
# [ 16.73117756   3.50063702  -0.55281571   0.47568049  -0.32343297
#   -2.75010462 -26.45844313   2.47333578   0.92252693   0.37722052
#    2.31363903   1.53966147]
# carbon number for each pathways
n_carbon = np.array([6, 6, 6, 6, 6, 1])
Biomass_index = 1
sub_index = 0

Bhyd_our_cb = Cybernetic_Functions.Cybernetic_Model('CB model for Bhyd')
Bhyd_our_cb.Smz = Smz
Bhyd_our_cb.x0 = initial_x0
Bhyd_our_cb.kmax = kmax
Bhyd_our_cb.K = K
Bhyd_our_cb.ke = ke
Bhyd_our_cb.alpha = alpha
Bhyd_our_cb.beta = beta
Bhyd_our_cb.n_carbon = n_carbon
Bhyd_our_cb['sub_index'] = sub_index
Bhyd_our_cb['Biomass_index'] = Biomass_index

CB_model = Bhyd_our_cb


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
        # x[3] / (K[5] + x[3])
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

    if CB_model.name in ['CB model for Bhyd', 'CB model for Ecoli reduced matrix ',
                         'CB model for Ecoli iML1515 matrix ']:
        u[-1] = 1.0
        v[-1] = 1.0
    return u, v


sol = Cybernetic_Functions.cb_model_simulate(CB_model, tspan, draw=True)
np.savetxt('sol_Bhyd.txt', sol, delimiter=',')

metabObj = ['fru', 'biomass', 'ac', 'for', 'lac']

# %% < Fiting >
CB_model['metas_names'] = ['fru', 'biomass', 'ac', 'for', 'lac']
para_to_fit = {'kmax': [0, 1, 2, 3, 4, 5], 'K': []}

# experiment data
experiment_data_df = pd.read_csv('Bhyd_experiment_data.txt', delimiter='\t', header=0)

experiment_data_to_fit = experiment_data_df[['time'] + metabObj].iloc[0:-1]
experiment_data_to_fit['biomass'] = experiment_data_to_fit['biomass'] * 5e-3
minimum = Cybernetic_Functions.parameters_fitting(CB_model, experiment_data_to_fit, para_to_fit, tspan, draw=True)

# %% <plot cybernetic model result>
# experiment data

fig = plt.figure(figsize=(6, 2.5))
ax = fig.add_subplot(111)
colors = ['blue', 'teal', 'tab:red', 'tab:orange']
color_list = plt.cm.tab10(np.linspace(0, 1, 12))
experiment_data_df = pd.read_csv('Bhyd_experiment_data.txt', delimiter='\t', header=0)
experiment_data_df['biomass'] = experiment_data_df['biomass'] * 5e-3
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
