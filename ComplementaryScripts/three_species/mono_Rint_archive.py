#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2/17/20

"""mono_Rint_archive.py
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
import gemstool
from scipy.spatial import ConvexHull
import os

os.chdir('../../ComplementaryData/three_species/')

Rint = cobra.io.load_json_model('Rint_growth.json')

# %%
Rint.reactions.get_by_id('EX_fru_e').bounds = (-10, 1000)
Rint.reactions.get_by_id('EX_ac_e').bounds = (-10, 1000)
Rint.objective = 'biomass'
print(Rint.optimize())

Rint.reactions.get_by_id('EX_glc__D_e').bounds = (0, 1000)
Rint.reactions.get_by_id('EX_fru_e').bounds = (0, 1000)
Rint.reactions.get_by_id('EX_ac_e').bounds = (0, 1000)
Rint.reactions.get_by_id('EX_for_e').bounds = (0, 1000)
Rint.reactions.get_by_id('EX_for_e').bounds = (0, 1000)
Rint.reactions.get_by_id('EX_4abut_e').bounds = (0, 1000)
model = Rint.copy()

# %% <all modes for fru>
model = Rint.copy()
production_rea_ids_x = ['biomass']
production_rea_ids_y = ['EX_ac_e', 'EX_for_e', 'EX_lac__L_e', 'EX_4abut_e']
carbon_source_rea_id = 'EX_fru_e'  # 'EX_ac_e',

steps = 10
carbon_uptake_direction = -1
model.reactions.get_by_id(carbon_source_rea_id).bounds = (-10, -0.001)

yield_normalized_df_fru, fluxes_all, hull_index_all = GEM2pathways.get_yield_space_multi(model, production_rea_ids_x,
                                                                                         production_rea_ids_y,
                                                                                         carbon_source_rea_id,
                                                                                         steps=steps,
                                                                                         carbon_uptake_direction=carbon_uptake_direction,
                                                                                         draw=True)
yield_normalized_df_hull_fru = yield_normalized_df_fru[yield_normalized_df_fru.columns[hull_index_all]]

yield_normalized_df_hull_fru.to_csv('Rint_yield_normalized_df_hull_fru.csv', sep=',', index=True)

# %% <all modes for ac>
model = Rint.copy()
production_rea_ids_x = ['biomass']
production_rea_ids_y = ['EX_for_e', 'EX_lac__L_e', 'EX_4abut_e']
carbon_source_rea_id = 'EX_ac_e'  # 'EX_ac_e',

steps = 10
carbon_uptake_direction = -1
model.reactions.get_by_id(carbon_source_rea_id).bounds = (-10, -0.001)

yield_normalized_df_ac, fluxes_all, hull_index_all = GEM2pathways.get_yield_space_multi(model, production_rea_ids_x,
                                                                                        production_rea_ids_y,
                                                                                        carbon_source_rea_id,
                                                                                        steps=steps,
                                                                                        carbon_uptake_direction=carbon_uptake_direction,
                                                                                        draw=True)
yield_normalized_df_hull_ac = yield_normalized_df_ac[yield_normalized_df_ac.columns[hull_index_all]]

yield_normalized_df_hull_ac.to_csv('Rint_yield_normalized_df_hull_ac.csv', sep=',', index=True)

# %%
yield_normalized_df_hull_fru_ = pd.read_csv('Rint_yield_normalized_df_hull_fru.csv',
                                            index_col=0)

yield_normalized_df_hull_ac_ = pd.read_csv('Rint_yield_normalized_df_hull_ac.csv',
                                           index_col=0)

#
# yield_normalized_df_hull_ac_.loc['EX_fru_e'] = 0
# yield_normalized_df_hull_ = pd.concat([yield_normalized_df_hull_fru_,yield_normalized_df_hull_ac_],axis = 1,sort =True)
# yield_normalized_df_hull_ = yield_normalized_df_hull_.loc[['EX_fru_e','biomass','EX_ac_e','EX_for_e', 'EX_4abut_e','EX_lac__L_e', ]]
#
#
# our_all_points = yield_normalized_df_hull_.values.T
# # our_all_points = our_all_points[abs(our_all_points[:, 0]) > 1e-10, :]
# # our_all_points = our_all_points[:, 1:]  # exclude the carbon colume
#
#
# # our_all_points[:, 1] = our_all_points[:, 1]
#
# print('Our method:', our_all_points.shape)


# %% <MYA fur>

our_all_points = yield_normalized_df_hull_fru_.values.T
our_all_points = our_all_points[abs(our_all_points[:, 0]) > 1e-10, :]
our_all_points = our_all_points[:, 1:]  # exclude the carbon colume

print('\n---------- Loading experiment data ... ---------- ')
# NOTE: according to the our_estimated_datas, set the biomacc_mas_coefficient = 5.0 it means 1 in model , 1 *5.0 in experiment data
coefficient = 1e+9
experiment_data_df = pd.read_csv('Rint_experiment_data2.txt', delimiter='\t', header=0)
data_cloum_name = ['fru', 'biomass', 'ac', 'for', 'but', 'lac', ]  # biomass
experiment_data_df['biomass'] = experiment_data_df['biomass'] / coefficient
experiment_data_df_trimed = experiment_data_df[data_cloum_name]
experiment_data_df_trimed_values = experiment_data_df_trimed.values[:, :] - experiment_data_df_trimed.values[0, :]
experiment_data_df_trimed_values = experiment_data_df_trimed_values[1:, :]
experiment_data_df_trimed_values = experiment_data_df_trimed_values.T
experiment_data_df_trimed_values = experiment_data_df_trimed_values[:, :] / abs(experiment_data_df_trimed_values[0, :])
experiment_data_df_trimed_values = experiment_data_df_trimed_values.T  # TODO
experiment_data_df_trimed_values = experiment_data_df_trimed_values[:, 1:]

experiment_datas = []  # TODO experiment data!!!
for i in range(0, experiment_data_df_trimed_values.shape[0]):
    experiment_data = list(experiment_data_df_trimed_values[i, :])
    # experiment_data[0] = ''
    experiment_datas.append(experiment_data)

qhull_options = 'QJ Qx A0.9999999'
cutoff_persent = 0.99
our_indexes, our_weights, our_estimated_datas, our_in_hulls = ConvexHull_yield.pipeline_mya(our_all_points,
                                                                                            experiment_datas,
                                                                                            qhull_options=qhull_options,
                                                                                            method=1,
                                                                                            cutoff_persent=cutoff_persent)

# our_hull = ConvexHull(our_all_points[:, 1:], qhull_options=qhull_options)
# ConvexHull_yield.point_in_hull(experiment_datas[-1][1:], our_hull, tolerance=1e-12)


# hull_cutoff_index_99 = [ 0,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 13, 14, 15, 16, 17, 18,
#         19, 20, 21, 22, 26, 34]
# hull_active_index = our_indexes[-1][-1]  # list(set(our_indexes[-1][-1]) | set(our_indexes[-1][-2]) | set(our_indexes[-1][-3]))
#
#
# our_all_points_hull_1 = our_all_points[our_indexes[0],:]
# our_all_points_hull_99 = our_all_points[hull_cutoff_index_99,:]
# our_all_points_hull_act = our_all_points[hull_active_index, :]


# %% <part MYA>

coefficient = 1e+11
experiment_data_df = pd.read_csv('Rint_experiment_data2.txt', delimiter='\t', header=0)
data_cloum_name = ['fru', 'biomass', 'ac', 'for', 'but', 'lac', ]  # biomass
experiment_data_df['biomass'] = experiment_data_df['biomass'] / coefficient
experiment_data_df_trimed = experiment_data_df[data_cloum_name]
experiment_data_df_trimed_values = experiment_data_df_trimed.values[:, :] - experiment_data_df_trimed.values[0, :]
experiment_data_df_trimed_values = experiment_data_df_trimed_values[1:, :]
experiment_data_df_trimed_values = experiment_data_df_trimed_values.T
experiment_data_df_trimed_values = experiment_data_df_trimed_values[:, :] / abs(experiment_data_df_trimed_values[0, :])
experiment_data_df_trimed_values = experiment_data_df_trimed_values.T * 0.7  # TODO
# experiment_data_df_trimed_values = experiment_data_df_trimed_values[:, 1:]

experiment_datas = []  # TODO experiment data!!!
for i in range(0, experiment_data_df_trimed_values.shape[0]):
    experiment_data = list(experiment_data_df_trimed_values[i, :])
    # experiment_data[0] = ''
    experiment_datas.append(experiment_data)

qhull_options = 'QJ Qx A0.9999999'
cutoff_persent = 0.99
our_indexes, our_weights, our_estimated_datas, our_in_hulls = ConvexHull_yield.pipeline_mya(our_all_points,
                                                                                            experiment_datas,
                                                                                            qhull_options=qhull_options,
                                                                                            method=1,
                                                                                            cutoff_persent=cutoff_persent)

# our_hull = ConvexHull(our_all_points[:, 1:], qhull_options=qhull_options)
# ConvexHull_yield.point_in_hull(experiment_datas[-1][1:], our_hull, tolerance=1e-12)


hull_cutoff_index_99 = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18,
                        19, 20, 21, 22, 26, 34]
hull_active_index = our_indexes[-1][
    -1]  # list(set(our_indexes[-1][-1]) | set(our_indexes[-1][-2]) | set(our_indexes[-1][-3]))

our_all_points_hull_1 = our_all_points[our_indexes[0], :]
our_all_points_hull_99 = our_all_points[hull_cutoff_index_99, :]
our_all_points_hull_act = our_all_points[hull_active_index, :]

# %% <plot initial yield>: TODO not divided by carbon, so not yield

yield_rea_ids_name = ['Acetate', 'Formate', 'Butyrate', 'Lacate']
markersize_list = [10, 10, 10]
lw_sizelist = [2, 2, 2]


def trim_axs(axs, N):
    """little helper to massage the axs list to have correct length..."""
    axs = axs.flat
    for ax in axs[N:]:
        ax.remove()
    return axs[:N]


cols = 3
rows = len(yield_rea_ids_name) // cols + 1
figsize = (10, 8)
fig, axs = plt.subplots(cols, rows, figsize=figsize)
axs = trim_axs(axs, len(yield_rea_ids_name))

for ax, exmet_reaction in zip(axs, yield_rea_ids_name):
    index = yield_rea_ids_name.index(exmet_reaction) + 2
    xy_our = our_all_points[:, [1, index]]

    colors_list = ['tab:blue', 'tab:orange', 'tab:red', 'tab:orange']

    points_our = ax.plot(xy_our[:, 0], xy_our[:, 1], '.', color=colors_list[0],
                         alpha=0.8, markersize=10)

    points_exp = ax.plot(experiment_datas[-1][0], experiment_datas[-1][index], '^', color='red',
                         alpha=1, markersize=15)
    points_exp = ax.plot(experiment_datas[-2][0], experiment_datas[-2][index], '^', color='red',
                         alpha=1, markersize=15)
    points_exp = ax.plot(experiment_datas[-3][0], experiment_datas[-3][index], '^', color='red',
                         alpha=1, markersize=15)
    # points_FBAem = ax.plot(xy_FBAem[:,0], xy_FBAem[:,1],    's', color=colors_list[3],
    #                 alpha=1,  markersize=8)

    draw_hull = True

    if draw_hull:

        xy_our_hull = our_all_points_hull_1[:, [1, index]]
        ax.plot(xy_our_hull[:, 0], xy_our_hull[:, 1], 'x', markerfacecolor='none', color=colors_list[0],
                alpha=0.8, markersize=markersize_list[0])

        hull_our = ConvexHull(xy_our_hull, qhull_options=qhull_options)

        for simplex in hull_our.simplices:
            line_our_1 = ax.plot(xy_our_hull[simplex, 0], xy_our_hull[simplex, 1], 'x-', markerfacecolor='none',
                                 color=colors_list[0],
                                 alpha=0.5, markersize=10, lw=lw_sizelist[0])

        xy_our_hull = our_all_points_hull_99[:, [1, index]]
        ax.plot(xy_our_hull[:, 0], xy_our_hull[:, 1], 'x', markerfacecolor='none', color=colors_list[1],
                alpha=0.8, markersize=markersize_list[1])

        hull_our = ConvexHull(xy_our_hull, qhull_options=qhull_options)

        for simplex in hull_our.simplices:
            line_our_99 = ax.plot(xy_our_hull[simplex, 0], xy_our_hull[simplex, 1], 'x-', markerfacecolor='none',
                                  color=colors_list[1],
                                  alpha=0.5, markersize=10, lw=lw_sizelist[1])

        xy_our_hull = our_all_points_hull_act[:, [1, index]]
        ax.plot(xy_our_hull[:, 0], xy_our_hull[:, 1], 'x', markerfacecolor='none', color=colors_list[2],
                alpha=0.8, markersize=markersize_list[2])

        hull_our = ConvexHull(xy_our_hull, qhull_options=qhull_options)
        for simplex in hull_our.simplices:
            line_our_ac = ax.plot(xy_our_hull[simplex, 0], xy_our_hull[simplex, 1], 'x-', markerfacecolor='none',
                                  color=colors_list[2],
                                  alpha=0.5, markersize=10, lw=lw_sizelist[2])


    else:
        # line_EFMs = points_EFMs
        # line_EFVs = points_EFVs
        # line_FBAem = points_FBAem
        line_our = points_our

    ax.set_ylabel(yield_rea_ids_name[index - 2] + ' / Glc', fontsize=18)

ax.set_xlabel('Yield Biomass/Glucose', fontsize=18)
# fig.legend((line_EFMs[0],line_EFVs[0],line_our[0],line_FBAem[0]),('EFMs','EFVs','This study','FBA modes'),bbox_to_anchor=(0.6, 0.3), loc='upper left', borderaxespad=0.,prop={'size': 18})
fig.legend((points_our[0], points_exp[0], line_our_1[0], line_our_99[0], line_our_ac[0]),
           ('All pathways', 'Experiment data', 'Convex hull', 'Convex hull %99', 'Convex hull with data'),
           bbox_to_anchor=(0.55, 0.35), loc='upper left', borderaxespad=0., prop={'size': 20})
# fig.savefig('../ComplementaryData/Case3_iML1515/4.png',format ='png')
fig.show()

# %% <Cybernetic model simulations>:

tStart = 0.0  # DefineTime
tStop = 8.5
tStep = 0.1
tspan = np.linspace(tStart, tStop, (tStop - tStart) / tStep)

# matrix Z: reactions x pathways; and Smz metabolites x pathways , S metabolites x reactions : Smz = Sm @ Z
final_index = hull_active_index
Smz = yield_normalized_df_hull_.values[:, final_index]
Smz[1, :] = Smz[1, :] * coefficient
path_for = np.array([0, 0, 0, -1, 0, 0, 0])  # Note !!! this pathway is nessary  to simulate the for experimentdata
# Smz = np.column_stack((Smz,path_for))
Smz = np.insert(Smz, Smz.shape[1], values=path_for, axis=1)
# metabolites and pathways number
(n_mets, n_path) = Smz.shape
np.savetxt('temp.txt', Smz, delimiter=',')
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
    1.8137e-04,
    5.0561e-01,
    3.2933e+00,
    5.6628e+00,
    1.7679e+01,
    1.7825e-01,
    7.1085e+00,
    2.4708e+02
])

# K : n_path
K = np.array([
    4.8584e-05,
    4.9439e-07,
    7.7285e-02,
    3.4581e-04,
    1.2480e+01,
    3.5425e-01,
    1.9143e-03,
    8.5745e+01,
])

# carbon number for each pathways
n_carbon = np.array([6, 6, 6, 6, 6, 6, 6, 0])
Biomass_index = 1
sub_index = 0

ecoli_iML1515_our_cb = Cybernetic_Functions.Cybernetic_Model('CB model for Ecoli iML1515 matrix ')
ecoli_iML1515_our_cb.Smz = Smz
ecoli_iML1515_our_cb.x0 = initial_x0
ecoli_iML1515_our_cb.kmax = kmax
ecoli_iML1515_our_cb.K = K
ecoli_iML1515_our_cb.ke = ke
ecoli_iML1515_our_cb.alpha = alpha
ecoli_iML1515_our_cb.beta = beta
ecoli_iML1515_our_cb.n_carbon = n_carbon
ecoli_iML1515_our_cb['sub_index'] = sub_index
ecoli_iML1515_our_cb['Biomass_index'] = Biomass_index

CB_model = ecoli_iML1515_our_cb


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
        x[3] / (K[7] + x[3])]

    rM = r_kin_basic * x[n_mets:] * kmax

    rE = ke * r_kin_basic

    rG = Smz[Biomass_index, :] * rM[:]  # 11: biomass index

    return rM, rE, rG


sol = Cybernetic_Functions.cb_model_simulate(CB_model, tspan, draw=False)

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
