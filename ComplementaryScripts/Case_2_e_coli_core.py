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

import ConvexHull_yield
import GEM2pathways
import Cybernetic_Functions
import cobra
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull

os.chdir('../ComplementaryData/')

# %% <method EFMs > can not calculate
print('\n---------- Loading EFMs ... ---------- ')
core_EFMs_z = np.genfromtxt('Case2_1_ecoli_core/e_coli_core_EFMs_standardized.csv', delimiter = ',')
EFMs_all_points = core_EFMs_z.T
EFMs_all_points = EFMs_all_points[:,1:]   # modes x reas
experiment_datas = []
print('EFMs:', EFMs_all_points.shape)

# %% <method EFVs > can not calculate
print('\n---------- Loading EFVs ... ---------- ')
core_EFVs_z = np.genfromtxt('Case2_1_ecoli_core/e_coli_core_EFVs_standardized.csv', delimiter = ',')
EFVs_all_points = core_EFVs_z.T
EFVs_all_points = EFVs_all_points[:,1:]   # modes x reas
experiment_datas = []
print('EFVs:', EFVs_all_points.shape)

#%% <FBA mode>
print('\n---------- Loading FBA modes ... ---------- ')
core_FBAMs_z = np.genfromtxt('Case2_1_ecoli_core/ecoli_core_FBAMs_standardized.csv', delimiter = ',')

FBAMs_all_points = core_FBAMs_z.T  # modes x rea
FBAMs_all_points = FBAMs_all_points[:, 1:]
experiment_datas = []
print('FBA models:', FBAMs_all_points.shape)

#%% <Our method>
print('\n---------- Caculating Yield space by this method ... ---------- ')
e_coli_core = cobra.io.read_sbml_model('../ComplementaryData/Case2_1_ecoli_core/e_coli_core.xml')
#e_coli_core.reactions.get_by_id('EX_o2_e').bounds = (0.0,1000.0)
# glc as sub
model = e_coli_core.copy()

production_rea_ids_x = ['BIOMASS_Ecoli_core_w_GAM', ]
production_rea_ids_y = ['EX_ac_e', 'EX_for_e', 'EX_etoh_e', 'EX_lac__D_e', 'EX_succ_e', ]
carbon_source_rea_id = 'EX_glc__D_e'

steps = 10
carbon_uptake_direction = -1
model.reactions.get_by_id(carbon_source_rea_id).bounds = (-10, -0.001)

# all modes
yield_normalized_df, fluxes_all, hull_index_all = GEM2pathways.get_yield_space_multi(model, production_rea_ids_x,
                                                                                     production_rea_ids_y,
                                                                                     carbon_source_rea_id,
                                                                                     steps=steps,
                                                                                     carbon_uptake_direction=carbon_uptake_direction,
                                                                                     draw=True)
yield_normalized_df.to_csv('Case2_1_ecoli_core/ecoli_core_our_yield_normalized_df.csv', index=True, sep=',')
yield_normalized_df_ = pd.read_csv('Case2_1_ecoli_core/ecoli_core_our_yield_normalized_df.csv', sep=',',
                                   index_col=0, )
yield_normalized_df_hull = yield_normalized_df_[yield_normalized_df_.columns[hull_index_all]]
yield_normalized_df_hull.to_csv('Case2_1_ecoli_core/ecoli_core_our_yield_normalized_df_hull.csv', index=True,
                                sep=',')
yield_normalized_df_hull_ = pd.read_csv('Case2_1_ecoli_core/ecoli_core_our_yield_normalized_df_hull.csv', sep=',',
                                        index_col=0, )

our_all_points = yield_normalized_df_hull_.values.T
our_all_points = our_all_points[abs(our_all_points[:, 0]) > 1e-10, :]
our_all_points = our_all_points[:, 1:]  # exclude the carbon colume
print('Our method:', our_all_points.shape)

# ac as sub: Note limtied by experiment data
model = e_coli_core.copy()

model.reactions.get_by_id('EX_glc__D_e').bounds = (0, 1000)
model.reactions.get_by_id('EX_ac_e').bounds = (-10, -0.001)
p1_yield_value_max, fluxes_opt_max = GEM2pathways.get_yield_opt(model, 'BIOMASS_Ecoli_core_w_GAM', 'EX_ac_e',
                                                                max_or_min='max',
                                                                carbon_uptake_direction=carbon_uptake_direction)
p1_yield_value_min, fluxes_opt_min = GEM2pathways.get_yield_opt(model, 'BIOMASS_Ecoli_core_w_GAM', 'EX_ac_e',
                                                                max_or_min='min',
                                                                carbon_uptake_direction=-1)

additional_points1 = fluxes_opt_max[['EX_glc__D_e', 'BIOMASS_Ecoli_core_w_GAM', 'EX_ac_e']].values
additional_points2 = fluxes_opt_min[['EX_glc__D_e', 'BIOMASS_Ecoli_core_w_GAM', 'EX_ac_e']].values
additional_points = np.insert(additional_points1, 3, additional_points2, axis=0)
additional_points = additional_points.reshape((2, 3)).T
additional_points = additional_points[:, :] / abs(additional_points[-1, :])

our_all_points = yield_normalized_df_hull_.values.T
our_all_points = our_all_points[ abs(our_all_points[:,0]) > 1e-10,:]
our_all_points = our_all_points[:,1:]       #exclude the carbon colume
print('Our method:', our_all_points.shape)


# %% <MYA >

print('\n---------- Loading experiment data ... ---------- ')
experiment_data_df = pd.read_csv('Case2_1_ecoli_core/Ecoli_aerobic_experiment_data.txt', delimiter='\t', header=0)

data_cloum_name = ['glc', 'biomass', 'ac']
experiment_data_df_trimed = experiment_data_df[data_cloum_name]
experiment_data_df_trimed_values = experiment_data_df_trimed.values[:, :] - experiment_data_df_trimed.values[0, :]
experiment_data_df_trimed_values = experiment_data_df_trimed_values[1:, :]
experiment_data_df_trimed_values = experiment_data_df_trimed_values.T
experiment_data_df_trimed_values = experiment_data_df_trimed_values[:, :] / abs(experiment_data_df_trimed_values[0, :])
experiment_data_df_trimed_values = experiment_data_df_trimed_values.T
experiment_data_df_trimed_values = experiment_data_df_trimed_values[:, 1:]

experiment_datas = []  # TODO experiment data!!!
for i in range(0, experiment_data_df_trimed_values.shape[0]):
    experiment_datas.append(list(experiment_data_df_trimed_values[i, :]) + [''] * 4)

qhull_options = 'QJ Qx A0.9999999'  # 'Qt QJ Pp Qw Qx'      #'QG0'mean expect point 0(index)
cutoff_persent = 1  # can't reduce much points because too many dimasions....
# Fixed (Fixme check the hull_all ??? 169 points??? to many !!! why ??? Decimals??)
# Note not Decimals problem try other reason
# decimals = 3
# EFMs_all_points_ = np.around(EFMs_all_points,decimals,)
# FBAMs_all_points_ = np.around(FBAMs_all_points,decimals,)
# our_all_points_ = np.around(our_all_points,decimals,)

# EFMs_indexes, EFMs_weights, EFMs_estimated_datas, EFMs_in_hulls = ConvexHull_yield.pipeline_mya(EFMs_all_points,
#                                                                                                 experiment_datas,
#                                                                                                 qhull_options=qhull_options,
#                                                                                                 method=1,
#                                                                                                 cutoff_persent=cutoff_persent)
#
# EFVs_indexes, EFVs_weights, EFVs_estimated_datas, EFVs_in_hulls = ConvexHull_yield.pipeline_mya(EFVs_all_points,
#                                                                                                 experiment_datas,
#                                                                                                 qhull_options=qhull_options,
#                                                                                                 method=1,
#                                                                                                 cutoff_persent=cutoff_persent)

# FBAMs_indexes, FBAMs_weights, FBAMs_estimated_datas, FBAMs_in_hulls = ConvexHull_yield.pipeline_mya(FBAMs_all_points,
#                                                                                                     experiment_datas,
#                                                                                                     qhull_options=qhull_options,
#                                                                                                     method=1,
#                                                                                                     cutoff_persent=cutoff_persent)

our_indexes, our_weights, our_estimated_datas, our_in_hulls = ConvexHull_yield.pipeline_mya(our_all_points,
                                                                                            experiment_datas,
                                                                                            qhull_options=qhull_options,
                                                                                            method=1,
                                                                                            cutoff_persent=cutoff_persent)


#%% <plot initial yield>: TODO not divided by carbon, so not yield


yield_rea_ids_name = ['Acetate', 'Formate', 'Ethanol', 'Lacate', 'Succinate']  #,


def trim_axs(axs, N):
    """little helper to massage the axs list to have correct length..."""
    axs = axs.flat
    for ax in axs[N:]:
        ax.remove()
    return axs[:N]

cols = 3
rows = len(yield_rea_ids_name) // cols + 1
figsize = (10, 8)
fig, axs = plt.subplots(cols, rows,figsize=figsize)
axs = trim_axs(axs, len(yield_rea_ids_name))

for ax , exmet_reaction in zip(axs,yield_rea_ids_name):
    index = yield_rea_ids_name.index(exmet_reaction)+1

    xy_EFMs = EFMs_all_points[:, [0,index]]
    xy_EFVs = EFVs_all_points[:, [0,index]]
    xy_FBAMs = FBAMs_all_points[:, [0, index]]
    xy_our = our_all_points[:, [0,index]]

    colors_list = ['blue','teal','tab:red','tab:orange']

    points_EFMs = ax.plot(xy_EFMs[:,0], xy_EFMs[:,1],       '.', markerfacecolor='none', color=colors_list[0],
                    alpha = 0.5, markersize=3)
    points_EFVs = ax.plot(xy_EFVs[:,0], xy_EFVs[:,1],       '.', markerfacecolor='none', color=colors_list[1],
                    alpha = 0.5, markersize=3)

    points_our = ax.plot(xy_our[:,0], xy_our[:,1],          '^', markerfacecolor='none', color=colors_list[2],
                         alpha=0.8, markersize=6)
    points_FBAMs = ax.plot(xy_FBAMs[:, 0], xy_FBAMs[:, 1], 's', color=colors_list[3],
                           alpha=1, markersize=8)

    draw_hull = True

    if draw_hull:
        hull_EFMs = ConvexHull(xy_EFMs,qhull_options = qhull_options )
        hull_EFVs = ConvexHull(xy_EFVs, qhull_options=qhull_options)
        hull_FBAMs = ConvexHull(xy_FBAMs, qhull_options=qhull_options)
        hull_our = ConvexHull(xy_our,qhull_options = qhull_options )

        for simplex in hull_EFMs.simplices:
            line_EFMs = ax.plot(xy_EFMs[simplex, 0], xy_EFMs[simplex, 1], 'x-', markerfacecolor='none', color=colors_list[0],
                            alpha= 0.5 ,  markersize=3,lw=2)

        for simplex in hull_EFVs.simplices:
            line_EFVs = ax.plot(xy_EFVs[simplex, 0], xy_EFVs[simplex, 1], 'v-', markerfacecolor='none', color=colors_list[1],
                            alpha = 0.5, markersize=3 ,lw=2)

        for simplex in hull_our.simplices:
            line_our = ax.plot(xy_our[simplex, 0], xy_our[simplex, 1], '^-',markerfacecolor='none', color=colors_list[2],
                            alpha = 0.5, markersize=5,lw=2)

        for simplex in hull_FBAMs.simplices:
            line_FBAMs = ax.plot(xy_FBAMs[simplex, 0], xy_FBAMs[simplex, 1], 'o-', markerfacecolor='none',
                                 color=colors_list[3],
                                 alpha=0.5, markersize=5, lw=2)

    else:
        line_EFMs = points_EFMs
        line_EFVs = points_EFVs
        line_our = points_our
        line_FBAMs = points_FBAMs

    ax.set_ylabel(yield_rea_ids_name[index-1]+' / Glc',fontsize = 18)

ax.set_xlabel('Yield Biomass/Glucose', fontsize=18)
fig.legend((line_EFMs[0], line_EFVs[0], line_our[0], line_FBAMs[0]), ('EFMs', 'EFVs', 'This study', 'FBA modes'),
           bbox_to_anchor=(0.6, 0.3), loc='upper left', borderaxespad=0., prop={'size': 18})
fig.show()

# fig.savefig('Case2_1_ecoli_core/4.png',format ='png')


# %% plot ac

yield_rea_ids_name = ['Acetate', 'Formate', 'Ethanol', 'Lacate', 'Succinate']
yield_rea_ids_name = ['Acetate']

fig = plt.figure()
ax = fig.add_subplot(111)

cols = 1
rows = len(yield_rea_ids_name) // cols + 1
figsize = (10, 8)
# ax= fig.subplots(111)

index = 1


xy_EFMs = EFMs_all_points[:, [0,index]]
xy_EFVs = EFVs_all_points[:, [0,index]]
xy_FBAMs = FBAMs_all_points[:, [0, index]]
xy_our = our_all_points[:, [0,index]]

colors_list = ['blue','teal','tab:red','tab:orange']

points_EFMs = ax.plot(xy_EFMs[:,0], xy_EFMs[:,1],       '.', markerfacecolor='none', color=colors_list[0],
                alpha = 0.5, markersize=3)
points_EFVs = ax.plot(xy_EFVs[:,0], xy_EFVs[:,1],       '.', markerfacecolor='none', color=colors_list[1],
                alpha = 0.5, markersize=3)

points_our = ax.plot(xy_our[:,0], xy_our[:,1],          '^', markerfacecolor='none', color=colors_list[2],
                     alpha=0.8, markersize=6)
points_FBAMs = ax.plot(xy_FBAMs[:, 0], xy_FBAMs[:, 1], 's', color=colors_list[3],
                       alpha=1, markersize=8)

draw_hull = True

if draw_hull:
    hull_EFMs = ConvexHull(xy_EFMs,qhull_options = qhull_options )
    hull_EFVs = ConvexHull(xy_EFVs, qhull_options=qhull_options)
    hull_FBAMs = ConvexHull(xy_FBAMs, qhull_options=qhull_options)
    hull_our = ConvexHull(xy_our,qhull_options = qhull_options )

    for simplex in hull_EFMs.simplices:
        line_EFMs = ax.plot(xy_EFMs[simplex, 0], xy_EFMs[simplex, 1], 'x-', markerfacecolor='none', color=colors_list[0],
                        alpha= 0.5 ,  markersize=3,lw=2)

    for simplex in hull_EFVs.simplices:
        line_EFVs = ax.plot(xy_EFVs[simplex, 0], xy_EFVs[simplex, 1], 'v-', markerfacecolor='none', color=colors_list[1],
                        alpha = 0.5, markersize=3 ,lw=2)

    for simplex in hull_our.simplices:
        line_our = ax.plot(xy_our[simplex, 0], xy_our[simplex, 1], '^-',markerfacecolor='none', color=colors_list[2],
                        alpha = 0.5, markersize=5,lw=2)

    for simplex in hull_FBAMs.simplices:
        line_FBAMs = ax.plot(xy_FBAMs[simplex, 0], xy_FBAMs[simplex, 1], 'o-', markerfacecolor='none',
                             color=colors_list[3],
                        alpha = 0.5, markersize=10,lw=2)

else:
    line_EFMs = points_EFMs
    line_EFVs = points_EFVs
    line_our = points_our
    line_FBAMs = points_FBAMs

ax.set_ylabel(yield_rea_ids_name[index-1]+' / Glc',fontsize = 18)

ax.set_xlabel('Yield Biomass/Glucose', fontsize=18)
fig.legend((line_EFMs[0], line_EFVs[0], line_our[0], line_FBAMs[0]), ('EFMs', 'EFVs', 'This study', 'FBA modes'),
           bbox_to_anchor=(0.6, 0.95), loc='upper left', borderaxespad=0., prop={'size': 16})
fig.show()

# fig.savefig('Case2_1_ecoli_core/0_ac.png',format ='png')


# %%  <Cybernetic model simulations>:


tStart = 0.0  # DefineTime
tStop = 10
tStep = 0.1
tspan = np.linspace(tStart, tStop, (tStop - tStart) / tStep)

# matrix Z: reactions x pathways; and Smz metabolites x pathways , S metabolites x reactions : Smz = Sm @ Z
final_index = our_indexes[-1][-1]
Smz = yield_normalized_df_hull_.values[:, final_index]
Smz = Smz[[0, 1, 2], :]
# Smz[0, :] = -Smz[0, :]
# path_ac = np.array([0, 0, -1])  # Note !!! this pathway is nessary  to simulate the for experimentdata
# Smz = np.column_stack((Smz,path_for))
Smz = np.insert(Smz, 3, values=additional_points.T, axis=1)
# metabolites and pathways number
(n_mets, n_path) = Smz.shape

# experiment data
metabObj = ['glc', 'biomass', 'ac']  # , 'for', 'etoh', 'lac', 'succ'

# initial metabolites at tome 0, t0: initial_mets.shape = (n_mets,)
initial_mets = experiment_data_df[metabObj].values[0, :]
# initial_mets = np.insert(initial_mets, 3, values=np.array([0.01,0.01,0.01,0.01]), axis=0)
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

    1.6641e+02,
    1.6801e+01,
    2.7765e+01,
    7.3644e+00,
    4.0808e-02,

])

# K : n_path
K = np.array([
    3.2575e-03,
    2.4931e+00,
    3.5281e+00,
    5.1913e-03,
    3.8855e-03, ])

# carbon number for each pathways
n_carbon = np.array([6, 6, 6, 2, 2])
Biomass_index = 1
sub_index = 0

ecoli_core_our_cb = Cybernetic_Functions.Cybernetic_Model('CB model for Ecoli core matrix ')
ecoli_core_our_cb.Smz = Smz
ecoli_core_our_cb.x0 = initial_x0
ecoli_core_our_cb.kmax = kmax
ecoli_core_our_cb.K = K
ecoli_core_our_cb.ke = ke
ecoli_core_our_cb.alpha = alpha
ecoli_core_our_cb.beta = beta
ecoli_core_our_cb.n_carbon = n_carbon
ecoli_core_our_cb['sub_index'] = sub_index
ecoli_core_our_cb['Biomass_index'] = Biomass_index

CB_model = ecoli_core_our_cb


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
        kmax[0] * x[sub_index] / (K[0] + x[sub_index]),
        kmax[1] * x[sub_index] / (K[1] + x[sub_index]),
        kmax[2] * x[sub_index] / (K[2] + x[sub_index]),
        kmax[3] * x[2] / (K[3] + x[2]),
        kmax[4] * x[2] / (K[4] + x[2])]

    rM = r_kin_basic * x[n_mets:]

    rE = ke * r_kin_basic / kmax

    rG = Smz[Biomass_index, :] * rM[:]  # 11: biomass index

    return rM, rE, rG


sol = Cybernetic_Functions.cb_model_simulate(ecoli_core_our_cb, tspan, draw=False)

# %% <plot cybernetic model result>ac

# experiment data

fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(111)

for metid in metabObj:
    experiment_p = ax.plot(experiment_data_df['time'], experiment_data_df[metid], 'o--',
                           alpha=0.5, )

for key in range(0, CB_model.Smz.shape[0]):
    model_line = ax.plot(tspan, sol[:, key], color="k", linewidth=2)
fig.show()
