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

'''
model:

-METEXT
GLU SUC FOR ACT LAC ETH B CO2 H2

R1 : GLU + PEP => G6P + PYR .
R2 : G6P + ATP => 2 T3P + ADP .
R3 : G6P + 6 NAD => T3P + 6 NADH .
R4 : T3P + NAD + ADP => PEP + NADH + ATP .
R5 : PEP + ADP => PYR + ATP .
R6 : PEP + CO2 + 2 NADH => SUC + 2 NAD .
R7 : PYR + CoA => AcCoA + FOR .
R8 : PYR + NADH => LAC + NAD .
R9 : AcCoA + ADP => ACT + CoA + ATP .
R10 : AcCoA + 2 NADH => ETH + CoA + 2 NAD .
R11 : FOR => CO2 + H2 .
R12 : 6.775 G6P + 82.2 ATP + 4.065 NADH => B + 82.2 ADP + 4.065 NAD .
'''

# %% <EFMs > caculated by matlab emftool. and standardized by step0
print('\n---------- Loading EFMs ... ---------- ')
em_z = np.genfromtxt('Case1_ecoli_reduced/ecoli_reduced_EFMs_standardized.csv', delimiter=',')
EFMs_all_points = em_z.T
EFMs_all_points = EFMs_all_points[:, 1:]  # modes x reas
print('EFMs:', EFMs_all_points.shape)

# %% <FBA modes> from reference and standardized by step0
print('\n---------- Loading FBA modes ... ---------- ')
FBAMs_z = np.genfromtxt('Case1_ecoli_reduced/ecoli_reduced_FBAMs_standardized.csv', delimiter=',')
FBAMs_all_points = FBAMs_z.T  # modes x rea
FBAMs_all_points = FBAMs_all_points[:, 1:]
print('FBA models:', FBAMs_all_points.shape)

# %% <Our method>
print('\n---------- Caculating Yield space by this method ... ---------- ')
# load a GEM
ecoli_reduced_model = cobra.io.read_sbml_model('Case1_ecoli_reduced/ecoli_reduced_model.xml')
model = ecoli_reduced_model.copy()

production_rea_ids_x = ['R12', ]
production_rea_ids_y = ['R9', 'EX_FOR', 'R10', 'R8', 'R6']
carbon_source_rea_id = 'R1'
model.reactions.get_by_id(carbon_source_rea_id).bounds = (0.001, 10)  # set carbon lb as -10
steps = 10
carbon_uptake_direction = 1

# all modes
yield_normalized_df, fluxes_all, hull_index_all = GEM2pathways.get_yield_space_multi(model, production_rea_ids_x,
                                                                                     production_rea_ids_y,
                                                                                     carbon_source_rea_id,
                                                                                     steps=steps,
                                                                                     carbon_uptake_direction=carbon_uptake_direction,
                                                                                     draw=True)

yield_normalized_df.to_csv('Case1_ecoli_reduced/ecoli_reduced_our_yield_normalized_df.csv', index=True, sep=',')
yield_normalized_df_ = pd.read_csv('Case1_ecoli_reduced/ecoli_reduced_our_yield_normalized_df.csv', sep=',',
                                   index_col=0, )
yield_normalized_df_hull = yield_normalized_df_[yield_normalized_df_.columns[hull_index_all]]
yield_normalized_df_hull.to_csv('Case1_ecoli_reduced/ecoli_reduced_our_yield_normalized_df_hull.csv', index=True,
                                sep=',')
yield_normalized_df_hull_ = pd.read_csv('Case1_ecoli_reduced/ecoli_reduced_our_yield_normalized_df_hull.csv', sep=',',
                                        index_col=0, )
our_all_points = yield_normalized_df_hull_.values.T
our_all_points = our_all_points[abs(our_all_points[:, 0]) > 1e-10, :]
our_all_points = our_all_points[:, 1:]  # exclude the carbon colume!!!
print('Our method:', our_all_points.shape)

# %% <MYA >
print('\n---------- Loading experiment data ... ---------- ')
experiment_data_df = pd.read_csv('Case1_ecoli_reduced/ecoli_reduce_experiment_data.txt', delimiter='\t', header=0)
data_cloum_name = ['glc', 'biomass', 'ac', 'for', 'etoh', 'lac', 'succ']
experiment_data_df_trimed = experiment_data_df[data_cloum_name]
experiment_data_df_trimed_values = experiment_data_df_trimed.values[:, :] - experiment_data_df_trimed.values[0, :]
experiment_data_df_trimed_values = experiment_data_df_trimed_values[1:, :]
experiment_data_df_trimed_values = experiment_data_df_trimed_values.T
experiment_data_df_trimed_values = experiment_data_df_trimed_values[:, :] / abs(experiment_data_df_trimed_values[0, :])
experiment_data_df_trimed_values = experiment_data_df_trimed_values.T
experiment_data_df_trimed_values = experiment_data_df_trimed_values[:, 1:]

experiment_datas = []  # TODO experiment data!!!
for i in range(0, experiment_data_df_trimed_values.shape[0]):
    experiment_datas.append(experiment_data_df_trimed_values[i, :])

qhull_options = 'QJ Qx A0.9999999'  # 'Qt QJ Pp Qw Qx'      #'QG0'mean expect point 0(index) TODO Qx will lose 7 points check it
cutoff_persent = 1  # can't reduce much points because too many dimasions....
# Fixme check the hull_all ??? 169points??? to many !!! why ??? Decimals??
# Note not Decimals problem try other reason
# decimals = 3
# EFMs_all_points_ = np.around(EFMs_all_points,decimals,)
# FBAMs_all_points_ = np.around(FBAMs_all_points,decimals,)
# our_all_points_ = np.around(our_all_points,decimals,)

EFMs_indexes, EFMs_weights, EFMs_estimated_datas, EFMs_in_hulls = ConvexHull_yield.pipeline_mya(EFMs_all_points,
                                                                                                experiment_datas,
                                                                                                qhull_options=qhull_options,
                                                                                                method=1,
                                                                                                cutoff_persent=cutoff_persent)
# # FBAMs not have enouth point to perform hull
# # FBAMs_indexes, FBAMs_weights, FBAMs_estimated_datas, FBAMs_in_hulls = ConvexHull_yield.pipeline_mya(FBAMs_all_points,
# #                                                                                             experiment_datas,
# #                                                                                             qhull_options=qhull_options,
# #                                                                                             method=1,
# #                                                                                             cutoff_persent=cutoff_persent)
#

# experiment_datas = [['', 0.70959994, 0.65218536, 0.75269061, 0.00598727, 0.10675876]]
our_indexes, our_weights, our_estimated_datas, our_in_hulls = ConvexHull_yield.pipeline_mya(our_all_points,
                                                                                            experiment_datas,
                                                                                            qhull_options=qhull_options,
                                                                                            method=1,
                                                                                            cutoff_persent=cutoff_persent)

# EFMs_hull = ConvexHull(EFMs_all_points,qhull_options = qhull_options )
# # FBAMs_hull = ConvexHull(FBAMs_all_points,qhull_options = qhull_options )
our_hull = ConvexHull(our_all_points[:, 1:], qhull_options=qhull_options)
ConvexHull_yield.point_in_hull(experiment_datas[-1][1:], our_hull, tolerance=1e-12)

# %% <plot initial yield>: TODO not divided by carbon, so not yield
cmap_list = ['Pastel1', 'Pastel2', 'Paired', 'Accent',
             'Dark2', 'Set1', 'Set2', 'Set3',
             'tab10', 'tab20', 'tab20b', 'tab20c']

plt.get_cmap('Pastel1')

yield_rea_ids_name = ['EX_ac_e', 'EX_for_e', 'EX_etoh_e', 'EX_lac__D_e', 'EX_succ_e', ]


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

experiment_points = np.array(experiment_datas)

for ax, exmet_reaction in zip(axs, yield_rea_ids_name):
    index = yield_rea_ids_name.index(exmet_reaction) + 1

    xy_EFMs = EFMs_all_points[:, [0, index]]
    xy_FBAMs = FBAMs_all_points[:, [0, index]]
    xy_our = our_all_points[:, [0, index]]
    xy_exp = experiment_points[-4:-1, [0, index]]

    hull_EFMs = ConvexHull(xy_EFMs, qhull_options='QJ')
    hull_FBAMs = ConvexHull(xy_FBAMs, qhull_options='QJ')
    hull_our = ConvexHull(xy_our, qhull_options='QJ')

    points_our = ax.plot(xy_our[:, 0], xy_our[:, 1], '^', markerfacecolor='none', color='tab:blue',
                         alpha=1, label='This_methd', markersize=5)
    for simplex in hull_our.simplices:
        line_our = ax.plot(xy_our[simplex, 0], xy_our[simplex, 1], '^--', markerfacecolor='none', color='tab:blue',
                           alpha=0.5, label='This_methd', markersize=5)

    points_FBAMs = ax.plot(xy_FBAMs[:, 0], xy_FBAMs[:, 1], 'o', markerfacecolor='none', color='black',
                           alpha=1, label='FBA_mode', markersize=10)
    for simplex in hull_FBAMs.simplices:
        line_FBAMs = ax.plot(xy_FBAMs[simplex, 0], xy_FBAMs[simplex, 1], 'o--', markerfacecolor='none', color='black',
                             alpha=0.5, label='FBA_mode', markersize=10)

    points_EFMs = ax.plot(xy_EFMs[:, 0], xy_EFMs[:, 1], 'x', markerfacecolor='none', color='tab:orange',
                          alpha=1, label='EFMs', markersize=10)
    for simplex in hull_EFMs.simplices:
        line_EFMs = ax.plot(xy_EFMs[simplex, 0], xy_EFMs[simplex, 1], 'x-', markerfacecolor='none', color='tab:orange',
                            alpha=0.5, label='EFMs', markersize=10)

    points_exp = ax.plot(xy_exp[:, 0], xy_exp[:, 1], '^', color='tab:red',
                         alpha=1, label='experiment data', markersize=11)

    ax.set_ylabel(yield_rea_ids_name[index - 1] + '/Glucose', fontsize=12)

ax.set_xlabel('Yield Biomass/Glucose', fontsize=12)
fig.legend((line_EFMs[0], line_FBAMs[0], line_our[0]), ('EFMs', 'FBA modes', 'This study'), bbox_to_anchor=(0.55, 0.25),
           loc='upper left', borderaxespad=0.)
fig.show()
fig.show()

# %% <Cybernetic model simulations>:


tStart = 0.0  # DefineTime
tStop = 8.5
tStep = 0.1
tspan = np.linspace(tStart, tStop, (tStop - tStart) / tStep)

# matrix Z: reactions x pathways; and Smz metabolites x pathways , S metabolites x reactions : Smz = Sm @ Z
final_index = our_indexes[-1][-1]
Smz = yield_normalized_df_hull_.values[:, final_index]
Smz[0, :] = -Smz[0, :]
path_for = np.array([0, 0, 0, -1, 0, 0, 0])  # Note !!! this pathway is nessary  to simulate the for experimentdata
# Smz = np.column_stack((Smz,path_for))
Smz = np.insert(Smz, 5, values=path_for, axis=1)
# metabolites and pathways number
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
    1.0895e+01,
    2.6936e+01,
    4.1210e-09,
    2.2697e+00,
    8.4682e-07,
    7.9862e+00
])

# K : n_path
K = np.array([
    1.5235e-5,
    1.5181e+01,
    3.5111e+00,
    2.4034e-5,
    6.7265e-5,
    7.8,
])

# carbon number for each pathways
n_carbon = np.array([6, 6, 6, 6, 6, 0])
Biomass_index = 1
sub_index = 0

ecoli_reduced_our_cb = Cybernetic_Functions.Cybernetic_Model('CB model for Ecoli reduced matrix ')
ecoli_reduced_our_cb.Smz = Smz
ecoli_reduced_our_cb.x0 = initial_x0
ecoli_reduced_our_cb.kmax = kmax
ecoli_reduced_our_cb.K = K
ecoli_reduced_our_cb.ke = ke
ecoli_reduced_our_cb.alpha = alpha
ecoli_reduced_our_cb.beta = beta
ecoli_reduced_our_cb.n_carbon = n_carbon
ecoli_reduced_our_cb['sub_index'] = sub_index
ecoli_reduced_our_cb['Biomass_index'] = Biomass_index

CB_model = ecoli_reduced_our_cb


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


    rM = [
        kmax[0] * x[0 + n_mets] * x[sub_index] / (K[0] + x[sub_index]),
        kmax[1] * x[1 + n_mets] * x[sub_index] / (K[1] + x[sub_index]),
        kmax[2] * x[2 + n_mets] * x[sub_index] / (K[2] + x[sub_index]),
        kmax[3] * x[3 + n_mets] * x[sub_index] / (K[3] + x[sub_index]),
        kmax[4] * x[4 + n_mets] * x[sub_index] / (K[4] + x[sub_index]),
        kmax[5] * (x[3] ** 2) / ((K[5] ** 2) + (x[3] ** 2)),
    ]

    rE = ke * rM / kmax / np.array(list(x[n_mets:-1]) + [1])

    rG = Smz[Biomass_index, :] * rM[:]  # 11: biomass index

    return rM, rE, rG


sol = Cybernetic_Functions.cb_model_simulate(ecoli_reduced_our_cb, tspan, draw=False)

# %% <plot cybernetic model result>

# experiment data

fig = plt.figure()
ax = fig.add_subplot(111)

for metid in metabObj:
    experiment_p = ax.plot(experiment_data_df['time'], experiment_data_df[metid], 'o--',
                           alpha=0.5, )

for key in range(0, CB_model.Smz.shape[0]):
    model_line = ax.plot(tspan, sol[:, key], color="k", linewidth=2)
fig.show()
