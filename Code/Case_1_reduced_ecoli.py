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

import cobra
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull

import ConvexHull_yield
import Cybernetic_Functions
import GEM2pathways
import seaborn as sns
import matplotlib

matplotlib.rc('font', family="Arial")
matplotlib.rcParams["font.family"] = 'Arial'  # 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
os.chdir('../Data/')

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

print('\n---------- Caculating Yield space by this method ... ---------- ')
# load a GEM
write_files = False
ecoli_reduced_model = cobra.io.read_sbml_model('Case1_ecoli_reduced/ecoli_reduced_model.xml')
model = ecoli_reduced_model.copy()
production_rea_ids_x = ['R12', ]
production_rea_ids_y = ['R9', 'EX_FOR', 'R10', 'R8', 'R6']
carbon_source_rea_id = 'R1'
model.reactions.get_by_id(carbon_source_rea_id).bounds = (0.1, 10)
# model.reactions.get_by_id('EX_FOR').bounds = (0, 1000)
steps = 10
carbon_uptake_direction = 1

# %% ac and all metabolites, yield space/pathways calculation
fluxes_2d, hull_index_2d = GEM2pathways.get_yield_space_2d(model,
                                                           production_rea_ids_2d=['R12', 'R9'],
                                                           carbon_source_rea_ids_2d=['R1', 'R1'], steps=steps,
                                                           carbon_uptake_direction=-1, draw=True)

# all modes
yield_normalized_df, fluxes_all, hull_index_all = GEM2pathways.get_yield_space_multi(model, production_rea_ids_x,
                                                                                     production_rea_ids_y,
                                                                                     carbon_source_rea_id,
                                                                                     steps=steps,
                                                                                     carbon_uptake_direction=carbon_uptake_direction,
                                                                                     draw=True)
if write_files:
    fluxes_2d.to_csv('Case1_ecoli_reduced/ecoli_reduced_our_fluxes_2d_df.csv', index=True, sep=',')
    yield_normalized_df.to_csv('Case1_ecoli_reduced/ecoli_reduced_our_yield_normalized_df.csv', index=True, sep=',')
yield_normalized_df_ = pd.read_csv('Case1_ecoli_reduced/ecoli_reduced_our_yield_normalized_df.csv', sep=',',
                                   index_col=0, )
yield_normalized_df_hull = yield_normalized_df_[yield_normalized_df_.columns[hull_index_all]]
if write_files:
    yield_normalized_df_hull.to_csv('Case1_ecoli_reduced/ecoli_reduced_our_yield_normalized_df_hull.csv', index=True,
                                    sep=',')
yield_normalized_df_hull_ = pd.read_csv('Case1_ecoli_reduced/ecoli_reduced_our_yield_normalized_df_hull.csv', sep=',',
                                        index_col=0, )

our_key_points = yield_normalized_df_hull_.values.T
our_key_points = our_key_points[abs(our_key_points[:, 0]) > 1e-10, :]
our_key_points = our_key_points[:, 1:]  # exclude the carbon colume!!!
print('Our method:', our_key_points.shape)

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

# %% ConvexHull
qhull_options = 'QJ Qx A0.9999999'  # 'Qt QJ Pp Qw Qx'      #'QG0'mean expect point 0(index) TODO Qx will lose 7 points check it
cutoff_persent = 1  # can't reduce much points because too many dimasions....
# Fixme check the hull_all ??? 169points??? to many !!! why ??? Decimals??
# Note not Decimals problem try other reason
# decimals = 3
# EFMs_all_points_ = np.around(EFMs_all_points,decimals,)
# FBAMs_all_points_ = np.around(FBAMs_all_points,decimals,)
# our_all_points = np.around(our_key_points,decimals,)

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
our_indexes, our_weights, our_estimated_datas, our_in_hulls = ConvexHull_yield.pipeline_mya(our_key_points,
                                                                                            experiment_datas,
                                                                                            qhull_options=qhull_options,
                                                                                            method=1,
                                                                                            cutoff_persent=cutoff_persent)

# EFMs_hull = ConvexHull(EFMs_all_points,qhull_options = qhull_options )
# # FBAMs_hull = ConvexHull(FBAMs_all_points,qhull_options = qhull_options )
our_hull = ConvexHull(our_key_points[:, 1:], qhull_options=qhull_options)
ConvexHull_yield.point_in_hull(experiment_datas[-1][1:], our_hull, tolerance=1e-12)

# %% data for plot
points_2d = fluxes_2d.values
points_2d = points_2d / points_2d[0, :]
points_2d = points_2d[[8, 11], :]
points_2d = points_2d.T
our_all_points = yield_normalized_df_.values.T
our_all_points = our_all_points[:, 1:]
our_key_points = our_key_points

# %%fig 2a
figsize = (4, 2.5)
fig, ax = plt.subplots(figsize=figsize)
colors = sns.color_palette("Set2")
sns.set_style("ticks", )
xy_our = points_2d

hull_our = ConvexHull(xy_our, qhull_options='QJ Qx A0.9999999')
points_our_all = ax.plot(xy_our[:, 1], xy_our[:, 0], 'o', markerfacecolor='none', color='black',
                         alpha=0.3, label='This_methd', markersize=3)
for simplex in hull_our.simplices:
    points_our_key = ax.plot(xy_our[simplex, 1], xy_our[simplex, 0], '+', color='tab:blue',
                             alpha=0.8, label='This_methd', markersize=8, linewidth=.5, )
    line_our = ax.plot(xy_our[simplex, 1], xy_our[simplex, 0], '--', color='tab:blue',
                       alpha=0.8, label='This_methd', markersize=8, linewidth=.5, )

xy_EFMs = EFMs_all_points[:, [0, 1]]
points_EFMs = ax.plot(xy_EFMs[:, 0], xy_EFMs[:, 1], 'x', markerfacecolor='none', color='tab:orange',
                      alpha=0.8, label='EFMs', markersize=7)

# all points:
# xy_our = our_all_points[:, [0, 1]]
# points_our_all = ax.plot(xy_our[:, 0], xy_our[:, 1], 'o', markerfacecolor='none', color='black',
#                          alpha=0.5, label='This_methd', markersize=3)
ax.set_xlim((0.01, 0.03))
ax.set_ylim((-0.1, 1))
ax.set_ylabel('Yield: Acetate/Glucose', fontsize=10, family='Arial', )
ax.set_xlabel('Yield: Biomass/Glucose', fontsize=10, family='Arial', )

plt.yticks(fontname="Arial", fontsize=8)
plt.xticks(fontname="Arial", fontsize=8)

fig.legend((points_our_all[0], points_our_key[0], points_EFMs[0], line_our[0]),
           ('All points', 'Key points', 'EFMs', 'ConvexHull'),
           # bbox_to_anchor=(0.65, -0.12),ncol=2, loc=8,
           loc=4, bbox_to_anchor=(0.93, 0.2),
           # loc=2, bbox_to_anchor=(0.15, 0.93),
           prop={'family': 'Arial', 'size': 8})

fig.tight_layout()
fig.savefig('Case1_ecoli_reduced/fig2a_ac_biomass.pdf', bbox_inches='tight')
fig.savefig('Case1_ecoli_reduced/fig2a_ac_biomass.jpg', bbox_inches='tight')

fig.show()

# %% fig 2b
yield_rea_ids_name = ['EX_etoh_e', 'EX_lac__D_e', 'EX_succ_e', 'EX_for_e', ]
yield_rea_ids_lables = ['Ethanol', 'Lactate', 'Succinate', 'Formate', ]


def trim_axs(axs, N):
    """little helper to massage the axs list to have correct length..."""
    axs = axs.flat
    for ax in axs[N:]:
        ax.remove()
    return axs[:N]


cols = 2
rows = 2
figsize = (4, 2.75)
fig, axs = plt.subplots(cols, rows, figsize=figsize)
axs = trim_axs(axs, len(yield_rea_ids_name))

# experiment_points = np.array(experiment_datas)
for ax, exmet_reaction in zip(axs, yield_rea_ids_name):
    index = yield_rea_ids_name.index(exmet_reaction) + 1
    lable = yield_rea_ids_lables[index - 1]
    if lable == 'Ethanol':
        index = 3
    elif lable == 'Lactate':
        index = 4
    elif lable == 'Succinate':
        index = 5
    elif lable == 'Formate':
        index = 2

    xy_EFMs = EFMs_all_points[:, [0, index]]
    xy_our_all = our_all_points[:, [0, index]]
    xy_our_key = our_key_points[:, [0, index]]
    # xy_exp = experiment_points[-4:-1, [0, index]]

    hull_EFMs = ConvexHull(xy_EFMs, qhull_options='QJ')
    hull_our = ConvexHull(xy_our_all, qhull_options='QJ')

    points_our_all = ax.plot(xy_our_all[:, 0], xy_our_all[:, 1], 'o', markerfacecolor='none', color='black',
                             alpha=0.3, label='This_methd', markersize=2)

    points_our_key = ax.plot(xy_our_key[:, 0], xy_our_key[:, 1], '+', color='tab:blue',
                             alpha=0.8, label='This_methd', markersize=6)

    points_EFMs = ax.plot(xy_EFMs[:, 0], xy_EFMs[:, 1], 'x', color='tab:orange',
                          alpha=0.8, label='EFMs', markersize=5)
    # outline
    for simplex in hull_EFMs.simplices:
        line_EFMs = ax.plot(xy_EFMs[simplex, 0], xy_EFMs[simplex, 1], ':', markerfacecolor='none', color='tab:orange',
                            alpha=0.8, label='EFMs', markersize=5, linewidth=.5, )

    for simplex in hull_our.simplices:
        line_our = ax.plot(xy_our_all[simplex, 0], xy_our_all[simplex, 1], '--', color='tab:blue',
                           alpha=0.8, label='This_methd', markersize=3, linewidth=.5, )

    # points_exp = ax.plot(xy_exp[:, 0], xy_exp[:, 1], '^', color='tab:red',
    #                      alpha=1, label='experiment data', markersize=11)

    ax.set_ylabel(lable + '/Glucose', fontsize=8, family='Arial', labelpad=1)
    ax.set_xlim((0.01, 0.03))
    ax.set_ylim((-0.2, 2))
    if index == 5:
        ax.set_ylim((-0.1, 1))
        ax.set_yticks(np.arange(0, 1.2, 0.2))
        ax.set_yticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

    ax.tick_params(axis="x", labelsize=6, pad=0.5)
    ax.tick_params(axis="y", labelsize=6, pad=0.5),

# fig.legend((points_our_all[0], points_our_key[0], points_EFMs[0], line_our[0]),
#            ('All points', 'Key points', 'EFMs', 'ConvexHull'),
#            bbox_to_anchor=(0.8, -0.22),
#            loc=8, ncol=1, prop={'family': 'Arial', 'size': 8})

fig.tight_layout(pad=1, h_pad=0.5, w_pad=-0.2)
fig.savefig('Case1_ecoli_reduced/fig2b_all_mets_biomass.pdf', bbox_inches='tight', transparent=True)
fig.savefig('Case1_ecoli_reduced/fig2b_all_mets_biomass.svg', bbox_inches='tight', transparent=True)

fig.show()

# %% <plot initial yield>:
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
figsize = (8, 6)
fig, axs = plt.subplots(cols, rows, figsize=figsize)
axs = trim_axs(axs, len(yield_rea_ids_name))

experiment_points = np.array(experiment_datas)

for ax, exmet_reaction in zip(axs, yield_rea_ids_name):
    index = yield_rea_ids_name.index(exmet_reaction) + 1

    xy_EFMs = EFMs_all_points[:, [0, index]]
    xy_FBAMs = FBAMs_all_points[:, [0, index]]
    xy_our = our_key_points[:, [0, index]]
    xy_exp = experiment_points[-4:-1, [0, index]]

    hull_EFMs = ConvexHull(xy_EFMs, qhull_options='QJ')
    hull_FBAMs = ConvexHull(xy_FBAMs, qhull_options='QJ')
    hull_our = ConvexHull(xy_our, qhull_options='QJ')

    points_our = ax.plot(xy_our[:, 0], xy_our[:, 1], '^', markerfacecolor='none', color='tab:blue',
                         alpha=1, label='This_methd', markersize=5)
    for simplex in hull_our.simplices:
        line_our = ax.plot(xy_our[simplex, 0], xy_our[simplex, 1], '^--', markerfacecolor='none', color='tab:blue',
                           alpha=0.5, label='This_methd', markersize=5)

    # points_FBAMs = ax.plot(xy_FBAMs[:, 0], xy_FBAMs[:, 1], 'o', markerfacecolor='none', color='black',
    #                        alpha=1, label='FBA_mode', markersize=10)
    # for simplex in hull_FBAMs.simplices:
    #     line_FBAMs = ax.plot(xy_FBAMs[simplex, 0], xy_FBAMs[simplex, 1], 'o--', markerfacecolor='none', color='black',
    #                          alpha=0.5, label='FBA_mode', markersize=10)

    points_EFMs = ax.plot(xy_EFMs[:, 0], xy_EFMs[:, 1], 'x', markerfacecolor='none', color='tab:orange',
                          alpha=1, label='EFMs', markersize=10)
    for simplex in hull_EFMs.simplices:
        line_EFMs = ax.plot(xy_EFMs[simplex, 0], xy_EFMs[simplex, 1], 'x-', markerfacecolor='none', color='tab:orange',
                            alpha=0.5, label='EFMs', markersize=10)

    # points_exp = ax.plot(xy_exp[:, 0], xy_exp[:, 1], '^', color='tab:red',
    #                      alpha=1, label='experiment data', markersize=11)

    ax.set_ylabel(yield_rea_ids_name[index - 1] + '/Glucose', fontsize=12)

ax.set_xlabel('Yield Biomass/Glucose', fontsize=12)
# fig.legend((line_EFMs[0], line_FBAMs[0], line_our[0]), ('EFMs', 'FBA modes', 'This study'), bbox_to_anchor=(0.55, 0.25),
#            loc='upper left', borderaxespad=0.)
fig.legend((line_EFMs[0], line_our[0]), ('EFMs', 'This study'), bbox_to_anchor=(0.55, 0.25),
           loc='upper left', borderaxespad=0.)
fig.show()

# %% <Cybernetic model simulations>:

yield_normalized_df_hull_ = pd.read_csv('Case1_ecoli_reduced/ecoli_reduced_our_yield_normalized_df_hull.csv', sep=',',
                                        index_col=0, )
final_index = [4, 12, 13, 14, 19]

Smz = yield_normalized_df_hull_.values[:, final_index]
tStart = 0.0  # DefineTime
tStop = 8.5
tStep = 0.1
# tspan = np.linspace(tStart, tStop, (tStop - tStart) / tStep)
tspan = np.linspace(tStart, tStop, int(((tStop - tStart) / tStep) + 1))
# matrix Z: reactions x pathways; and Smz metabolites x pathways , S metabolites x reactions : Smz = Sm @ Z
final_index = our_indexes[-1][-1]
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
initial_enzyme = np.array([0.9] * (n_path - 1) + [1])

# initial data x0 :initial_x0.shape  = (n_mets + n_path,)
initial_x0 = np.concatenate((initial_mets, initial_enzyme))

# Enzyme Rate Parameters: alpha,beta,ke : de/dt =  alpha + rE(ke) * u - (beta + mu) * e
alpha = np.array([0.04] * n_path)
beta = np.array([0.05] * n_path)
ke = np.array([0.620342] * n_path)  # or 0.5

# Metabolites rate Parameters kmax , Ki : dm/dt =  Smz @ V @ rM(kmax,K) * c
# kmax : n_path
kmax = np.array([
    10.9834382998037,
    22.0252764422883,
    2.51869911959873e-09,
    8.63617956141342e-06,
    1.84082444427781e-07,
    28.4203763528915
])

# K : n_path
K = np.array([
    0.00242654588166614,
    9.65979532962462,
    2.91194363535861e-07,
    3.39429549655414,
    4.32918435262398,
    4.75053521070238,
])

# carbon number for each pathways
n_carbon = np.array([6, 6, 6, 6, 6, 0])
Biomass_index = 1
sub_index = 0


# Construct the model:
class Cybernetic_Model_basic(Cybernetic_Functions.Cybernetic_Model):
    pass


ecoli_reduced_our_cb = Cybernetic_Model_basic('CB model for Ecoli reduced matrix ')
ecoli_reduced_our_cb.Smz = Smz
# ecoli_reduced_our_cb.x0 = initial_x0
ecoli_reduced_our_cb.initial_mets = initial_mets
ecoli_reduced_our_cb.initial_enzymes = initial_enzyme
ecoli_reduced_our_cb.kmax = kmax
ecoli_reduced_our_cb.K = K
ecoli_reduced_our_cb.ke = ke
ecoli_reduced_our_cb.alpha = alpha
ecoli_reduced_our_cb.beta = beta
ecoli_reduced_our_cb.n_carbon = n_carbon
ecoli_reduced_our_cb.sub_index = sub_index
ecoli_reduced_our_cb.Biomass_index = Biomass_index
# ecoli_reduced_our_cb.mets_name = mets_name
# ecoli_reduced_our_cb.experiment_data_df = experiment_data_df_1

CB_model = ecoli_reduced_our_cb


def rate_def(self, x):
    name = self.name
    Smz = self.Smz
    kmax = self.kmax
    K = self.K
    ke = self.ke
    Biomass_index = self.Biomass_index
    sub_index = self.sub_index

    (n_mets, n_path) = Smz.shape

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

    if type(Biomass_index) != int and len(Biomass_index) == n_path:
        rG = rM[:] * np.array([Smz[Biomass_index[i], i] for i in np.arange(0, n_path)])
    else:
        rG = rM[:] * Smz[Biomass_index, :]  # biomass index

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

    if CB_model.name in ['CB model for Ecoli reduced matrix ']:
        # print(123)
        u[-1] = 1.0
        v[-1] = 1.0
    return u, v


setattr(Cybernetic_Model_basic, 'rate_def', rate_def)
setattr(Cybernetic_Model_basic, 'cybernetic_var_def', cybernetic_var_def)

sol = Cybernetic_Functions.cb_model_simulate(CB_model, tspan, draw=False)

np.savetxt('Case1_ecoli_reduced/ecoli_reduced_sol.csv', sol, delimiter=',')
sol_ = np.loadtxt('Case1_ecoli_reduced/ecoli_reduced_sol.csv', delimiter=',')
# %% <plot cybernetic model result>

# experiment data

fig = plt.figure(figsize=(5.5, 2.5))
ax = fig.add_subplot(111)
# colors = ['blue', 'teal', 'tab:red', 'tab:orange']
color_list = plt.cm.tab10(np.linspace(0, 1, 11))

for index in range(0, CB_model.Smz.shape[0]):
    if index == 1:
        ax1 = ax.twinx()
        ax1.plot(tspan, sol[:, index], color=color_list[index], linewidth=1, label=metabObj[index])

        experiment_p = ax1.plot(experiment_data_df['time'], experiment_data_df[metabObj[index]], '*',
                                color=color_list[index], alpha=0.8,
                                linewidth=1)
    else:
        ax.plot(tspan, sol[:, index], color=color_list[index], linewidth=1, label=metabObj[index])

        experiment_p = ax.plot(experiment_data_df['time'], experiment_data_df[metabObj[index]], '*',
                               color=color_list[index], alpha=0.8,
                               linewidth=1)

ax.set_xlabel('Time (h)', fontsize=10, family='Arial', )
ax.set_ylabel('Concentration (mM)', fontsize=10, family='Arial', )
ax1.set_ylabel('Biomass (g/L)', fontsize=10, family='Arial', )
plt.yticks(fontname="Arial", fontsize=8)
plt.xticks(fontname="Arial", fontsize=8)
L = fig.legend(loc='lower left', bbox_to_anchor=(1, 0.2), ncol=1, fontsize=8, )
plt.setp(L.texts, family='Arial')
fig.tight_layout()
fig.savefig('Case1_ecoli_reduced/fig2c_simulate.pdf', bbox_inches='tight')
fig.savefig('Case1_ecoli_reduced/fig2c_simulate.jpg', bbox_inches='tight')

fig.show()
