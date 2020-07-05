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

os.chdir('../ComplementaryData/')

iML1515 = cobra.io.read_sbml_model('../ComplementaryData/Case3_iML1515/iML1515.xml')

iML1515.reactions.get_by_id('EX_o2_e').bounds = (0.0, 1000.0)
# % Constrain the phosphotransferase system
# model = changeRxnBounds(model, 'GLCabcpp', -1000, 'l');
# model = changeRxnBounds(model, 'GLCptspp', -1000, 'l');
# model = changeRxnBounds(model, 'GLCabcpp', 1000, 'u');
# model = changeRxnBounds(model, 'GLCptspp', 1000, 'u');
# model = changeRxnBounds(model, 'GLCt2pp', 0, 'b');
iML1515.reactions.get_by_id('GLCabcpp').bounds = (-1000.0,1000.0)
iML1515.reactions.get_by_id('GLCptspp').bounds = (-1000.0,1000.0)
iML1515.reactions.get_by_id('GLCt2pp').bounds = (0.0,0.0)
# iML1515.reactions.get_by_id('EX_o2_e').bounds = (0.0,1000.0)
solution = iML1515.optimize()
print(solution)


model = iML1515.copy()
production_rea_ids_x = ['BIOMASS_Ec_iML1515_core_75p37M', ]
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
# yield_normalized_df_hull = yield_normalized_df[yield_normalized_df.columns[hull_index_all]]
#
#
# yield_normalized_df_hull.to_csv('../ComplementaryData/Case3_iML1515/iML151_yield_normalized_df_hull.csv', sep=',',index=True)


yield_normalized_df_hull_ = pd.read_csv('../ComplementaryData/Case3_iML1515/iML151_yield_normalized_df_hull.csv',
                                        index_col=0)

our_all_points = yield_normalized_df_hull_.values.T
our_all_points = our_all_points[abs(our_all_points[:, 0]) > 1e-10, :]
our_all_points = our_all_points[:, 1:]  # exclude the carbon colume

# NOTE: according to the our_estimated_datas, set the biomacc_mas_coefficient = 5.0 it means 1 in model , 1 *5.0 in experiment data
coefficient = 3

our_all_points[:, 0] = our_all_points[:, 0] * coefficient

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
    experiment_data = list(experiment_data_df_trimed_values[i, :])
    # experiment_data[0] = ''
    experiment_datas.append(experiment_data)

qhull_options = 'QJ Qx A0.9999999'
cutoff_persent = 1
our_indexes, our_weights, our_estimated_datas, our_in_hulls = ConvexHull_yield.pipeline_mya(our_all_points,
                                                                                            experiment_datas,
                                                                                            qhull_options=qhull_options,
                                                                                            method=1,
                                                                                            cutoff_persent=cutoff_persent)

our_hull = ConvexHull(our_all_points[:, 1:], qhull_options=qhull_options)
ConvexHull_yield.point_in_hull(experiment_datas[-1][1:], our_hull, tolerance=1e-12)

hull_cutoff_index_99 = [0, 2, 3, 9, 10, 11, 16, 17, 18, 19, 20, 23, 26, 30, 33, 34, 35]
hull_active_index = our_indexes[-1][
    -2]  # list(set(our_indexes[-1][-1]) | set(our_indexes[-1][-2]) | set(our_indexes[-1][-3]))

# hull_active_index = [3, 17, 20, 30, 32, 33, 34]
our_all_points_hull_1 = our_all_points[our_indexes[0],:]
our_all_points_hull_99 = our_all_points[hull_cutoff_index_99,:]
our_all_points_hull_act = our_all_points[hull_active_index, :]


#%% <plot initial yield>: TODO not divided by carbon, so not yield

yield_rea_ids_name = ['Acetate', 'Formate', 'Ethanol', 'Lacate', 'Succinate']
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
fig, axs = plt.subplots(cols, rows,figsize=figsize)
axs = trim_axs(axs, len(yield_rea_ids_name))

for ax , exmet_reaction in zip(axs,yield_rea_ids_name):
    index = yield_rea_ids_name.index(exmet_reaction)+1

    # xy_EFMs = EFMs_all_points[:, [0,index]]
    # xy_EFVs = EFVs_all_points[:, [0,index]]
    # xy_FBAem = FBAem_all_points[:, [0,index]]
    xy_our = our_all_points[:, [0,index]]

    colors_list = ['tab:blue','tab:orange','tab:red','tab:orange']

    # points_EFMs = ax.plot(xy_EFMs[:,0], xy_EFMs[:,1],       '.', markerfacecolor='none', color=colors_list[0],
    #                 alpha = 0.5, markersize=3)
    # points_EFVs = ax.plot(xy_EFVs[:,0], xy_EFVs[:,1],       '.', markerfacecolor='none', color=colors_list[1],
    #                 alpha = 0.5, markersize=3)

    points_our = ax.plot(xy_our[:,0], xy_our[:,1],          '.', color=colors_list[0],
                    alpha=0.8,  markersize=10)

    points_exp = ax.plot(experiment_datas[-1][0], experiment_datas[-1][index], '^', color='red',
                         alpha=1, markersize=15)
    points_exp = ax.plot(experiment_datas[-2][0], experiment_datas[-2][index], '^', color='red',
                         alpha=1, markersize=15)
    points_exp = ax.plot(experiment_datas[-3][0], experiment_datas[-3][index], '^', color='red',
                    alpha=1,  markersize=15)
    # points_FBAem = ax.plot(xy_FBAem[:,0], xy_FBAem[:,1],    's', color=colors_list[3],
    #                 alpha=1,  markersize=8)

    draw_hull = True

    if draw_hull:

        xy_our_hull = our_all_points_hull_1[:, [0,index]]
        ax.plot(xy_our_hull[:,0], xy_our_hull[:,1],          'x', markerfacecolor='none', color=colors_list[0],
                    alpha=0.8,  markersize=markersize_list[0])

        hull_our = ConvexHull(xy_our_hull,qhull_options = qhull_options )

        for simplex in hull_our.simplices:
            line_our_1 = ax.plot(xy_our_hull[simplex, 0], xy_our_hull[simplex, 1], 'x-',markerfacecolor='none', color=colors_list[0],
                            alpha = 0.5, markersize=10,lw=lw_sizelist[0])


        xy_our_hull = our_all_points_hull_99[:, [0,index]]
        ax.plot(xy_our_hull[:,0], xy_our_hull[:,1],          'x', markerfacecolor='none', color=colors_list[1],
                    alpha=0.8,  markersize=markersize_list[1])

        hull_our = ConvexHull(xy_our_hull,qhull_options = qhull_options )

        for simplex in hull_our.simplices:
            line_our_99 = ax.plot(xy_our_hull[simplex, 0], xy_our_hull[simplex, 1], 'x-',markerfacecolor='none', color=colors_list[1],
                            alpha = 0.5, markersize=10,lw=lw_sizelist[1])

        xy_our_hull = our_all_points_hull_act[:, [0, index]]
        ax.plot(xy_our_hull[:,0], xy_our_hull[:,1],          'x', markerfacecolor='none', color=colors_list[2],
                    alpha=0.8,  markersize=markersize_list[2])

        hull_our = ConvexHull(xy_our_hull,qhull_options = qhull_options )
        for simplex in hull_our.simplices:
            line_our_ac = ax.plot(xy_our_hull[simplex, 0], xy_our_hull[simplex, 1], 'x-',markerfacecolor='none', color=colors_list[2],
                            alpha = 0.5, markersize=10,lw=lw_sizelist[2])


    else:
        # line_EFMs = points_EFMs
        # line_EFVs = points_EFVs
        # line_FBAem = points_FBAem
        line_our = points_our

    ax.set_ylabel(yield_rea_ids_name[index-1]+' / Glc',fontsize = 18)

ax.set_xlabel('Yield Biomass/Glucose',fontsize = 18)
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
# tspan = np.linspace(tStart, tStop, (tStop - tStart) / tStep)
tspan = np.linspace(tStart, tStop, int(((tStop - tStart) / tStep) + 1))

# matrix Z: reactions x pathways; and Smz metabolites x pathways , S metabolites x reactions : Smz = Sm @ Z
final_index = hull_active_index
Smz = yield_normalized_df_hull_.values[:, final_index]
Smz[1, :] = Smz[1, :] * coefficient
path_for = np.array([0, 0, 0, -1, 0, 0, 0])  # Note !!! this pathway is nessary  to simulate the for experimentdata
# Smz = np.column_stack((Smz,path_for))
Smz = np.insert(Smz, Smz.shape[1], values=path_for, axis=1)
# metabolites and pathways number
(n_mets, n_path) = Smz.shape
# np.savetxt('Case3_iML1515/iML1515_Smz.txt', Smz, delimiter=',') # note： check the result with saved file
Smz = np.loadtxt('Case3_iML1515/iML1515_Smz.txt', delimiter=',')
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

# %%
# minimums = []
# for i in range(0,7):
#     CB_model.weight_of_method_t_f = i
#     # cb_model_1.weights_of_mets = weights_of_mets/1000
#     cb_models = [CB_model, ]
#
#     experiment_data_dfs = [experiment_data_df, ]
#     para_to_fit = {'kmax': np.array([[0], [1], [2], [3],[4],[5], [6], [7],]),
#                    'K': np.array([])}
#
#     minimum = Cybernetic_Functions.parameters_fitting(cb_models, experiment_data_dfs,
#                                                         para_to_fit, tspan, draw=True,
#                                                         full_output=False,
#                                                         options={'xtol': 5, 'ftol': 0.01, 'maxiter': 100,
#                                                                  'disp': True})
#
#
#
#     minimums.append(minimum)
#
# for minimum in minimums:
#     CB_models = Cybernetic_Functions.update_paras_func(minimum.x, [CB_model], para_to_fit,
#                                                        retern_model_or_paras='model')
#     sol3 = Cybernetic_Functions.cb_model_simulate(CB_models[0], tspan, draw=True)


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
