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

import cobra
import numpy as np

import ConvexHull_yield

os.chdir('../../ComplementaryData/Branch_work')

yeast = cobra.io.read_sbml_model('yeastGEM.xml')
# s_cer = cobra.io.read_sbml_model('Saccharomyces_cerevisiae.xml')
# p_sti = cobra.io.read_sbml_model('yeastGEM.xml')
# k_mar = cobra.io.read_sbml_model('Kluyveromyces_marxianus.xml')

# r_1714	D-glucose exchange	D-glucose[e] <=>
# r_1718	D-xylose exchange	D-xylose[e] =>
# r_1715	D-mannose exchange	D-mannose[e] =>
# r_1710	D-galactose exchange	D-galactose[e] =>
# r_1761	ethanol exchange	ethanol[e] =>
# r_2111	growth	biomass[c] =>

# % Constrain the phosphotransferase system

solution = yeast.optimize()
print(solution)
# # %% < glc >
# model = yeast.copy()
# production_rea_ids_x = ['r_2111', ]
# production_rea_ids_y = ['r_1761', ]
# carbon_source_rea_id = 'r_1714'
#
# steps = 10
# carbon_uptake_direction = -1
# model.reactions.get_by_id(carbon_source_rea_id).bounds = (-10, -0.1)
#
# # all modes
# yield_normalized_df, fluxes_all, hull_index_all = GEM2pathways.get_yield_space_multi(model, production_rea_ids_x,
#                                                                                      production_rea_ids_y,
#                                                                                      carbon_source_rea_id,
#                                                                                      steps=steps,
#                                                                                      carbon_uptake_direction=carbon_uptake_direction,
#                                                                                      draw=True)
# yield_normalized_df_hull = yield_normalized_df[yield_normalized_df.columns[hull_index_all]]
#
# yield_normalized_df_hull_ = yield_normalized_df_hull
# our_all_points = yield_normalized_df_hull_.values.T
# our_all_points = our_all_points[abs(our_all_points[:, 0]) > 1e-10, :]
# our_all_points = our_all_points[:, 1:]  # exclude the carbon colume
#
# # NOTE: according to the our_estimated_datas, set the biomacc_mas_coefficient = 5.0 it means 1 in model , 1 *5.0 in experiment data
# coefficient = 1
# molar_mass = np.array([44.03 * 26.04, 46, 180, 150, 180, 180])  # unit conversion : 0.01g/l = 1000*0.01/900 mM(mmol/l)
#
# our_all_points[:, 0] = our_all_points[:, 0] * coefficient
#
# print('Our method:', our_all_points.shape)
#
# # %% <MYA >
#
# print('\n---------- Loading experiment data ... ---------- ')
#
# experiment_data_df_trimed_values = np.array([[0.105*180/900,0.403*180/46]])
#
# experiment_datas = []  # TODO experiment data!!!
# for i in range(0, experiment_data_df_trimed_values.shape[0]):
#     experiment_data = list(experiment_data_df_trimed_values[i, :])
#     # experiment_data[0] = ''
#     experiment_datas.append(experiment_data)
#
# qhull_options = 'QJ Qx A0.9999999'
# cutoff_persent = 1
# our_indexes, our_weights, our_estimated_datas, our_in_hulls = ConvexHull_yield.pipeline_mya(our_all_points,
#                                                                                             experiment_datas,
#                                                                                             qhull_options=qhull_options,
#                                                                                             method=1,
#                                                                                             cutoff_persent=cutoff_persent)
# # hull_active_index = 	 [0, 4, 3]
# # our_hull = ConvexHull(our_all_points[:, :], qhull_options=qhull_options)
# # ConvexHull_yield.point_in_hull(experiment_datas[-1][1:], our_hull, tolerance=1e-12)
#
# # hull_cutoff_index_99 = [0, 2, 3, 9, 10, 11, 16, 17, 18, 19, 20, 23, 26, 30, 33, 34, 35]
# hull_active_index = our_indexes[-1]  # list(set(our_indexes[-1][-1]) | set(our_indexes[-1][-2]) | set(our_indexes[-1][-3]))
#
# our_all_points_hull_act = our_all_points[hull_active_index, :]
#
# m1_glc = our_all_points_hull_act
#
#
# # %% <MYA >
#
# print('\n---------- Loading experiment data ... ---------- ')
#
# experiment_data_df_trimed_values = np.array([[0.129*180/900,0.40*180/46]])
#
# experiment_datas = []  # TODO experiment data!!!
# for i in range(0, experiment_data_df_trimed_values.shape[0]):
#     experiment_data = list(experiment_data_df_trimed_values[i, :])
#     # experiment_data[0] = ''
#     experiment_datas.append(experiment_data)
#
# qhull_options = 'QJ Qx A0.9999999'
# cutoff_persent = 1
# our_indexes, our_weights, our_estimated_datas, our_in_hulls = ConvexHull_yield.pipeline_mya(our_all_points,
#                                                                                             experiment_datas,
#                                                                                             qhull_options=qhull_options,
#                                                                                             method=1,
#                                                                                             cutoff_persent=cutoff_persent)
# hull_active_index = our_indexes[-1]  # list(set(our_indexes[-1][-1]) | set(our_indexes[-1][-2]) | set(our_indexes[-1][-3]))
#
# our_all_points_hull_act = our_all_points[hull_active_index, :]
#
# m2_glc = our_all_points_hull_act
#
#
#
# # %% <MYA >
#
# print('\n---------- Loading experiment data ... ---------- ')
#
# experiment_data_df_trimed_values = np.array([[0.140*180/900,0.421*180/46]])
#
# experiment_datas = []  # TODO experiment data!!!
# for i in range(0, experiment_data_df_trimed_values.shape[0]):
#     experiment_data = list(experiment_data_df_trimed_values[i, :])
#     # experiment_data[0] = ''
#     experiment_datas.append(experiment_data)
#
# qhull_options = 'QJ Qx A0.9999999'
# cutoff_persent = 1
# our_indexes, our_weights, our_estimated_datas, our_in_hulls = ConvexHull_yield.pipeline_mya(our_all_points,
#                                                                                             experiment_datas,
#                                                                                             qhull_options=qhull_options,
#                                                                                             method=1,
#                                                                                             cutoff_persent=cutoff_persent)
# hull_active_index = our_indexes[-1]  # list(set(our_indexes[-1][-1]) | set(our_indexes[-1][-2]) | set(our_indexes[-1][-3]))
#
# our_all_points_hull_act = our_all_points[hull_active_index, :]
#
# m3_glc = our_all_points_hull_act

# %% all

# our_all_points_z = []
# for carbon_i in ['r_1714', 'r_1718', 'r_1715', 'r_1710']:
#     model = yeast.copy()
#     # model.solver = 'cplex'
#     model.reactions.get_by_id('r_1714').bounds = (0, 1000)
#     production_rea_ids_x = ['r_2111', ]
#     production_rea_ids_y = ['r_1761', ]
#     carbon_source_rea_id = carbon_i
#
#     steps = 10
#     carbon_uptake_direction = -1
#     model.reactions.get_by_id(carbon_source_rea_id).bounds = (-10, -0.5)
#
#     # all modes
#     yield_normalized_df, fluxes_all, hull_index_all = GEM2pathways.get_yield_space_multi(model, production_rea_ids_x,
#                                                                                          production_rea_ids_y,
#                                                                                          carbon_source_rea_id,
#                                                                                          steps=steps,
#                                                                                          carbon_uptake_direction=carbon_uptake_direction,
#                                                                                          draw=True)
#     yield_normalized_df_hull = yield_normalized_df[yield_normalized_df.columns[hull_index_all]]
#
#     yield_normalized_df_hull_ = yield_normalized_df_hull
#     our_all_points = yield_normalized_df_hull_.values.T
#     our_all_points = our_all_points[abs(our_all_points[:, 0]) > 1e-10, :]
#     our_all_points = our_all_points[:, 1:]  # exclude the carbon colume
#
#     coefficient = 1
#
#     our_all_points[:, 0] = our_all_points[:, 0] * coefficient
#     our_all_points_z.append(our_all_points)
#
#     print('Our method:', our_all_points.shape)
#
# file = open('our_all_points_z', 'wb')
# pickle.dump(our_all_points_z, file)
# file.close()

our_all_points_z = pickle.load(open('our_all_points_z', 'rb'))

etoh_Ys = [0.403, 0.409, 0.421]
biomass_Ys = [0.105, 0.129, 0.140]
Bio_mass_s = []

for i in range(0, 3):
    model = yeast.copy()
    etoh_Y = etoh_Ys[i] * 180 / 46
    biomass_Y = biomass_Ys[i] * 180
    model.reactions.get_by_id('r_1761').bounds = (etoh_Y, etoh_Y)
    solution = model.optimize()
    Bio_mass = biomass_Y / solution.objective_value
    Bio_mass_s.append(Bio_mass)

Bio_mass = np.mean(Bio_mass_s)

m1_datas = [[0.105 * 180 / Bio_mass, 0.403 * 180 / 46],
            [0, 0],
            [0.203 * 180 / Bio_mass, 0.364 * 180 / 46],
            [0.105 * 150 / Bio_mass, 0.340 * 180 / 46]]

m2_datas = [[0.129 * 180 / Bio_mass, 0.409 * 180 / 46],
            [0.147 * 150 / Bio_mass, 0.404 * 150 / 46],
            [0.156 * 180 / Bio_mass, 0.372 * 180 / 46],
            [0.146 * 180 / Bio_mass, 0.366 * 180 / 46]]

m3_datas = [[0.140 * 180 / Bio_mass, 0.421 * 180 / 46],
            [0.241 * 150 / Bio_mass, 0.251 * 150 / 46],
            [0.174 * 180 / Bio_mass, 0.426 * 180 / 46],
            [0.126 * 180 / Bio_mass, 0.412 * 180 / 46]]

datas = [m1_datas, m2_datas, m3_datas]
# %%

all_points_hull_acts = []
for model_i in range(0, 3):
    for carbon_i in range(0, 4):
        our_all_points = our_all_points_z[carbon_i]
        experiment_data_df_trimed_values = np.array([datas[model_i][carbon_i]])
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
        # hull_active_index = 	 [0, 4, 3]
        # our_hull = ConvexHull(our_all_points[:, :], qhull_options=qhull_options)
        # ConvexHull_yield.point_in_hull(experiment_datas[-1][1:], our_hull, tolerance=1e-12)

        hull_active_index = our_indexes[
            -1]  # list(set(our_indexes[-1][-1]) | set(our_indexes[-1][-2]) | set(our_indexes[-1][-3]))

        our_all_points_hull_act = our_all_points[hull_active_index, :]
        all_points_hull_acts.append(our_all_points_hull_act)

# %%
for i in range(0, 4):
    temp_z = np.zeros((6, all_points_hull_acts[i].shape[0]))
    temp_z[i, :] = -1
    temp_z[[4, 5], :] = all_points_hull_acts[i].T
    if i == 0:
        model_1_smz = temp_z
    elif i == 1:
        continue
    else:
        model_1_smz = np.append(model_1_smz, temp_z, axis=1)

# model_1_smz = np.around(model_1_smz, decimals=5, )
print(model_1_smz)
print(model_1_smz.shape)
for i in range(4, 8):
    temp_z = np.zeros((6, all_points_hull_acts[i].shape[0]))
    temp_z[i - 4, :] = -1
    temp_z[[4, 5], :] = all_points_hull_acts[i].T
    if i - 4 == 0:
        model_2_smz = temp_z
    else:
        model_2_smz = np.append(model_2_smz, temp_z, axis=1)
# model_2_smz = np.around(model_2_smz, decimals=5, )
print(model_2_smz)
print(model_2_smz.shape)

for i in range(8, 12):
    temp_z = np.zeros((6, all_points_hull_acts[i].shape[0]))
    temp_z[i - 8, :] = -1
    temp_z[[4, 5], :] = all_points_hull_acts[i].T
    if i - 8 == 0:
        model_3_smz = temp_z
    else:
        model_3_smz = np.append(model_3_smz, temp_z, axis=1)
# model_3_smz = np.around(model_3_smz, decimals=5, )
print(model_3_smz)
print(model_3_smz.shape)

model_smz = [model_1_smz, model_2_smz, model_3_smz]
file = open('model_smz', 'wb')
pickle.dump(model_smz, file)
file.close()
