#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 7/8/20

"""Bhyd_get_pathways.py
:description : script
:param : 
:returns: 
:rtype: 
"""
import os

import cobra
import numpy as np
import pandas as pd

import ConvexHull_yield
import GEM2pathways

os.chdir('../../Data/three_species/')

# %% <load model>
print('\n---------- Loading model ... ---------- ')
Bhyd_model = cobra.io.load_json_model('GEM_from_templates_B_hyd_refined.json')
initial_model = Bhyd_model.copy()
model = initial_model.copy()

# %% <load experiment data>
experiment_df = pd.read_excel('experiment_data/experiment_data_trimmed.xlsx', sheet_name='BH_14', index_col=0, header=0)
data_cloum_name = ['fru', 'R.i', 'F.p', 'B.h', 'ac', 'for', 'but', 'lac']

### yield of experiment data
experiment_df_trimmed = experiment_df[['fru', 'B.h', 'ac']]
experiment_df_trimmed_values = experiment_df_trimmed.values[:, :] - experiment_df_trimmed.values[0, :]  # - time0 value
experiment_df_trimmed_values = experiment_df_trimmed_values[1:, :].T  # remove time 0
experiment_df_trimmed_values = experiment_df_trimmed_values[:, :] / abs(
    experiment_df_trimmed_values[0, :])  # / carbon uptake
experiment_df_trimmed_values = experiment_df_trimmed_values.T[:, 1:]  # remove carbon (all are -1 )

final_yield_row = experiment_df_trimmed_values[-5:, :]  # stable yield range
final_yield_row = final_yield_row.mean(axis=0)

### check experiment result
model.reactions.get_by_id('EX_o2_e').bounds = (0.0, 1000.0)
model.reactions.get_by_id('EX_glc__D_e').bounds = (0, 1000)
model.reactions.get_by_id('EX_fru_e').bounds = (-16, 1000)
model.reactions.get_by_id('EX_ac_e').bounds = (20, 1000)
model.reactions.get_by_id('EX_for_e').bounds = (-40, 1000)
model.reactions.get_by_id('EX_but_e').bounds = (0, 1000)
model.objective = 'BIOMASS'
solution = model.optimize()  # check growth
print(solution)  # 1.042


###  model reset
model.reactions.get_by_id('EX_fru_e').bounds = (-10, 1000)
model.reactions.get_by_id('EX_ac_e').bounds = (0, 1000)
model.reactions.get_by_id('EX_for_e').bounds = (0, 1000)
model.reactions.get_by_id('EX_but_e').bounds = (0, 1000)
model.objective = 'BIOMASS'
solution = model.optimize()  # check growth
print(solution)  # 0.721

### set GEM2pathways prams
production_rea_ids_x = ['BIOMASS', ]
production_rea_ids_y = ['EX_ac_e', ]
carbon_source_rea_id = 'EX_fru_e'

steps = 10
carbon_uptake_direction = -1
model.reactions.get_by_id(carbon_source_rea_id).bounds = (-10, -3)  # at least uptake 3 carbon_source, max 10

yield_row = final_yield_row[1:]
constrains = {'x': [0.1],
              'y_no_obj': [yield_row * 0.84, yield_row * 1.1],
              'y_obj': [yield_row * 0.2, yield_row * 2]}
## all modes
yield_normalized_df, fluxes_all, hull_index_all = GEM2pathways.get_yield_space_multi(model, production_rea_ids_x,
                                                                                     production_rea_ids_y,
                                                                                     carbon_source_rea_id,
                                                                                     steps=steps,
                                                                                     carbon_uptake_direction=carbon_uptake_direction,
                                                                                     draw=True, constrains=constrains)
yield_normalized_df_hull = yield_normalized_df[yield_normalized_df.columns[hull_index_all]]
#
# yield_normalized_df_hull.to_csv(model.id+'_yield_normalized_df_hull_constrain.csv', sep=',',index=True)

yield_normalized_df_hull_ = pd.read_csv(model.id + '_yield_normalized_df_hull_constrain.csv',
                                        index_col=0)

our_all_points = yield_normalized_df_hull_.values.T
our_all_points = our_all_points[abs(our_all_points[:, 0]) > 1e-10, :]
our_all_points = our_all_points[:, 1:]  # exclude the carbon colume

# %% <biomass coefficient>  convert biomass unites, g/mol/ml --> n
model = initial_model.copy()
model.reactions.get_by_id('EX_o2_e').bounds = (0.0, 1000.0)
model.reactions.get_by_id('EX_glc__D_e').bounds = (0, 1000)
model.reactions.get_by_id('EX_fru_e').bounds = (-16, 1000)
model.reactions.get_by_id('EX_ac_e').bounds = (20, 1000)
model.reactions.get_by_id('EX_for_e').bounds = (-40, 1000)
model.reactions.get_by_id('EX_but_e').bounds = (0, 1000)
model.objective = 'BIOMASS'
solution = model.optimize()  # check growth
print(solution)  # 1.042
# in model biomass yield = 0.863g/10 from 10 fru and 10 ac; in experiment data final_yield_row[0]: 1.27

coefficient = final_yield_row[0] / (solution.objective_value / 16)
print(coefficient == 39.69474611037711)
coefficient = 39.69474611037711

our_all_points[:, 0] = our_all_points[:, 0] * coefficient

print('Our method points shape:', our_all_points.shape)

# additional process: because of no only EX_fru_e as only carbon, we guss the 0.79 from fur, rest from ac.
temp_dir = 'initial_data/templates/'
iML1515 = cobra.io.load_json_model(temp_dir + 'iML1515_standlized.json')
iML1515.reactions.get_by_id('EX_ac_e').bounds = (0, 1000)
iML1515.reactions.get_by_id('EX_glc__D_e').bounds = (-16, 1000)
f1 = iML1515.optimize()
iML1515.reactions.get_by_id('EX_for_e').bounds = (-40, 1000)
iML1515.reactions.get_by_id('EX_glc__D_e').bounds = (-16, 1000)
f2 = iML1515.optimize()
rate_glc_ac = f1.objective_value / f2.objective_value
rate_glc_ac = 0.8627442009424916

pathway_of_ac = final_yield_row * (1 - rate_glc_ac)
# for next step!!  array([0.35480667, 0.19953176])

experiment_df_trimmed_values = experiment_df_trimmed_values * rate_glc_ac
experiment_datas = []  # TODO experiment data!!!
for i in range(0, experiment_df_trimmed_values.shape[0]):
    experiment_data = list(experiment_df_trimmed_values[i, :])
    # experiment_data[0] = ''
    experiment_datas.append(experiment_data)
experiment_datas.append(list(final_yield_row * rate_glc_ac))

# %% <MYA >
print('\n---------- Loading experiment data ... ---------- ')

qhull_options = 'QJ Qx A0.9999999'
cutoff_persent = 1
our_indexes, our_weights, our_estimated_datas, our_in_hulls = ConvexHull_yield.pipeline_mya(our_all_points,
                                                                                            experiment_datas,
                                                                                            qhull_options=qhull_options,
                                                                                            method=1,
                                                                                            cutoff_persent=cutoff_persent)

# our_hull = ConvexHull(our_all_points[:, 1:], qhull_options=qhull_options)
# ConvexHull_yield.point_in_hull(experiment_datas[-1][1:], our_hull, tolerance=1e-12)

hull_cutoff_index_99 = our_indexes[1]
hull_cutoff_index_99 = [0, 1, 2, 7, 8, 9, 10]
hull_active_index = our_indexes[-1][-1]
hull_active_index = [3, 1, 4]
# weights =	 [0.341917 0.       0.352333 0.       0.30575 ]


# our_all_points_hull_1 = our_all_points[our_indexes[0], :]
# our_all_points_hull_99 = our_all_points[hull_cutoff_index_99, :]
# our_all_points_hull_act = our_all_points[hull_active_index, :]

final_index = hull_active_index
Smz = yield_normalized_df_hull_.values[:, final_index]
Smz[1, :] = Smz[1, :] * coefficient
data_cloum_name = ['fru', 'R.i', 'F.p', 'B.h', 'ac', 'for', 'but', 'lac']

other_biomass = np.array([0] * len(final_index))  # the row for 'R.i', 'F.p', 'B.h', 'ac', 'for', 'but',
Smz = np.insert(Smz, 1, values=other_biomass, axis=0)  # 'R.i',
Smz = np.insert(Smz, 2, values=other_biomass, axis=0)  # 'F.p',
Smz = np.insert(Smz, 5, values=other_biomass, axis=0)  # 'for',
Smz = np.insert(Smz, 6, values=other_biomass, axis=0)  # 'but',

pathway_of_ac_path = np.array([-0.001, 0, 0, pathway_of_ac[0], pathway_of_ac[1], -1, 0])  # for
Smz = np.insert(Smz, Smz.shape[1], values=pathway_of_ac_path, axis=1)
# Note !!! this pathway is nessary  to simulate the for experimentdata
# Smz = np.column_stack((Smz,path_for))
# metabolites and pathways number
(n_mets, n_path) = Smz.shape

np.savetxt(model.id + '_Smz_constrain.csv', Smz, delimiter=',')  # note： check the result with saved file
# Smz_ = np.genfromtxt(model.id + '_Smz_constrain.csv', delimiter=',')
print(Smz_ == Smz)
