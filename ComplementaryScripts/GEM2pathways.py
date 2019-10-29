#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 9/24/19

"""GEM2pathways.py
:description : script to caculate the pathways from GEMs by Yield optimal
:param : 
:returns: 
:rtype: 
"""

import cobra
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# %%
def get_yield_opt(model2, production_rea_id, carbon_source_rea_id, max_or_min='max', carbon_uptake_direction=-1):
    '''
    Get the max or min Yield , Y = production / carbon

    Math :
    𝑌max⁡ = 𝑌𝑡𝑒𝑚𝑝
    Wℎ𝑒𝑛:max⁡(𝑟𝑝 −𝑌𝑡𝑒𝑚𝑝∙𝑟𝑠) = 0
    𝑌𝑡𝑒𝑚𝑝 ≥⁡ 0
    𝑆∙r = 0
    rlb ≤ 𝑟 ≤ rub
    𝑟𝑠>0

    param :
    :param model2:                  A GEM; cobrapy model
    :param production_rea_id:       production reaction id; str
    :param carbon_source_rea_id:    carbon source reaction id, the bounds can't be 0!;  int
    :param max_or_min:              'max' or 'min' .opt direction;  str
    :param carbon_uptake_direction: 1 or -1 the flux direction when uptake the source;  int
    :return:
    yield_opt,                      opted value;    flot
    fluxes_opt                      fluexs distribuation;   dataframe

    Examples:
    yield_value,fluxes = get_yield_opt(e_coli_core,'EX_ac_e','EX_glc__D_e',max_or_min = 'max',carbon_uptake_direction = -1)
    # print(yield_value)
    '''

    model = model2.copy()  # not change initial model
    model.objective.direction = max_or_min  # max or min

    # first round, yield >0, so test from 0, (when yield_temp = 0 , the method is same as FBA)
    yield_temp = 0
    model.objective = {model.reactions.get_by_id(production_rea_id): 1,
                       model.reactions.get_by_id(carbon_source_rea_id): carbon_uptake_direction * (-yield_temp)}
    selution = model.optimize()
    # print(selution.objective_value)

    # Next round,  when max⁡(𝑟𝑝 −𝑌𝑡𝑒𝑚𝑝∙𝑟𝑠) = 0  𝑌max⁡ = 𝑌𝑡𝑒𝑚𝑝, 1e-6~ 0, a cut off
    while selution.objective_value > 1e-6:
        yield_temp = selution.fluxes[production_rea_id] / (carbon_uptake_direction * selution.fluxes[
            carbon_source_rea_id])  # get fluxes, then get new yield_temp

        model.objective = {model.reactions.get_by_id(production_rea_id): 1,
                           model.reactions.get_by_id(carbon_source_rea_id): carbon_uptake_direction * (
                               -yield_temp)}  # add constraints
        selution = model.optimize()

    # get yield_opt and fluxes
    fluxes_opt = selution.fluxes
    yield_opt = fluxes_opt[production_rea_id] / (carbon_uptake_direction * fluxes_opt[carbon_source_rea_id])

    return yield_opt, fluxes_opt


# yield_value,fluxes = get_yield_opt(e_coli_core,'EX_ac_e','EX_glc__D_e',max_or_min = 'max',carbon_uptake_direction = -1)
# print(yield_value)

# %% (TODO parallel!!) (tried, but Filed (computer)

def _get_yield_opt(tasks):
    return get_yield_opt(*tasks)


def get_yield_space_2d(model2, production_rea_ids_2d, carbon_source_rea_ids_2d, steps, carbon_uptake_direction=-1,
                       draw=False):
    '''
    Get the Yield space of x =p1/s1 y = p2/s2

    :param model2:                  A GEM
    :param production_rea_ids:       production reaction ids;   list
    :param carbon_source_rea_id:    carbon source reaction ids, the bounds can't be 0!; list
    :param steps:                   steps ;    int
    :param carbon_uptake_direction: 1 or -1 the flux direction when uptake the source;  int
    :param draw:                    True or False get fig or not ;  bool
    :return:
    fluxes_2d,                     pd.dataframe of fluxes

    Examples:
    fluxes_2d = get_yield_space_2d(e_coli_core, production_rea_ids_2d = ['BIOMASS_Ecoli_core_w_GAM','EX_ac_e']
                                                ,carbon_source_rea_ids_2d = ['EX_glc__D_e','EX_glc__D_e'],steps = 20 ,
                                                carbon_uptake_direction = -1,draw = True)
    '''
    # <Step1 check input data and find the max and min value of p1 yield >

    model = model2.copy()

    if type(carbon_source_rea_ids_2d) == str:  # test input type and size
        carbon_source_rea_ids_2d_normalized = [carbon_source_rea_ids_2d] * 2
    elif type(carbon_source_rea_ids_2d) == list:
        if len(carbon_source_rea_ids_2d) == 1:
            carbon_source_rea_ids_2d_normalized = carbon_source_rea_ids_2d * 2
        else:
            carbon_source_rea_ids_2d_normalized = carbon_source_rea_ids_2d

    p1_rea_id = production_rea_ids_2d[0]
    p2_rea_id = production_rea_ids_2d[1]
    s1_rea_id = carbon_source_rea_ids_2d_normalized[0]
    s2_rea_id = carbon_source_rea_ids_2d_normalized[1]

    model.objective = p1_rea_id
    # solution = model.optimize()
    # print('model.optimize():\t', solution.objective_value)

    p1_yield_value_max, _ = get_yield_opt(model, p1_rea_id, s1_rea_id, max_or_min='max',
                                          carbon_uptake_direction=carbon_uptake_direction)
    p1_yield_value_min, _ = get_yield_opt(model, p1_rea_id, s1_rea_id, max_or_min='min',
                                          carbon_uptake_direction=carbon_uptake_direction)
    print(p1_rea_id, 'yield\t:\t range from \t', p1_yield_value_min, '\tto\t', p1_yield_value_max)

    # < Step2 get yield dataFrame (matrix/ points for MYA/Next function)>

    fluxes_2d = pd.DataFrame()

    for step in np.arange(0, steps + 1):
        persent_i = step / steps

        yield_k = p1_yield_value_min + (p1_yield_value_max - p1_yield_value_min) * persent_i  # split yield_p1

        # set p1 yield constraints
        model = model2.copy()
        yield_point_flux = model.problem.Constraint(
            model.reactions.get_by_id(
                p1_rea_id).flux_expression - carbon_uptake_direction * yield_k * model.reactions.get_by_id(
                s1_rea_id).flux_expression,
            lb=0,
            ub=0)
        model.add_cons_vars(yield_point_flux)

        # Get all cores     multiprocessing but filed...
        # cores =2    # multiprocessing.cpu_count()
        #
        # pool = multiprocessing.Pool(processes=cores,)
        #
        # tasks = [(model,p2_rea_id,s2_rea_id,'max',carbon_uptake_direction),
        #          (model,p2_rea_id,s2_rea_id,'min',carbon_uptake_direction)]
        #
        # fluxes_list = []
        # for i in pool.imap(_get_yield_opt,tasks):
        #     fluxes = i[1]
        #     fluxes_list = fluxes_list+[fluxes]
        #
        # pool.close()
        # pool.join()
        #
        # fluxes_2d[p2_rea_id + '_max_' + str(persent_i) + '_x_' + p1_rea_id] = fluxes_list[0]
        # fluxes_2d[p2_rea_id + '_min_' + str(persent_i) + '_x_' + p1_rea_id] = fluxes_list[1]

        _, fluxes_max = get_yield_opt(model, p2_rea_id, s2_rea_id, max_or_min='max',
                                      carbon_uptake_direction=carbon_uptake_direction)
        _, fluxes_min = get_yield_opt(model, p2_rea_id, s2_rea_id, max_or_min='min',
                                      carbon_uptake_direction=carbon_uptake_direction)

        fluxes_2d[p2_rea_id + '_max_' + str(persent_i) + '_x_' + p1_rea_id] = fluxes_max
        fluxes_2d[p2_rea_id + '_min_' + str(persent_i) + '_x_' + p1_rea_id] = fluxes_min

    # < Step3 draw >
    if draw:
        fig = plt.figure()
        ax = fig.add_subplot(111)

        temp_columns_max = [column for column in fluxes_2d.columns if ('max' in column)]
        temp_columns_min = [column for column in fluxes_2d.columns if ('min' in column)]

        x = fluxes_2d.loc[p1_rea_id, temp_columns_max].values / abs(fluxes_2d.loc[s1_rea_id, temp_columns_max].values)
        y_max = fluxes_2d.loc[p2_rea_id, temp_columns_max].values / abs(
            fluxes_2d.loc[s2_rea_id, temp_columns_max].values)
        y_min = fluxes_2d.loc[p2_rea_id, temp_columns_min].values / abs(
            fluxes_2d.loc[s2_rea_id, temp_columns_min].values)

        ax.plot([x[0], x[0]], [y_max[0], y_min[0]], 'x-', color='tab:blue', alpha=0.8)
        ax.plot(x, y_max, 'x-', color='tab:blue', alpha=0.8)
        ax.plot(x, y_min, 'x-', color='tab:blue', alpha=0.8)

        ax.set_ylabel(p1_rea_id + '/' + s1_rea_id)
        ax.set_xlabel(p1_rea_id + '/' + s2_rea_id)
        fig.show()

    return fluxes_2d


# fluxes_2d = get_yield_space_2d(e_coli_core, production_rea_ids_2d = ['BIOMASS_Ecoli_core_w_GAM','EX_ac_e'],
#                                          carbon_source_rea_ids_2d = ['EX_glc__D_e','EX_glc__D_e'],steps = 20 ,
#                                          carbon_uptake_direction = -1,draw = True)

# %% TODO: real multi-dimension?
def get_yield_space_multi(model2, production_rea_ids_x, production_rea_ids_y, carbon_source_rea_id,
                          steps, carbon_uptake_direction=-1, draw=False):
    '''
    Get the Yield space of x1 =p1/s y1 = p2/s2 ; x2 =p1/s y2 = p2/s
    only for single carbon_source

    :param model2:                  A GEM
    :param production_rea_ids:       production reaction ids;   list
    :param carbon_source_rea_id:    carbon source reaction ids, the bounds can't be 0!; list
    :param steps:                   steps ;    int
    :param carbon_uptake_direction: 1 or -1 the flux direction when uptake the source;  int
    :param draw:                    True or False get fig or not ;  bool
    :return:
    fluxes_2d,                     dataframe of fluxes

    Examples:
    fluxes_2d = get_yield_space_2d(e_coli_core, production_rea_ids_2d = ['BIOMASS_Ecoli_core_w_GAM','EX_ac_e']
                                                ,carbon_source_rea_ids_2d = ['EX_glc__D_e','EX_glc__D_e'],steps = 20 ,
                                                carbon_uptake_direction = -1,draw = True)
    '''

    # <Step1 get production_rea_ids_2d>

    production_rea_ids_2ds = []
    for x_rea_i in production_rea_ids_x:
        for y_rea_i in production_rea_ids_y:
            production_rea_ids_2ds.append([x_rea_i, y_rea_i])

    # < Step2 get yield dataFrame(matrix/ points for MYA/Next function)>
    fluxes_all = pd.DataFrame()
    fluxes_all_list = []
    for production_rea_ids_2d in production_rea_ids_2ds:
        model = model2.copy()
        fluxes_2d = get_yield_space_2d(model, production_rea_ids_2d, carbon_source_rea_id, steps=steps,
                                       carbon_uptake_direction=carbon_uptake_direction, draw=False)
        fluxes_all = pd.concat([fluxes_all, fluxes_2d], axis=1)
        fluxes_all_list = fluxes_all_list + [fluxes_2d]

    index_list = [carbon_source_rea_id] + production_rea_ids_x
    for i in production_rea_ids_y:
        if i not in index_list:
            index_list.append(i)

    yield_df = fluxes_all.loc[index_list, :]  # select reactions
    yield_normalized_df = yield_df.copy()
    # yield_normalized_df = yield_normalized_df.T.drop_duplicates().T

    yield_df_values  = yield_normalized_df.values      #normalize and sort
    yield_df_values = yield_df_values/abs(yield_df_values[0,:])

    yield_normalized_df.loc[:,:] = yield_df_values
    # yield_normalized_df = yield_normalized_df.sort_values(by = index_list[1:],axis = 1,)

    # < Step3 draw >
    for x_i in production_rea_ids_x:

        yield_normalized_df_temp = yield_normalized_df.loc[:,
                                   [column for column in yield_normalized_df.columns if ('_x_' + x_i in column)]]

        if draw:
            def trim_axs(axs, N):
                """little helper to massage the axs list to have correct length..."""
                axs = axs.flat
                for ax in axs[N:]:
                    ax.remove()
                return axs[:N]

            cols = 3
            rows_i = 1
            if len(production_rea_ids_y) % cols == 0:
                rows_i = 0
            rows = (len(production_rea_ids_y)) // cols + rows_i
            figsize = (10, 8)
            fig, axs = plt.subplots(cols, rows, figsize=figsize)
            axs = trim_axs(axs, len(production_rea_ids_y))

            for ax, exmet_reaction in zip(axs, production_rea_ids_y):
                temp_columns_max = [column for column in yield_normalized_df_temp.columns if
                                    (exmet_reaction in column) and ('max' in column)]
                temp_columns_min = [column for column in yield_normalized_df_temp.columns if
                                    (exmet_reaction in column) and ('min' in column)]
                x = yield_normalized_df_temp.loc[x_i, temp_columns_max].values
                y_max = yield_normalized_df_temp.loc[exmet_reaction, temp_columns_max]
                y_min = yield_normalized_df_temp.loc[exmet_reaction, temp_columns_min]

                # lines:
                ax.plot([x[0], x[0]], [y_max[0], y_min[0]], 'x-', color='tab:blue', alpha=0.5)
                ax.plot([x[-1], x[-1]], [y_max[-1], y_min[-1]], 'x-', color='tab:blue', alpha=0.5)
                ax.plot(x, y_max, 'x-', color='tab:blue', alpha=0.5)
                ax.plot(x, y_min, 'x-', color='tab:blue', alpha=0.5)

                # points
                x_points = yield_normalized_df.loc[x_i, :].values
                y_points = yield_normalized_df.loc[exmet_reaction, :].values
                ax.plot(x_points, y_points, 'x', color='tab:blue', alpha=0.5)

                ax.set_ylabel(exmet_reaction + '/' + carbon_source_rea_id)
            ax.set_xlabel(x_i + '/' + carbon_source_rea_id)

            fig.show()

    return yield_normalized_df


# %%
if __name__ == '__main__':
    # load model

    e_coli_core = cobra.io.read_sbml_model('../ComplementaryData/Case2_1_ecoli_core/e_coli_core.xml')
    e_coli_core.reactions.get_by_id('EX_o2_e').bounds = (0.0, 1000.0)
    e_coli_core.reactions.get_by_id('EX_glc__D_e').bounds = (-10.0, -0.01)
    steps = 20
    model = e_coli_core
    draw = True

    # %% <get_yield_opt> model, biomass_rea_id,carbon_source_rea_ids,yield_rea_ids,step_of_biomass
    print('----------   get_yield_opt   ----------\n')
    yield_opt, fluxes = get_yield_opt(model, 'EX_ac_e', 'EX_glc__D_e', max_or_min='max', carbon_uptake_direction=-1)
    print(yield_opt)

    # %% < get_yield_space_2d>
    print('----------get_yield_space_2d----------\n')
    fluxes_2d = get_yield_space_2d(model, production_rea_ids_2d=['BIOMASS_Ecoli_core_w_GAM', 'EX_ac_e'],
                                   carbon_source_rea_ids_2d=['EX_glc__D_e', 'EX_glc__D_e'], steps=steps,
                                   carbon_uptake_direction=-1, draw=draw)

    # %% < get_yield_space_multi >
    print('----------get_yield_space_multi----------\n')
    production_rea_ids_x = ['BIOMASS_Ecoli_core_w_GAM', 'EX_ac_e', ]
    production_rea_ids_y = ['BIOMASS_Ecoli_core_w_GAM', 'EX_ac_e', 'EX_for_e', 'EX_etoh_e', 'EX_lac__D_e', 'EX_succ_e']
    carbon_source_rea_id = 'EX_glc__D_e'

    yield_normalized_df = get_yield_space_multi(model, production_rea_ids_x, production_rea_ids_y, carbon_source_rea_id,
                                                steps=steps, carbon_uptake_direction=-1, draw=True)
