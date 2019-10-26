#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-09-12

"""Step1_iML1515_test.py
:description : script
:param : 
:returns: 
:rtype: 
"""
import cobra
import My_def
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import statistics
from cobra.flux_analysis import production_envelope
from cobra.flux_analysis import flux_variability_analysis
import pandas as pd
import ConvexHull_yield
import GEM2pathways
from scipy.spatial import ConvexHull


iML1515 = cobra.io.read_sbml_model('../ComplementaryData/Case3_iML1515/iML1515.xml')


iML1515.reactions.get_by_id('EX_o2_e').bounds = (0.0,1000.0)
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
biomass_rea_id = 'BIOMASS_Ec_iML1515_core_75p37M'
carbon_rea_id = 'EX_glc__D_e'
yield_rea_ids = ['EX_ac_e','EX_for_e','EX_etoh_e','EX_lac__D_e','EX_succ_e',]


step_of_biomass = 20
carbon_uptake_direction = -1
model.reactions.get_by_id(carbon_rea_id).bounds = (-10,0.001)

# all modes
# yield_normalized_df = GEM2pathways.get_yield_space(model, biomass_rea_id,carbon_rea_id,yield_rea_ids,step_of_biomass,carbon_uptake_direction = carbon_uptake_direction,draw = True)
#
# yield_normalized_df = yield_normalized_df.loc[[carbon_rea_id] + [biomass_rea_id] + yield_rea_ids ,:]
#
# yield_normalized_df.to_csv('../ComplementaryData/Case3_iML1515/iML151_yield_normalized_df.csv', sep=',',index=True)


yield_normalized_df = pd.read_csv('../ComplementaryData/Case3_iML1515/iML151_yield_normalized_df.csv',index_col=0)

our_all_points = yield_normalized_df.values.T
our_all_points = our_all_points[ abs(our_all_points[:,0]) > 1e-10,:]
our_all_points = our_all_points[:,1:]       #exclude the carbon colume


our_all_points_2 = our_all_points.copy()
our_all_points = our_all_points_2.copy()
# our_all_points_2 = np.around(our_all_points,decimals=3)

initial_cons = [25.8461538462,0.0833333333,1.3986254296,0.6888564811,0.5165630156,0.0056365136,0.0547254725]
end_cons = [0.01,0.7583333333,19.7319587629,17.5388178777,19.9631932986,0.1603244432,2.8129612961]
initial_cons_2 = np.array(initial_cons)
end_cons_2 = np.array(end_cons)

yield_experiment = end_cons_2 - initial_cons_2
yield_experiment[:] = yield_experiment[:]/abs(yield_experiment[0])
yield_experiment = yield_experiment[1:]
yield_experiment[0] = yield_experiment[0] * 0.0070/ yield_experiment[0]

      # TODO experiment data!!!
experiment_datas = [ ]
experiment_datas =[yield_experiment]


qhull_options = 'Qt QJ Pp Qw Qx'    #'Qt QJ Pp Qw Qx'      #'QG0'mean expect point 0(index)
cutoff_persent = 1

our_indexes, our_hulls , our_weightss , our_estimated_datas, our_in_hulls = ConvexHull_yield.pipeline_mya(our_all_points, experiment_datas , qhull_options = qhull_options,methods = 4, cutoff_persent = cutoff_persent)
# [0.00934999 0.65153201 0.72955606 0.82466001 0.00598605 0.10897642]
# [0.0082365  0.68547344 0.68020769 0.76904425 0.00598605 0.10723878]

hull_cutoff_index_99 = 	 [2, 3, 4, 7, 8, 9, 12, 75, 77, 83, 87, 103, 120, 142, 150, 152, 155, 157, 190, 201]
hull_active_index_ac = 	 [  7  , 9  ,15 , 27 , 37 , 47 , 67 , 75 , 77, 130, 142]


our_all_points_hull_1 = our_all_points[our_indexes[0],:]
our_all_points_hull_99 = our_all_points[hull_cutoff_index_99,:]
our_all_points_hull_ac = our_all_points[hull_active_index_ac,:]

#%% <plot initial yield>: TODO not divided by carbon, so not yield

yield_rea_ids_name = ['Acetate', 'Formate', 'Ethanol', 'Lacate', 'Succinate']
markersize_list = [10,10,15]
lw_sizelist =  [2,2,10]
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

    points_exp = ax.plot(experiment_datas [0][0], experiment_datas [0][index],          '^', color='red',
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


        xy_our_hull = our_all_points_hull_ac[:, [0,index]]
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
fig.legend((points_our[0],points_exp[0],line_our_1[0],line_our_99[0] ,line_our_ac[0]),('All pathways','Experiment data','Convex hull','Convex hull %99','Convex hull with data'),bbox_to_anchor=(0.55, 0.35), loc='upper left', borderaxespad=0.,prop={'size': 20})

fig.show()

fig.savefig('../ComplementaryData/Case3_iML1515/4.png',format ='png')









