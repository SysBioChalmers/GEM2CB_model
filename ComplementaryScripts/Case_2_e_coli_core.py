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
import matplotlib.pyplot as plt

import ConvexHull_yield
import GEM2pathways
import scipy.io as sio
import os
import numpy as np
from scipy.spatial import ConvexHull

os.chdir('../ComplementaryData/')

#%% <method EFMs > can not caculate
core_EFMs_z = np.genfromtxt('Case2_1_ecoli_core/e_coli_core_EFMs_standardized.csv', delimiter = ',')

EFMs_all_points = core_EFMs_z.T
EFMs_all_points = EFMs_all_points[:,1:]   # modes x reas
experiment_datas = []


#%% <method EFVs > can not caculate
core_EFVs_z = np.genfromtxt('Case2_1_ecoli_core/e_coli_core_EFVs_standardized.csv', delimiter = ',')

EFVs_all_points = core_EFVs_z.T
EFVs_all_points = EFVs_all_points[:,1:]   # modes x reas
experiment_datas = []



#%% <FBA mode>

core_FBAMs_z = np.genfromtxt('Case2_1_ecoli_core/ecoli_core_FBAMs_standardized.csv', delimiter = ',')

FBAem_all_points = core_FBAMs_z.T     # modes x rea
FBAem_all_points = FBAem_all_points[:,1:]
experiment_datas = []


#%% <Our method>
e_coli_core = cobra.io.read_sbml_model('../ComplementaryData/Case2_1_ecoli_core/e_coli_core.xml')
#e_coli_core.reactions.get_by_id('EX_o2_e').bounds = (0.0,1000.0)
model = e_coli_core.copy()

biomass_rea_id = 'BIOMASS_Ecoli_core_w_GAM'
carbon_rea_id = 'EX_glc__D_e'
yield_rea_ids = ['EX_ac_e','EX_for_e','EX_etoh_e','EX_lac__D_e','EX_succ_e',]

step_of_biomass = 25
carbon_uptake_direction = -1
model.reactions.get_by_id(carbon_rea_id).bounds = (-10,-0.001)

# all modes
yield_normalized_df = GEM2pathways.get_yield_space(model, biomass_rea_id,carbon_rea_id,yield_rea_ids,step_of_biomass,carbon_uptake_direction = carbon_uptake_direction,draw = False)

#MYA hull#MYA experiment data
yield_normalized_df = yield_normalized_df.loc[[carbon_rea_id] + [biomass_rea_id] + yield_rea_ids ,:]

our_all_points = yield_normalized_df.values.T
our_all_points = our_all_points[ abs(our_all_points[:,0]) > 1e-10,:]
our_all_points = our_all_points[:,1:]       #exclude the carbon colume

experiment_datas = [ ]      # TODO experiment data!!!
qhull_options = 'Qt QJ Pp Qw Qx'      #'QG0'mean expect point 0(index)
cutoff_persent = 0.99


#%% <plot initial yield>: TODO not divided by carbon, so not yield

# yield_rea_ids_name = yield_rea_ids
yield_rea_ids_name = ['Acetate', 'Formate', 'Ethanol', 'Lacate', 'Succinate']


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
    xy_FBAem = FBAem_all_points[:, [0,index]]
    xy_our = our_all_points[:, [0,index]]

    colors_list = ['blue','teal','tab:red','tab:orange']

    points_EFMs = ax.plot(xy_EFMs[:,0], xy_EFMs[:,1],       '.', markerfacecolor='none', color=colors_list[0],
                    alpha = 0.5, markersize=3)
    points_EFVs = ax.plot(xy_EFVs[:,0], xy_EFVs[:,1],       '.', markerfacecolor='none', color=colors_list[1],
                    alpha = 0.5, markersize=3)

    points_our = ax.plot(xy_our[:,0], xy_our[:,1],          '^', markerfacecolor='none', color=colors_list[2],
                    alpha=0.8,  markersize=6)
    points_FBAem = ax.plot(xy_FBAem[:,0], xy_FBAem[:,1],    's', color=colors_list[3],
                    alpha=1,  markersize=8)

    draw_hull = True

    if draw_hull:
        hull_EFMs = ConvexHull(xy_EFMs,qhull_options = qhull_options )
        hull_EFVs = ConvexHull(xy_EFVs,qhull_options = qhull_options )
        hull_FBAem = ConvexHull(xy_FBAem,qhull_options = qhull_options )
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

        for simplex in hull_FBAem.simplices:
            line_FBAem = ax.plot(xy_FBAem[simplex, 0], xy_FBAem[simplex, 1], 'o-', markerfacecolor='none', color=colors_list[3],
                            alpha = 0.5, markersize=10,lw=8)

    else:
        line_EFMs = points_EFMs
        line_EFVs = points_EFVs
        line_our = points_our
        line_FBAem = points_FBAem

    ax.set_ylabel(yield_rea_ids_name[index-1]+' / Glc',fontsize = 18)

ax.set_xlabel('Yield Biomass/Glucose',fontsize = 18)
fig.legend((line_EFMs[0],line_EFVs[0],line_our[0],line_FBAem[0]),('EFMs','EFVs','This study','FBA modes'),bbox_to_anchor=(0.6, 0.3), loc='upper left', borderaxespad=0.,prop={'size': 18})
fig.show()

fig.savefig('Case2_1_ecoli_core/4.png',format ='png')


# %% ac

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
xy_FBAem = FBAem_all_points[:, [0,index]]
xy_our = our_all_points[:, [0,index]]

colors_list = ['blue','teal','tab:red','tab:orange']

points_EFMs = ax.plot(xy_EFMs[:,0], xy_EFMs[:,1],       '.', markerfacecolor='none', color=colors_list[0],
                alpha = 0.5, markersize=3)
points_EFVs = ax.plot(xy_EFVs[:,0], xy_EFVs[:,1],       '.', markerfacecolor='none', color=colors_list[1],
                alpha = 0.5, markersize=3)

points_our = ax.plot(xy_our[:,0], xy_our[:,1],          '^', markerfacecolor='none', color=colors_list[2],
                alpha=0.8,  markersize=6)
points_FBAem = ax.plot(xy_FBAem[:,0], xy_FBAem[:,1],    's', color=colors_list[3],
                alpha=1,  markersize=8)

draw_hull = True

if draw_hull:
    hull_EFMs = ConvexHull(xy_EFMs,qhull_options = qhull_options )
    hull_EFVs = ConvexHull(xy_EFVs,qhull_options = qhull_options )
    hull_FBAem = ConvexHull(xy_FBAem,qhull_options = qhull_options )
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

    for simplex in hull_FBAem.simplices:
        line_FBAem = ax.plot(xy_FBAem[simplex, 0], xy_FBAem[simplex, 1], 'o-', markerfacecolor='none', color=colors_list[3],
                        alpha = 0.5, markersize=10,lw=2)

else:
    line_EFMs = points_EFMs
    line_EFVs = points_EFVs
    line_our = points_our
    line_FBAem = points_FBAem

ax.set_ylabel(yield_rea_ids_name[index-1]+' / Glc',fontsize = 18)

ax.set_xlabel('Yield Biomass/Glucose',fontsize = 18)
fig.legend((line_EFMs[0],line_EFVs[0],line_our[0],line_FBAem[0]),('EFMs','EFVs','This study','FBA modes'),bbox_to_anchor=(0.6, 0.95), loc='upper left', borderaxespad=0.,prop={'size': 16})
fig.show()

fig.savefig('Case2_1_ecoli_core/0_ac.png',format ='png')