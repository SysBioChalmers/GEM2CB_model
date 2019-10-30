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
import cobra
import matplotlib.pyplot as plt
import numpy as np
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


#%% <EFMs > caculated by matlab emftool. and standardized by step0

em_z = np.genfromtxt('Case1_ecoli_reduced/ecoli_reduced_EFMs_standardized.csv', delimiter = ',')
EFMs_all_points = em_z.T
EFMs_all_points = EFMs_all_points[:,1:]   # modes x reas


#%% <FBA modes> from reference and standardized by step0

FBAMs_z = np.genfromtxt('Case1_ecoli_reduced/ecoli_reduced_FBAMs_standardized.csv', delimiter = ',')
FBAem_all_points = FBAMs_z.T     # modes x rea
FBAem_all_points = FBAem_all_points[:,1:]


#%% <Our method>

# load a GE
ecoli_reduced_model = cobra.io.read_sbml_model('Case1_ecoli_reduced/ecoli_reduced_model.xml')
model = ecoli_reduced_model.copy()

production_rea_ids_x = ['R12']
production_rea_ids_y = ['R9', 'EX_FOR', 'R10', 'R8', 'R6']
carbon_source_rea_id = 'R1'
model.reactions.get_by_id(carbon_source_rea_id).bounds = (0.001, 10)  # set carbon lb as -10
steps = 20
carbon_uptake_direction = 1

# all modes
yield_normalized_df, fluxes_all = GEM2pathways.get_yield_space_multi(model, production_rea_ids_x, production_rea_ids_y,
                                                                     carbon_source_rea_id,
                                                                     steps=steps,
                                                                     carbon_uptake_direction=carbon_uptake_direction,
                                                                     draw=True)

our_all_points = yield_normalized_df.values.T
our_all_points = our_all_points[ abs(our_all_points[:,0]) > 1e-10,:]
our_all_points = our_all_points[:, 1:]  # exclude the carbon colume!!!

# %% <MYA >

experiment_datas = [ ]      # TODO experiment data!!!
qhull_options = 'Qt QJ Pp Qw Qx'      #'QG0'mean expect point 0(index)
cutoff_persent = 0.99

our_indexes, our_hulls , our_weightss , our_estimated_datas, our_in_hulls = ConvexHull_yield.pipeline_mya(our_all_points, experiment_datas , qhull_options = qhull_options,methods = 4, cutoff_persent = cutoff_persent)

len(our_indexes[1])
# EFMs_indexes, EFMs_hulls , EFMs_weightss , EFMs_estimated_datas, EFMs_in_hulls = \
#     ConvexHull_yield.pipeline_mya(EFM_all_points, experiment_datas , qhull_options = 'Qt QJ Pp Qw Qx', cutoff_persent = cutoff_persent )

# EFMs_hull = ConvexHull(EFMs_all_points,qhull_options = qhull_options )
# FBAem_hull = ConvexHull(FBAem_all_points,qhull_options = qhull_options )
our_hull = ConvexHull(our_all_points,qhull_options = qhull_options )



#%% <plot initial yield>: TODO not divided by carbon, so not yield

yield_rea_ids_name = ['EX_ac_e','EX_for_e','EX_etoh_e','EX_lac__D_e','EX_succ_e',]

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
    xy_FBAem = FBAem_all_points[:, [0,index]]
    xy_our = our_all_points[:, [0,index]]

    hull_EFMs = ConvexHull(xy_EFMs,qhull_options = qhull_options )
    hull_FBAem = ConvexHull(xy_FBAem,qhull_options = qhull_options )
    hull_our = ConvexHull(xy_our,qhull_options = qhull_options )

    points_our = ax.plot(xy_our[:,0], xy_our[:,1],                 '^',markerfacecolor='none',color = 'tab:blue',
                      alpha = 1,label = 'This_methd',markersize=5)
    for simplex in hull_our.simplices: 
        line_our = ax.plot(xy_our[simplex, 0], xy_our[simplex, 1], '^--',markerfacecolor='none',color = 'tab:blue',
                      alpha = 0.5,label = 'This_methd',markersize=5)
    
    points_FBAem = ax.plot(xy_FBAem[:,0], xy_FBAem[:,1],                'o', markerfacecolor='none', color='black',
                      alpha=1, label='FBA_mode', markersize=10)
    for simplex in hull_FBAem.simplices: 
        line_FBAem = ax.plot(xy_FBAem[simplex, 0], xy_FBAem[simplex, 1], 'o--', markerfacecolor='none', color='black',
                      alpha=0.5, label='FBA_mode', markersize=10)
    
    points_EFMs = ax.plot(xy_EFMs[:,0], xy_EFMs[:,1],                'x', markerfacecolor='none', color='tab:red',
                      alpha=1, label='EFMs', markersize=10)
    for simplex in hull_EFMs.simplices: 
        line_EFMs = ax.plot(xy_EFMs[simplex, 0], xy_EFMs[simplex, 1], 'x-', markerfacecolor='none', color='tab:red',
                      alpha=0.5, label='EFMs', markersize=10)
    
    
    ax.set_ylabel(yield_rea_ids_name[index-1]+'/Glucose',fontsize = 12)

ax.set_xlabel('Yield Biomass/Glucose',fontsize = 12)
fig.legend((line_EFMs[0],line_FBAem[0],line_our[0]),('EFMs','FBA modes','This study'),bbox_to_anchor=(0.55, 0.25), loc='upper left', borderaxespad=0.)
fig.show()













