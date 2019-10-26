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

yield_rea_ids = ['EX_ac_e','EX_for_e','EX_etoh_e','EX_lac__D_e','EX_succ_e',]

'''
reaids = ['R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12']
yield_rea_ids = ['R1','R12','R9','R7','R10','R8','R6']
yield_indexes = [reaids.index(i) for i in yield_rea_ids ]

#%% <method EFM + MYA> caculated by matlab emftool. and standardized by step0

em_z = np.genfromtxt('Case1_ecoli_reduced/ecoli_reduced_EFMs_standardized.csv', delimiter = ',')

# MYA hull & #MYA experiment data
EFM_all_points = em_z.T
EFM_all_points = EFM_all_points[:,1:]   # modes x reas
experiment_datas = []

#EFM_indexes, EFM_hulls , EFM_weightss , EFM_estimated_datas, EFM_in_hulls = ConvexHull_yield.pipeline_mya(EFM_all_points, experiment_datas , qhull_options = 'Qt QJ Pp Qw Qx', cutoff_persent = 0.99,methods = 1)


#%% <FBA modes> from reference and standardized by step0

FBAMs_z = np.genfromtxt('Case1_ecoli_reduced/ecoli_reduced_FBAMs_standardized.csv', delimiter = ',')

FBAem_all_points = FBAMs_z.T     # modes x rea
FBAem_all_points = FBAem_all_points[:,1:]
experiment_datas = []
# FBAem_indexes, FBAem_hulls , FBAem_weightss , FBAem_estimated_datas, FBAem_in_hulls = ConvexHull_yield.pipeline_mya(\
#     FBAem_all_points, experiment_datas , qhull_options = 'Qt QJ Pp Qw Qx', cutoff_persent = 0.99,methods = 4)


#%% <Our method>

# load a GEM
ecoli_reduced_model = cobra.io.read_sbml_model('Case1_ecoli_reduced/ecoli_reduced_model.xml')
model = ecoli_reduced_model.copy()

# cobra.io.write_sbml_model(ecoli_reduced_model,'Case1_ecoli_reduced/ecoli_reduced_model.xml')
biomass_rea_id = 'R12'
carbon_rea_id = 'R1'
yield_rea_ids = ['R9','EX_FOR','R10','R8','R6']

step_of_biomass = 20
carbon_uptake_direction = 1
# ecoli_reduced_model.solver = 'glpk'
model.reactions.get_by_id(carbon_rea_id).bounds = (0.001,10)     #set carbon lb as -10

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

# our_indexes, our_hulls , our_weightss , our_estimated_datas, our_in_hulls = ConvexHull_yield.pipeline_mya(our_all_points, experiment_datas , qhull_options = 'Qt QJ Pp Qw Qx', cutoff_persent = 0.99)


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

    points3 = ax.plot(our_all_points[:,0], our_all_points[:,index],'^',markerfacecolor='none',color = 'tab:blue',
                      alpha = 1,label = 'This_methd',markersize=5)
    points2 = ax.plot(FBAem_all_points[:, 0], FBAem_all_points[:, index], 'o', markerfacecolor='none', color='black',
                      alpha=1, label='FBA_mode', markersize=10)
    points1 = ax.plot(EFM_all_points[:, 0], EFM_all_points[:, index], 'x', markerfacecolor='none', color='tab:red',
                      alpha=1, label='EFM', markersize=10)
    ax.set_ylabel(yield_rea_ids_name[index-1]+'/Glucose',fontsize = 12)

ax.set_xlabel('Yield Biomass/Glucose',fontsize = 12)
fig.legend((points1[0],points2[0],points3[0]),('EFMs','FBA modes','This study'),bbox_to_anchor=(0.55, 0.25), loc='upper left', borderaxespad=0.)
fig.show()













