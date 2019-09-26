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

-METINT
G6P T3P PEP PYR AcCoA ATP ADP NADH NAD CoA

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

reaids = ['R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12']
yield_rea_ids = ['EX_ac_e','EX_for_e','EX_etoh_e','EX_lac__D_e','EX_succ_e',]
yield_rea_ids = ['R1','R12','R9','R7','R10','R8','R6']
'''

#%% <method EFM + MYA> caculated by matlab emftool. y1: biomass y2 acet
#EFM    from
reaids = ['R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12']
yield_rea_ids = ['R1','R12','R9','R7','R10','R8','R6']
yield_indexes = [reaids.index(i) for i in yield_rea_ids ]

ecoli_reduced_z = sio.loadmat('ecoli_reduced_z.mat')
ecoli_reduced_efms = ecoli_reduced_z['reduce_z']

em_z  = ecoli_reduced_efms[yield_indexes,:]   #normalization by carbon
em_z[:,1:]  = em_z[:,1:]/abs(em_z[0,1:])

#MYA hull & #MYA experiment data
EFM_all_points = em_z[1:,:]
EFM_all_points = EFM_all_points.T
EFM_all_points = EFM_all_points[1:,:]
experiment_datas = []
#EFM_indexes, EFM_hulls , EFM_weightss , EFM_estimated_datas, EFM_in_hulls = ConvexHull_yield.pipeline_mya(EFM_all_points, experiment_datas , qhull_options = 'Qt QJ Pp Qw Qx', cutoff_persent = 0.99,methods = 1)





#%% <FBA mode>

#load a GEM
# reaids = ['R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12']
# reaids_all = reaids +['EX_GLU','EX_SUC','EX_FOR','EX_ACT','EX_LAC','EX_ETH','EX_B','EX_CO2','EX_H2']
ecoli_reduced_model = cobra.io.read_sbml_model('ecoli_reduced_model.xml')

# load modes
ecoli_reduced_FBA_ems = np.genfromtxt('FBA_em.csv', delimiter=',')
FBAem_z = ecoli_reduced_FBA_ems[:,yield_indexes]
FBAem_z = FBAem_z.T
FBAem_z[:,:-1] = FBAem_z[:,:-1]/FBAem_z[0,:-1]
FBAem_z = FBAem_z.T

FBAem_all_points = FBAem_z[:,1:]
FBAem_all_points = FBAem_all_points[:-1,:]
experiment_datas = []
# FBAem_indexes, FBAem_hulls , FBAem_weightss , FBAem_estimated_datas, FBAem_in_hulls = ConvexHull_yield.pipeline_mya(\
#     FBAem_all_points, experiment_datas , qhull_options = 'Qt QJ Pp Qw Qx', cutoff_persent = 0.99,methods = 4)


#%% <Our method>

#load a GEM
# reaids = ['R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12']
# reaids_all = reaids +['EX_GLU','EX_SUC','EX_FOR','EX_ACT','EX_LAC','EX_ETH','EX_B','EX_CO2','EX_H2']
ecoli_reduced_model = cobra.io.read_sbml_model('ecoli_reduced_model.xml')


biomass_rea_id = 'EX_B'
carbon_rea_ids = ['EX_GLU']
yield_rea_ids = ['EX_ACT','EX_FOR','EX_ETH','EX_LAC','EX_SUC']
# yield_rea_ids = ['EX_ac_e','EX_for_e','EX_etoh_e','EX_lac__D_e','EX_succ_e',]
step_of_biomass = 10
for carbon_rea_id in carbon_rea_ids:        #set carbon lb as -10
    ecoli_reduced_model.reactions.get_by_id(carbon_rea_id).bounds = (-10.0,0)
model = ecoli_reduced_model
# all modes
yield_normalized_df = GEM2pathways.get_fba_yield_df(model, biomass_rea_id,carbon_rea_ids,yield_rea_ids,step_of_biomass)
#yield_normalized_df = yield_normalized_df.fillna(0)


#MYA hull#MYA experiment data
yield_normalized_df = yield_normalized_df.loc[['EX_GLU','EX_B','EX_ACT','EX_FOR','EX_ETH','EX_LAC','EX_SUC'],:]
our_all_points = yield_normalized_df.values.T
our_all_points = our_all_points[:,1:]       #exclude the carbon colume

experiment_datas = [ ]      # TODO experiment data!!!
qhull_options = 'Qt QJ Pp Qw Qx'      #'QG0'mean expect point 0(index)
cutoff_persent = 0.99

# our_indexes, our_hulls , our_weightss , our_estimated_datas, our_in_hulls = ConvexHull_yield.pipeline_mya(our_all_points, experiment_datas , qhull_options = 'Qt QJ Pp Qw Qx', cutoff_persent = 0.99)


#%% <plot initial yield>: TODO not divided by carbon, so not yield


fig = plt.figure()
ax1 = fig.add_subplot(111)

index = 5
ax1.plot(EFM_all_points[:,0], EFM_all_points[:,index],'x',markerfacecolor='none',color = 'black',label = 'EFM',markersize=12)
ax1.plot(FBAem_all_points[:,0], FBAem_all_points[:,index],'o',markerfacecolor='none',color = 'tab:red',alpha = 0.8,label = 'FBA_mode',markersize=12)
ax1.plot(our_all_points[:,0], our_all_points[:,index],'v',markerfacecolor='none',color = 'tab:blue',alpha = 0.8,label = 'This_methd',markersize=12)
ax1.legend()
ax1.set_xlabel('Yield Biomass/Glucose',fontsize = 12)
ax1.set_ylabel(yield_rea_ids[index-1]+'/Glucose',fontsize = 12)
fig.show()





#
# def trim_axs(axs, N):
#     """little helper to massage the axs list to have correct length..."""
#     axs = axs.flat
#     for ax in axs[N:]:
#         ax.remove()
#     return axs[:N]
#
# cols = 3
# rows = len(yield_rea_ids) // cols + 1
# figsize = (10, 8)
# fig, axs = plt.subplots(cols, rows,figsize=figsize)
# axs = trim_axs(axs, len(yield_rea_ids))
#
# for ax , exmet_reaction in zip(axs,yield_rea_ids):
#     temp_columns_max = [column for column in yield_normalized_df.columns if (exmet_reaction in column) and ('max' in column)]
#     temp_columns_min = [column for column in yield_normalized_df.columns if (exmet_reaction in column) and ('min' in column)]
#     x = yield_normalized_df.loc[biomass_rea_id,temp_columns_max].values
#     y_max = yield_normalized_df.loc[exmet_reaction,temp_columns_max]
#     y_min = yield_normalized_df.loc[exmet_reaction,temp_columns_min]
#
#     ax.plot(x, y_max,'x-',color = 'black',alpha=0.5)
#     ax.plot(x, y_min,'x-',color = 'black',alpha=0.5)
#     # temp_df.plot.area(x='biomass', y=['maximum','minimum'],label=['max', 'max'],color=['r', 'w'],color=['tab:blue', 'b'],stacked=False);
#     ax.set_ylabel(exmet_reaction)
# ax.set_xlabel('biomass')
# fig.show()












