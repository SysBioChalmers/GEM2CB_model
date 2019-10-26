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

#%% <method EFMs + MYA> can not caculate
core_EFMs_z = np.genfromtxt('Case2_1_ecoli_core/e_coli_core_EFMs_standardized.csv', delimiter = ',')

EFM_all_points = core_EFMs_z.T
EFM_all_points = EFM_all_points[:,1:]   # modes x reas
experiment_datas = []


#%% <method EFVs + MYA> can not caculate
core_EFVs_z = np.genfromtxt('Case2_1_ecoli_core/e_coli_core_EFVs_standardized.csv', delimiter = ',')

EFV_all_points = core_EFVs_z.T
EFV_all_points = EFV_all_points[:,1:]   # modes x reas
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

step_of_biomass = 20
carbon_uptake_direction = -1
model.reactions.get_by_id(carbon_rea_id).bounds = (-10,0.001)

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

yield_rea_ids_name = yield_rea_ids

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

    points1 = ax.plot(EFM_all_points[:, 0], EFM_all_points[:, index], 'x', markerfacecolor='none', color='black',
                      alpha=0.5, label='EFM', markersize=3)
    points4 = ax.plot(EFV_all_points[:,0], EFV_all_points[:,index],'^',markerfacecolor='none',color = 'tab:blue',
                      alpha = 0.5,label = 'This_methd',markersize=3)
    points3 = ax.plot(our_all_points[:, 0], our_all_points[:, index], '^', markerfacecolor='none', color='tab:red',
                      alpha=0.8, label='EFM', markersize=5)
    points2 = ax.plot(FBAem_all_points[:, 0], FBAem_all_points[:, index], '.', color='tab:green',
                      alpha=1, label='FBA_mode', markersize=10)



    ax.set_ylabel(yield_rea_ids_name[index-1]+'/Glucose',fontsize = 12)

ax.set_xlabel('Yield Biomass/Glucose',fontsize = 12)
fig.legend((points1[0],points4[0],points2[0],points3[0]),('EFMs','EFVs','FBA modes','This study'),bbox_to_anchor=(0.55, 0.25), loc='upper left', borderaxespad=0.)
fig.show()


# fig = plt.figure()
# ax1 = fig.add_subplot(111)
# index = 5
# ax1.plot(EFM_all_points[:,0], EFM_all_points[:,index],'x',markerfacecolor='none',color = 'black',label = 'EFM',markersize=5)
# ax1.plot(EFV_all_points[:,0], EFV_all_points[:,index],'^',markerfacecolor='none',color = 'tab:green',alpha = 0.8,label = 'FBA_mode',markersize=5)
# ax1.plot(FBAem_all_points[:,0], FBAem_all_points[:,index],'o',markerfacecolor='none',color = 'tab:red',alpha = 0.8,label = 'FBA_mode',markersize=12)
#
# #ax1.plot(our_all_points[:,0], our_all_points[:,index],'v',markerfacecolor='none',color = 'tab:blue',alpha = 0.8,label = 'This_methd',markersize=12)
# ax1.legend()
# ax1.set_xlabel('Yield Biomass/Glucose',fontsize = 12)
# ax1.set_ylabel(yield_rea_ids[index-1]+'/Glucose',fontsize = 12)
# fig.show()



