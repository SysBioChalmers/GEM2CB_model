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
import My_def
os.chdir('../ComplementaryData/')

#%% <method EFMs + MYA>
ECC2comp_EFVs_noO2 = sio.loadmat('ECC2comp/ECC2comp_EFVs_noO2.mat')

ECC2comp_EFVs_noO2_z = ECC2comp_EFVs_noO2['elmoden']

# rea_ids = ECC2comp_EFVs_noO2['mode_rates_names']
# biomass_rea_id = 'Growth'
# carbon_rea_ids = 'R_GlcUp'
# yield_rea_ids = ['R_AcEx','R_FormEx','R_EthEx','R_LacEx','R_SuccEx',]
#
# data_cloumes = [carbon_rea_ids] + [biomass_rea_id] + yield_rea_ids
# data_cloumes_indexs = data_cloumes
#
# for i in range(len(rea_ids)):
#     for ii in range(len(data_cloumes)):
#         if type(data_cloumes[ii]) == str:
#             if data_cloumes[ii] in rea_ids[i]:
#                 data_cloumes_indexs[ii] = i

data_cloumes_indexs = [29, 0, 2, 21, 17, 37, 59]






#%% <method EFVs + MYA>  # EFVs

ECC2comp_EFVs_noO2 = sio.loadmat('ECC2comp/ECC2comp_EFVs_noO2.mat')

ECC2comp_EFVs_noO2_z = ECC2comp_EFVs_noO2['elmoden']

# rea_ids = ECC2comp_EFVs_noO2['mode_rates_names']
# biomass_rea_id = 'Growth'
# carbon_rea_ids = 'R_GlcUp'
# yield_rea_ids = ['R_AcEx','R_FormEx','R_EthEx','R_LacEx','R_SuccEx',]
#
# data_cloumes = [carbon_rea_ids] + [biomass_rea_id] + yield_rea_ids
# data_cloumes_indexs = data_cloumes
#
# for i in range(len(rea_ids)):
#     for ii in range(len(data_cloumes)):
#         if type(data_cloumes[ii]) == str:
#             if data_cloumes[ii] in rea_ids[i]:
#                 data_cloumes_indexs[ii] = i

data_cloumes_indexs = [29, 0, 2, 21, 17, 37, 59]



# normlize
em_z  = ECC2comp_EFVs_noO2_z.T[data_cloumes_indexs,:]
em_z[:,1:]  = em_z[:,1:]/abs(em_z[0,1:])

#MYA hull & #MYA experiment data
EFM_all_points = em_z[1:,:]
EFM_all_points = EFM_all_points.T
EFM_all_points = EFM_all_points[1:,:]


#%% <Our method>

'''
ECC2comp = cobra.io.read_sbml_model('../ComplementaryData/ECC2comp/ECC2_glucose_standard.xml')
#e_coli_core_mat = cobra.io.load_matlab_model('../ComplementaryData/ECC2comp/ECC2_glucose_standard.mat')
#My_def.io_file.model2txt(ECC2comp,'../ComplementaryData/ECC2comp/ECC2_glucose_standard_xml.txt')

ECC2comp.reactions.AcUp.bounds = (0.0,0.0)
ECC2comp.reactions.GlycUp.bounds = (0.0,0.0)
ECC2comp.reactions.SuccUp.bounds = (0.0,0.0)

ECC2comp.reactions.GlcUp.bounds = (0.0,10.0)
ECC2comp.reactions.O2Up.bounds = (0.0,5.0)
ECC2comp.reactions.ATPM.bounds = (3.15,1000.0)

cobra.io.write_sbml_model(ECC2comp,'../ComplementaryData/ECC2comp/ECC2_standard.xml')

'''
#  standard:


ECC2comp = cobra.io.read_sbml_model('../ComplementaryData/ECC2comp/ECC2_standard.xml')
ECC2comp.reactions.O2Up.bounds = (0.0,0.0)

#My_def.io_file.model2txt(e_coli_core_mat,'../ComplementaryData/ECC2comp/ECC2_glucose_standard_mat.txt')

#e_coli_core.reactions.get_by_id('EX_o2_e').bounds = (0.0,1000.0)

biomass_rea_id = 'Growth'
carbon_rea_ids = 'GlcUp'
yield_rea_ids = ['EX_ac_ex','EX_for_ex','EX_etoh_ex','EX_lac_ex','EX_succ_ex',]

step_of_biomass = 20


# for carbon_rea_id in carbon_rea_ids:        #set carbon lb as -10
#     ECC2comp.reactions.get_by_id(carbon_rea_id).bounds = (-10.0,0.0)
model = ECC2comp


# all modes
yield_normalized_df = GEM2pathways.get_fba_yield_df(model, biomass_rea_id,carbon_rea_ids,yield_rea_ids,step_of_biomass)


#yield_normalized_df = yield_normalized_df.fillna(0)


#MYA hull#MYA experiment data
data_cloumes = [carbon_rea_ids] + [biomass_rea_id] + yield_rea_ids

yield_normalized_df = yield_normalized_df.loc[data_cloumes,:]
our_all_points = yield_normalized_df.values.T
our_all_points = our_all_points[:,1:]       #exclude the carbon colume


#%% <plot initial yield>: TODO not divided by carbon, so not yield


fig = plt.figure()
ax1 = fig.add_subplot(111)

index = 1
ax1.plot(EFM_all_points[:,0], EFM_all_points[:,index],'x',markerfacecolor='none',color = 'tab:blue',label = 'EFM',markersize=3,alpha = 0.8)
# ax1.plot(FBAem_all_points[:,0], FBAem_all_points[:,index],'o',markerfacecolor='none',color = 'tab:red',alpha = 0.8,label = 'FBA_mode',markersize=12)
ax1.plot(our_all_points[:,0], our_all_points[:,index],'v',markerfacecolor='none',color = 'tab:red',alpha = 0.8,label = 'This_methd',markersize=8)
ax1.legend()
ax1.set_xlabel('Yield Biomass/Glucose',fontsize = 12)
ax1.set_ylabel(yield_rea_ids[index-1]+'/Glucose',fontsize = 12)
fig.show()





