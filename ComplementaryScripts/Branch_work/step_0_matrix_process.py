#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 9/27/19

"""step0_process_initial_z.py
:description : script
:param : 
:returns:  All Z file will be rea x modes
:rtype: 
"""
import scipy.io as sio
import cobra
import os
import numpy as np

os.chdir('../../ComplementaryData/')


# %% case 1 Ecoli reduced 12 reactions

reaids = ['R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12']
yield_rea_ids = ['R1','R12','R9','R7','R10','R8','R6']
yield_indexes = [reaids.index(i) for i in yield_rea_ids ]

# EFMs  emftool

ecoli_reduced_z = sio.loadmat('Case1_ecoli_reduced/ecoli_reduced_z.mat')
ecoli_reduced_efms = ecoli_reduced_z['reduce_z']    # rea x modes

ecoli_reduced_efms[6,:] = ecoli_reduced_efms[6,:]-ecoli_reduced_efms[10,:]  # EX_for_e = R7 - R11

em_z  = ecoli_reduced_efms[yield_indexes,:]   # reas x modes    normalization by carbon:
em_z = em_z[:,abs(em_z[0,:]) > 1e-10]   # pathway must uptake glc
em_z[:,:]  = em_z[:,:]/abs(em_z[0,:])  # normalization

np.savetxt('Case1_ecoli_reduced/ecoli_reduced_EFMs_standardized.csv', em_z, delimiter=',')
reduced_EFMs = np.genfromtxt('Case1_ecoli_reduced/ecoli_reduced_EFMs_standardized.csv', delimiter = ',')


# FBA modes from reference

ecoli_reduced_FBA_ems = np.genfromtxt('Case1_ecoli_reduced/FBA_em_reduced.csv', delimiter=',')
ecoli_reduced_FBA_ems = ecoli_reduced_FBA_ems.T     # rea x modes

ecoli_reduced_FBA_ems[6,:] = ecoli_reduced_FBA_ems[6,:]-ecoli_reduced_FBA_ems[10,:]     # EX_for_e = R7 - R11

FBAem_z = ecoli_reduced_FBA_ems[yield_indexes,:]    # rea x modes
FBAem_z = FBAem_z[:,FBAem_z[0,:] > 1e-10]   # pathway must uptake glc
FBAem_z[:,:] = FBAem_z[:,:]/abs(FBAem_z[0,:])

np.savetxt('Case1_ecoli_reduced/ecoli_reduced_FBAMs_standardized.csv', FBAem_z, delimiter=',')
reduced_FBAMs = np.genfromtxt('Case1_ecoli_reduced/ecoli_reduced_FBAMs_standardized.csv', delimiter = ',')

# %% <experiment data process>












# %% case 2.1 Ecoli core model 95 reactions

# EFMs  CNA

ecoli_core_EFMs = sio.loadmat('Case2_1_ecoli_core/e_coli_core_EFMs.mat')
ecoli_core_EFMs_z = ecoli_core_EFMs['elmoden']  # modes x rea
rea_ids = ecoli_core_EFMs['mode_rates_names']

biomass_rea_id = 'BIOMASS_Ecoli_core_w_GAM'
carbon_rea_ids = 'EX_glc__D_e'
yield_rea_ids = ['EX_ac_e','EX_for_e','EX_etoh_e','EX_lac__D_e','EX_succ_e',]

data_cloumes = [carbon_rea_ids] + [biomass_rea_id] + yield_rea_ids
data_cloumes_indexs = data_cloumes

for i in range(len(rea_ids)):
    for ii in range(len(data_cloumes)):
        if type(data_cloumes[ii]) == str:
            if data_cloumes[ii] in rea_ids[i]:
                data_cloumes_indexs[ii] = i

# data_cloumes_indexs = [25, 12, 19, 24, 23, 29, 34]

em_z  = ecoli_core_EFMs_z.T[data_cloumes_indexs,:]  # rea x modes
em_z = em_z[:,abs(em_z[0,:]) > 1e-10]   # pathway must uptake glc
em_z[:,:]  = em_z[:,:]/abs(em_z[0,:])  # normalization

np.savetxt('Case2_1_ecoli_core/e_coli_core_EFMs_standardized.csv', em_z, delimiter=',')
core_EFMs = np.genfromtxt('Case2_1_ecoli_core/e_coli_core_EFMs_standardized.csv', delimiter = ',')

# EFVs  CNA

ecoli_core_EFVs = sio.loadmat('Case2_1_ecoli_core/e_coli_core_EFvs.mat')
ecoli_core_EFVs_z = ecoli_core_EFVs['elmoden']  # modes x rea
rea_ids = ecoli_core_EFVs['mode_rates_names']

biomass_rea_id = 'BIOMASS_Ecoli_core_w_GAM'
carbon_rea_ids = 'EX_glc__D_e'
yield_rea_ids = ['EX_ac_e','EX_for_e','EX_etoh_e','EX_lac__D_e','EX_succ_e',]

data_cloumes = [carbon_rea_ids] + [biomass_rea_id] + yield_rea_ids
data_cloumes_indexs = data_cloumes

for i in range(len(rea_ids)):
    for ii in range(len(data_cloumes)):
        if type(data_cloumes[ii]) == str:
            if data_cloumes[ii] in rea_ids[i]:
                data_cloumes_indexs[ii] = i

# data_cloumes_indexs = [25, 12, 19, 24, 23, 29, 34]

em_z  = ecoli_core_EFVs_z.T[data_cloumes_indexs,:]  # rea x modes
em_z = em_z[:,abs(em_z[0,:]) > 1e-10]   # pathway must uptake glc
em_z[:,:]  = em_z[:,:]/abs(em_z[0,:])  # normalization

np.savetxt('Case2_1_ecoli_core/e_coli_core_EFVs_standardized.csv', em_z, delimiter=',')
core_EFVs = np.genfromtxt('Case2_1_ecoli_core/e_coli_core_EFVs_standardized.csv', delimiter = ',')


# FBA modes from reference

ecoli_core_FBA_ems = np.genfromtxt('Case2_1_ecoli_core/FBA_em_core.csv', delimiter=',')
# ecoli_core_FBA_ems = ecoli_core_FBA_ems.T     # rea x modes

biomass_rea_id = 'BIOMASS_Ecoli_core_w_GAM'
carbon_rea_ids = 'EX_glc__D_e'
yield_rea_ids = ['EX_ac_e','EX_for_e','EX_etoh_e','EX_lac__D_e','EX_succ_e',]
# glc: 89 + 102
# biom: 105
# EX_ACT : 96 -97
# for 101
# eth 95
# LAC 99
# succ 98
glc = ecoli_core_FBA_ems[88,:] + ecoli_core_FBA_ems[101,:]
biom = ecoli_core_FBA_ems[104,:]
ac = ecoli_core_FBA_ems[95,:] - ecoli_core_FBA_ems[96,:]
ex_for = ecoli_core_FBA_ems[100,:]
ex_eth = ecoli_core_FBA_ems[94,:]
ex_lac = ecoli_core_FBA_ems[98,:]
ex_succ = ecoli_core_FBA_ems[97,:]

FBAem_z = np.array([glc,biom,ac,ex_for,ex_eth,ex_lac,ex_succ])      # rea x modes

FBAem_z = FBAem_z[:,FBAem_z[0,:] > 1e-10]   # pathway must uptake glc
FBAem_z[:,:] = FBAem_z[:,:]/abs(FBAem_z[0,:])

np.savetxt('Case2_1_ecoli_core/ecoli_core_FBAMs_standardized.csv', FBAem_z, delimiter=',')
core_FBAMs = np.genfromtxt('Case2_1_ecoli_core/ecoli_core_FBAMs_standardized.csv', delimiter = ',')



# %% case 2.2 ECC2comp















# ecoli_reduced = sio.loadmat('initial_data/ecoli_reduced_model.mat')
# ecoli_reduced_Z = ecoli_reduced['Z']   # rea x pat
# ecoli_reduced_S = ecoli_reduced['S']   # met x rea
# ecoli_reduced_reas_id = ecoli_reduced['reas_id']
# ecoli_reduced_mets_id = ecoli_reduced['mets_id']
#
# yeast_reduced = sio.loadmat('initial_data/yeast_reduced_model.mat')

# model = cobra.io.read_sbml_model('CNA_EColiCore2_compressed.xml')
# model = cobra.io.read_sbml_model('test_core2.xml')
# model = cobra.io.load_matlab_model('ECC2_from_CNA_cobra.mat')