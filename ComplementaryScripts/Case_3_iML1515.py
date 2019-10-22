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
yield_normalized_df = GEM2pathways.get_yield_space(model, biomass_rea_id,carbon_rea_id,yield_rea_ids,step_of_biomass,carbon_uptake_direction = carbon_uptake_direction,draw = True)

#MYA hull#MYA experiment data
yield_normalized_df = yield_normalized_df.loc[[carbon_rea_id] + [biomass_rea_id] + yield_rea_ids ,:]

our_all_points = yield_normalized_df.values.T
our_all_points = our_all_points[ abs(our_all_points[:,0]) > 1e-10,:]
our_all_points = our_all_points[:,1:]       #exclude the carbon colume

experiment_datas = [ ]      # TODO experiment data!!!
qhull_options = 'Qt QJ Pp Qw Qx'      #'QG0'mean expect point 0(index)
cutoff_persent = 0.99







# %%






