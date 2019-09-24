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





iML1515 = cobra.io.read_sbml_model('../ComplementaryData/iML1515.xml')




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



yiled_rea_ids = ['EX_ac_e','EX_for_e','EX_succ_e','EX_lac__D_e','EX_etoh_e']
yiled_rea = [iML1515.reactions.get_by_id(i) for i in yiled_rea_ids]
# %%

b = pd.DataFrame()
for i in range(1,10):
    iML1515.reactions.get_by_id('BIOMASS_Ec_iML1515_core_75p37M').bounds = (solution.objective_value/i,solution.objective_value/i)
    a = flux_variability_analysis(iML1515,yiled_rea)
    a['biomass'] = solution.objective_value/i
    b = b.append(a)

c = b.loc[ 'EX_etoh_e' , : ].plot(kind='line', x='biomass', y='maximum')


# c.plot(
#      kind='line', x='biomass', y='maximum')

plt.show()

# %%







