#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 9/26/19

"""fva_test.py
:description : script
:param : 
:returns: 
:rtype: 
"""
import cobra
from pandas import DataFrame
from optlang.symbolics import Zero
from numpy import zeros



ecoli_reduced_model = cobra.io.read_sbml_model('ecoli_reduced_model.xml')

model = ecoli_reduced_model.copy()
model.reactions.get_by_id('EX_GLU').bounds = (-10,-0.1)
yield_rea_ids = ['EX_ACT','EX_FOR','EX_ETH','EX_LAC','EX_SUC']
a = model.optimize().objective_value


model.reactions.get_by_id('EX_B').bounds = (a,a)

for rea_id in yield_rea_ids:
    model.objective = rea_id
    model.objective.direction = 'max'
    print(model.optimize())
    model.objective.direction = 'min'
    print(model.optimize())

def fva_self_def(model1, reaction_list=None, pfba=False,):
    model = model1.copy()

    if reaction_list is None:
        reaction_ids = [r.id for r in model.reactions]
    else:
        reaction_ids = reaction_list

    fva_value_result = DataFrame([])
    fva_fluxes_result = DataFrame([])

    for rea_id in reaction_ids:         #TODO pfba
        model.objective = rea_id
        model.objective.direction = 'max'

        f_max = model.optimize()
        value_max = f_max.objective_value
        fluxes_max = f_max.fluxes

        model.objective.direction = 'min'
        f_min = model.optimize()
        value_min = f_min.objective_value
        fluxes_min = f_min.fluxes


        fva_value_result.at[rea_id, 'max'] = value_max
        fva_fluxes_result[rea_id + 'max'] = fluxes_max
        fva_value_result.at[rea_id, 'min'] = value_min
        fva_fluxes_result[rea_id + 'min'] = fluxes_min

    return fva_value_result,fva_fluxes_result


fva_value_result,fva_fluxes_result = fva_self_def(model, reaction_list=yield_rea_ids, pfba=False,)

# %%

import matplotlib.pyplot as plt
import numpy as np
fig = plt.figure()
ax1 = fig.add_subplot(111)
np.random()

ax1.plot([0,1,], [],'x',markerfacecolor='none',color = 'black',label = 'EFM',markersize=12)
# ax1.plot(FBAem_all_points[:,0], FBAem_all_points[:,index],'o',markerfacecolor='none',color = 'tab:red',alpha = 0.8,label = 'FBA_mode',markersize=12)
# ax1.plot(our_all_points[:,0], our_all_points[:,index],'v',markerfacecolor='none',color = 'tab:blue',alpha = 0.8,label = 'This_methd',markersize=12)
ax1.legend()
ax1.set_xlabel('Yield Biomass/Glucose',fontsize = 12)
ax1.set_ylabel('Yield Acetate/Glucose',fontsize = 12)
fig.show()
