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
import os

os.chdir('../ComplementaryData/')



ecoli_reduced_model = cobra.io.read_sbml_model('ecoli_reduced_model.xml')

model = ecoli_reduced_model.copy()
model.reactions.get_by_id('EX_GLU').bounds = (-10,0.0) #((-10,-0.1)) 0 没有意义啊？？！！
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
model2 = cobra.io.read_sbml_model('e_coli_core.xml')
model2.optimize()
model2.reactions.get_by_id('EX_glc__D_e').bounds = (-10,-0.01)
production_reaid = 'BIOMASS_Ecoli_core_w_GAM'
carbon_rea_id = 'EX_glc__D_e'


def get_yield_m(model2,production_reaid,carbon_rea_id,max_or_min = 'max'):

    model = model2.copy()
    model.objective.direction = max_or_min
    k = 0
    model.objective = {model.reactions.get_by_id(production_reaid):1,model.reactions.get_by_id(carbon_rea_id):k}
    selution = model.optimize()

    while selution.objective_value > 1e-6:
        k = selution.fluxes[production_reaid]/(-selution.fluxes[carbon_rea_id])
        model.objective = {model.reactions.get_by_id(production_reaid):1,model.reactions.get_by_id(carbon_rea_id):k}
        selution = model.optimize()

    fluxes = selution.fluxes
    yield_value = fluxes[production_reaid]/(-fluxes[carbon_rea_id])

    return yield_value,fluxes

yield_value,fluxes = get_yield_m(model2,production_reaid,carbon_rea_id,max_or_min = 'max')




# %%
# model = model2.copy()
# model.objective.direction = 'min'
#
# k = 0
# model.objective = {model.reactions.get_by_id(production_reaid):1,model.reactions.get_by_id(carbon_rea_id):k}
# selution = model.optimize()
#
# while selution.objective_value > 1e-6:
#     k = selution.fluxes[production_reaid]/(-selution.fluxes[carbon_rea_id])
#     model.objective = {model.reactions.get_by_id(production_reaid):1,model.reactions.get_by_id(carbon_rea_id):k}
#     selution = model.optimize()
#
# fluxes_min = selution.fluxes
# yield_value_min = fluxes_min[production_reaid]/(-fluxes_min[carbon_rea_id])




model = model2.copy()
yield_point_flux = model.problem.Constraint(
    model.reactions.get_by_id(production_reaid).flux_expression + 0.08*model.reactions.get_by_id(carbon_rea_id).flux_expression,
    lb=0,
    ub=0)
model.add_cons_vars(yield_point_flux)
selution = model.optimize()
print(selution)





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
