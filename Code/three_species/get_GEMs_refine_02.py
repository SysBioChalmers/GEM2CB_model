#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2/18/20

"""get_GEMs_refine_02.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os

import cobra

os.chdir('../../Data/three_species/')
temp_dir = 'initial_data/templates/'
iNF517 = cobra.io.load_json_model(temp_dir + 'iNF517_standlized.json')
iBT721 = cobra.io.load_json_model(temp_dir + 'iBT721_standlized.json')
iML1515 = cobra.io.load_json_model(temp_dir + 'iML1515_standlized.json')
Lreuteri_530 = cobra.io.load_json_model(temp_dir + 'Lreuteri_530_standlized.json')
name_list = ['B_hyd', 'F_pra', 'R_int']  #
# for sp_name in name_list:
#     # sp_name = 'R_int'
#     model_i_path = 'GEM_from_templates_' + sp_name + '_refined.json'
#     model_i_refined = cobra.io.load_json_model(model_i_path)

Bhyd_initial = cobra.io.load_json_model('GEM_from_templates_' + name_list[0] + '_refined.json')
Fpra_initial = cobra.io.load_json_model('GEM_from_templates_' + name_list[1] + '_refined.json')
Rint_initial = cobra.io.load_json_model('GEM_from_templates_' + name_list[2] + '_refined.json')

# Bhyd_carveme = cobra.io.read_sbml_model('GEM_from_carveme_' + name_list[0] + '_LB.xml')
# Fpra_carveme = cobra.io.read_sbml_model('GEM_from_carveme_' + name_list[1] + '_LB.xml')
Rint_carveme = cobra.io.read_sbml_model('GEM_from_carveme_' + name_list[2] + '_LB.xml')
# %% <R_int>
model = Rint_initial.copy()
model.id = 'R_int'
model.optimize()
model.reactions.get_by_id('PFK').bounds = (0, 1000)
model.optimize()

model.reactions.get_by_id('EX_met__L_e').bounds = (-1, 1000)
model.reactions.get_by_id('EX_his__L_e').bounds = (-1, 1000)
model.reactions.get_by_id('EX_trp__L_e').bounds = (-1, 1000)
model.reactions.get_by_id('EX_tyr__L_e').bounds = (-1, 1000)
model.reactions.get_by_id('EX_gln__L_e').bounds = (-1, 1000)
model.reactions.get_by_id('EX_thr__L_e').bounds = (-1, 1000)
model.reactions.get_by_id('EX_leu__L_e').bounds = (-1, 1000)
model.reactions.get_by_id('EX_arg__L_e').bounds = (-1, 1000)
model.reactions.get_by_id('EX_asn__L_e').bounds = (-1, 1000)
model.reactions.get_by_id('EX_glu__L_e').bounds = (-1, 1000)
model.reactions.get_by_id('EX_ala__L_e').bounds = (-2.69, 1000)
model.reactions.get_by_id('EX_asp__L_e').bounds = (-3.16, 1000)
model.reactions.get_by_id('EX_cys__L_e').bounds = (-0.83, 1000)
model.reactions.get_by_id('EX_gly_e').bounds = (-2.33, 1000)
model.reactions.get_by_id('EX_ile__L_e').bounds = (-1.60, 1000)
model.reactions.get_by_id('EX_lys__L_e').bounds = (-2.68, 1000)
model.reactions.get_by_id('EX_pro__L_e').bounds = (-5.86, 1000)
model.reactions.get_by_id('EX_ser__L_e').bounds = (-3.24, 1000)

# other limations:
model.reactions.get_by_id('EX_etoh_e').bounds = (0, 1000)
model.reactions.get_by_id('EX_lac__L_e').bounds = (0, 1000)
model.reactions.get_by_id('EX_lac__L_e').bounds = (0, 1000)
model.reactions.get_by_id('EX_o2_e').bounds = (0, 1000)

model.reactions.get_by_id('EX_o2_e').bounds = (0, 1000)
model.reactions.get_by_id('EX_glc__D_e').bounds = (0, 1000)
model.reactions.get_by_id('EX_fru_e').bounds = (-10, 1000)
model.add_reaction(Rint_carveme.reactions.get_by_id('FRUpts'))
model.add_reaction(Rint_carveme.reactions.get_by_id('FRUpts2'))
model.optimize()

model.add_reaction(iML1515.reactions.get_by_id('PFL'))
model.reactions.get_by_id('EX_for_e').bounds = (0, 1000)
model.objective = 'EX_for_e'
model.optimize()

model.reactions.get_by_id('EX_ac_e').bounds = (0, 1000)
model.add_reaction(iML1515.reactions.get_by_id('SADT2'))
model.add_reaction(iML1515.reactions.get_by_id('ADSK'))
model.add_reaction(iML1515.reactions.get_by_id('PAPSR'))
model.add_reaction(iML1515.reactions.get_by_id('BPNT'))

model.objective = 'EX_ac_e'
model.optimize()

model.objective = 'EX_etoh_e'
model.optimize()

model.reactions.get_by_id('EX_lac__L_e').bounds = (0, 1000)
model.objective = 'EX_lac__L_e'
model.optimize()

model.add_reaction(Rint_carveme.reactions.get_by_id('ECOAH1'))
model.add_reaction(Rint_carveme.reactions.get_by_id('ACOAD1'))
model.add_reaction(Rint_carveme.reactions.get_by_id('ACOAD1f'))
model.add_reaction(Rint_carveme.reactions.get_by_id('PBUTT'))
model.add_reaction(Rint_carveme.reactions.get_by_id('BUTKr'))
model.add_reaction(Rint_carveme.reactions.get_by_id('BUTt2r'))
model.add_reaction(Rint_carveme.reactions.get_by_id('EX_but_e'))
model.reactions.get_by_id('BUTKr').bounds = (-1000, 1000)
model.objective = 'EX_but_e'
model.optimize()

# Rint_carveme.objective = 'EX_but_e'
# Rint_carveme.optimize()
# gemstool.io.solution2txt(iML1515.optimize(), iML1515, 'test2.txt')
# gemstool.io.gem2txt(Rint_carveme,'Rint_carveme.txt',True)
# iML1515.objective = 'EX_but_e'
# iML1515.optimize()

# for i in model.medium.keys():
#     if model.medium[i] < 1:
#         model.reactions.get_by_id(i).lower_bound = -1


Lreuteri_530.reactions.get_by_id('EX_glc__D_e').bounds = (0, 1000)
Lreuteri_530.reactions.get_by_id('EX_pyr_e').bounds = (0, 1000)
Lreuteri_530.reactions.get_by_id('EX_ac_e').bounds = (0, 1000)
Lreuteri_530.optimize()
# gemstool.io.solution2txt(Lreuteri_530.optimize(), Lreuteri_530, 'test2.txt')

model.reactions.get_by_id('EX_pyr_e').bounds = (0, 1000)
model.reactions.get_by_id('EX_fru_e').bounds = (0, 1000)

model.add_reaction(iML1515.reactions.get_by_id('CS'))
model.add_reaction(iML1515.reactions.get_by_id('ACONTb'))
model.add_reaction(iML1515.reactions.get_by_id('ACONTa'))
model.add_reaction(iML1515.reactions.get_by_id('POR5'))
model.add_reaction(iML1515.reactions.get_by_id('FLDR2'))
model.add_reaction(Lreuteri_530.reactions.get_by_id('GLTAL'))
model.add_reaction(Lreuteri_530.reactions.get_by_id('DALTAL'))
model.add_reaction(Lreuteri_530.reactions.get_by_id('kaasIII'))
model.add_reaction(Lreuteri_530.reactions.get_by_id('BTMAT1'))

model.objective = 'BIOMASS'
model.optimize()
# gemstool.io.solution2txt(model.optimize(), model, 'test2.txt')


# # %%

# gemstool.io.solution2txt(iML1515.optimize(), iML1515, 'test2.txt')
# # %%
# Lreuteri_530.objective = 'BIOMASS'
# Lreuteri_530.reactions.get_by_id('EX_ac_e').bounds = (0, 1000)
# Lreuteri_530.reactions.get_by_id('EX_glc__D_e').bounds = (-10, 1000)
# Lreuteri_530.reactions.get_by_id('EX_fru_e').bounds = (0, 1000)
# f1 = Lreuteri_530.optimize()
# Lreuteri_530.reactions.get_by_id('EX_ac_e').bounds = (-0, 1000)
# Lreuteri_530.reactions.get_by_id('EX_glc__D_e').bounds = (0, 1000)
# Lreuteri_530.reactions.get_by_id('EX_fru_e').bounds = (0, 1000)
# Lreuteri_530.reactions.get_by_id('EX_pyr_e').bounds = (-10, 1000)
# f2 = Lreuteri_530.optimize()
# gemstool.io.solution2txt(f2, Lreuteri_530, 'test3.txt')
#
#
# # %%
# model.reactions.get_by_id('EX_nh4_e').bounds = (-1000, 1000)
# model.add_reaction(iML1515.reactions.get_by_id('EX_nh4_e'))
# model.add_reaction(iML1515.reactions.get_by_id('AKGDH'))
# model.add_reaction(iML1515.reactions.get_by_id('ATPS4rpp'))
# model.add_reaction(iML1515.reactions.get_by_id('SUCDi'))
# model.add_reaction(iML1515.reactions.get_by_id('SUCOAS'))
#
#
# model.objective = 'BIOMASS'
# model.reactions.get_by_id('EX_ac_e').bounds = (0, 1000)
# model.reactions.get_by_id('EX_glc__D_e').bounds = (-10, 1000)
# model.reactions.get_by_id('EX_fru_e').bounds = (0, 1000)
# f1 = model.optimize()
# model.reactions.get_by_id('EX_ac_e').bounds = (-0, 1000)
# model.reactions.get_by_id('EX_glc__D_e').bounds = (0, 1000)
# model.reactions.get_by_id('EX_fru_e').bounds = (0, 1000)
# model.reactions.get_by_id('EX_pyr_e').bounds = (-10, 1000)
# f2 = model.optimize()
# gemstool.io.solution2txt(f2, model, 'test3.txt')
#
# rea_list = [i.id for i in model.reactions]
# set(f2[abs(f2.fluxes)>0.1].index) - set(rea_list)
# gemstool.io.solution2txt(iML1515.optimize(), iML1515, 'test2.txt')
# model.optimize()


# Growth = Rint_carveme.reactions.get_by_id('Growth')
# model.add_reaction(Growth)
# model.objective = Growth
# model.optimize()
iML1515.reactions.get_by_id('EX_ac_e').bounds = (0, 1000)
iML1515.reactions.get_by_id('EX_glc__D_e').bounds = (-10, 1000)
f1 = iML1515.optimize()
iML1515.reactions.get_by_id('EX_ac_e').bounds = (-10, 1000)
iML1515.reactions.get_by_id('EX_glc__D_e').bounds = (0, 1000)
f2 = iML1515.optimize()
rate_glc_ac = f1.objective_value / f2.objective_value
rate_glc_ac = 4.17
model.objective = 'BIOMASS'
model.reactions.get_by_id('EX_ac_e').bounds = (0, 1000)
model.reactions.get_by_id('EX_glc__D_e').bounds = (-10, 1000)
model.reactions.get_by_id('EX_fru_e').bounds = (0, 1000)
model.reactions.get_by_id('EX_pyr_e').bounds = (0, 1000)
model.optimize()

Suppose_rea_ac = cobra.Reaction('Suppose_rea_ac')
model.add_reaction(Suppose_rea_ac)
model.reactions.get_by_id('Suppose_rea_ac').reaction = '4.17 ac_e --> glc__D_e'
model.optimize()

cobra.io.save_json_model(model, 'GEM_from_templates_' + name_list[2] + '_refined.json', sort='True')

# %% <F_pra>
model = Fpra_initial.copy()
model.optimize()
model.reactions.get_by_id('PFK').bounds = (0, 1000)
model.optimize()

model.reactions.get_by_id('EX_met__L_e').bounds = (-1, 1000)
model.reactions.get_by_id('EX_his__L_e').bounds = (-1, 1000)
model.reactions.get_by_id('EX_trp__L_e').bounds = (-1, 1000)
model.reactions.get_by_id('EX_tyr__L_e').bounds = (-1, 1000)
model.reactions.get_by_id('EX_gln__L_e').bounds = (-1, 1000)
model.reactions.get_by_id('EX_thr__L_e').bounds = (-1, 1000)
model.reactions.get_by_id('EX_leu__L_e').bounds = (-1, 1000)
model.reactions.get_by_id('EX_arg__L_e').bounds = (-1, 1000)
model.reactions.get_by_id('EX_asn__L_e').bounds = (-1, 1000)
model.reactions.get_by_id('EX_glu__L_e').bounds = (-1, 1000)
model.reactions.get_by_id('EX_ala__L_e').bounds = (-2.69, 1000)
model.reactions.get_by_id('EX_asp__L_e').bounds = (-3.16, 1000)
model.reactions.get_by_id('EX_cys__L_e').bounds = (-0.83, 1000)
model.reactions.get_by_id('EX_gly_e').bounds = (-2.33, 1000)
model.reactions.get_by_id('EX_ile__L_e').bounds = (-1.60, 1000)
model.reactions.get_by_id('EX_lys__L_e').bounds = (-2.68, 1000)
model.reactions.get_by_id('EX_pro__L_e').bounds = (-5.86, 1000)
model.reactions.get_by_id('EX_ser__L_e').bounds = (-3.24, 1000)

model.reactions.get_by_id('EX_o2_e').bounds = (0, 1000)
model.reactions.get_by_id('EX_glc__D_e').bounds = (0, 1000)
model.reactions.get_by_id('EX_fru_e').bounds = (-10, 1000)
model.add_reaction(Rint_carveme.reactions.get_by_id('FRUpts'))
model.add_reaction(Rint_carveme.reactions.get_by_id('FRUpts2'))
model.optimize()

model.add_reaction(iML1515.reactions.get_by_id('PFL'))
model.reactions.get_by_id('EX_for_e').bounds = (0, 1000)
model.objective = 'EX_for_e'
model.optimize()

model.reactions.get_by_id('EX_ac_e').bounds = (0, 1000)
model.add_reaction(iML1515.reactions.get_by_id('SADT2'))
model.add_reaction(iML1515.reactions.get_by_id('ADSK'))
model.add_reaction(iML1515.reactions.get_by_id('PAPSR'))
model.objective = 'EX_ac_e'
model.optimize()

model.objective = 'EX_etoh_e'
model.optimize()

model.reactions.get_by_id('EX_lac__L_e').bounds = (0, 1000)
model.objective = 'EX_lac__L_e'
model.optimize()

model.add_reaction(Rint_carveme.reactions.get_by_id('ECOAH1'))
model.add_reaction(Rint_carveme.reactions.get_by_id('ACOAD1'))
model.add_reaction(Rint_carveme.reactions.get_by_id('ACOAD1f'))
model.add_reaction(Rint_carveme.reactions.get_by_id('PBUTT'))
model.add_reaction(Rint_carveme.reactions.get_by_id('BUTKr'))
model.add_reaction(Rint_carveme.reactions.get_by_id('BUTt2r'))
model.add_reaction(Rint_carveme.reactions.get_by_id('EX_but_e'))
model.reactions.get_by_id('BUTKr').bounds = (-1000, 1000)
model.objective = 'EX_but_e'
model.optimize()

model.objective = 'BIOMASS'
model.reactions.get_by_id('EX_ac_e').bounds = (0, 1000)
model.reactions.get_by_id('EX_glc__D_e').bounds = (-10, 1000)
model.reactions.get_by_id('EX_fru_e').bounds = (0, 1000)
model.reactions.get_by_id('EX_pyr_e').bounds = (0, 1000)
model.optimize()

Suppose_rea_ac = cobra.Reaction('Suppose_rea_ac')
model.add_reaction(Suppose_rea_ac)
model.reactions.get_by_id('Suppose_rea_ac').reaction = '4.17 ac_e --> glc__D_e'
model.optimize()

cobra.io.save_json_model(model, 'GEM_from_templates_' + name_list[1] + '_refined.json', sort='True')

# %% <B_hyd>
model = Bhyd_initial.copy()
model.optimize()
model.reactions.get_by_id('PFK').bounds = (0, 1000)
model.optimize()

model.reactions.get_by_id('EX_met__L_e').bounds = (-1, 1000)
model.reactions.get_by_id('EX_his__L_e').bounds = (-1, 1000)
model.reactions.get_by_id('EX_trp__L_e').bounds = (-1, 1000)
model.reactions.get_by_id('EX_tyr__L_e').bounds = (-1, 1000)
model.reactions.get_by_id('EX_gln__L_e').bounds = (-1, 1000)
model.reactions.get_by_id('EX_thr__L_e').bounds = (-1, 1000)
model.reactions.get_by_id('EX_leu__L_e').bounds = (-1, 1000)
model.reactions.get_by_id('EX_arg__L_e').bounds = (-1, 1000)
model.reactions.get_by_id('EX_asn__L_e').bounds = (-1, 1000)
model.reactions.get_by_id('EX_glu__L_e').bounds = (-1, 1000)
model.reactions.get_by_id('EX_ala__L_e').bounds = (-2.69, 1000)
model.reactions.get_by_id('EX_asp__L_e').bounds = (-3.16, 1000)
model.reactions.get_by_id('EX_cys__L_e').bounds = (-0.83, 1000)
model.reactions.get_by_id('EX_gly_e').bounds = (-2.33, 1000)
model.reactions.get_by_id('EX_ile__L_e').bounds = (-1.60, 1000)
model.reactions.get_by_id('EX_lys__L_e').bounds = (-2.68, 1000)
model.reactions.get_by_id('EX_pro__L_e').bounds = (-5.86, 1000)
model.reactions.get_by_id('EX_ser__L_e').bounds = (-3.24, 1000)

model.reactions.get_by_id('EX_o2_e').bounds = (0, 1000)
model.reactions.get_by_id('EX_glc__D_e').bounds = (0, 1000)
model.reactions.get_by_id('EX_fru_e').bounds = (-10, 1000)
model.add_reaction(Rint_carveme.reactions.get_by_id('FRUpts'))
model.add_reaction(Rint_carveme.reactions.get_by_id('FRUpts2'))
model.optimize()

model.add_reaction(iML1515.reactions.get_by_id('PFL'))
model.reactions.get_by_id('EX_for_e').bounds = (0, 1000)
model.objective = 'EX_for_e'
model.optimize()

model.reactions.get_by_id('EX_ac_e').bounds = (0, 1000)
model.add_reaction(iML1515.reactions.get_by_id('SADT2'))
model.add_reaction(iML1515.reactions.get_by_id('ADSK'))
model.add_reaction(iML1515.reactions.get_by_id('PAPSR'))
model.objective = 'EX_ac_e'
model.optimize()

model.objective = 'EX_etoh_e'
model.optimize()

model.reactions.get_by_id('EX_lac__L_e').bounds = (0, 1000)
model.objective = 'EX_lac__L_e'
model.optimize()

# model.add_reaction(Rint_carveme.reactions.get_by_id('ECOAH1'))
# model.add_reaction(Rint_carveme.reactions.get_by_id('ACOAD1'))
# model.add_reaction(Rint_carveme.reactions.get_by_id('ACOAD1f'))
# model.add_reaction(Rint_carveme.reactions.get_by_id('PBUTT'))
# model.add_reaction(Rint_carveme.reactions.get_by_id('BUTKr'))
# model.add_reaction(Rint_carveme.reactions.get_by_id('BUTt2r'))
model.add_reaction(Rint_carveme.reactions.get_by_id('EX_but_e'))
# model.reactions.get_by_id('BUTKr').bounds = (-1000, 1000)
model.objective = 'EX_but_e'
model.optimize()

model.reactions.get_by_id('EX_fru_e').bounds = (0, 1000)
model.reactions.get_by_id('EX_for_e').bounds = (-10, 1000)
model.objective = 'BIOMASS'
print(model.optimize())

cobra.io.save_json_model(model, 'GEM_from_templates_' + name_list[0] + '_refined.json', sort='True')
