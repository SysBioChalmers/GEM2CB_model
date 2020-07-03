#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2/18/20

"""get_GEMs_refine.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os

import cobra
import gemstool

os.chdir('../../ComplementaryData/three_species/')
iML1515 = cobra.io.read_sbml_model('../Case3_iML1515/iML1515.xml')
Ecoli_core = cobra.io.read_sbml_model('../Case2_1_ecoli_core/e_coli_core.xml')
Ecoli_core.reactions.get_by_id('EX_glc__D_e').bounds = (0, 1000)
Ecoli_core.reactions.get_by_id('EX_ac_e').bounds = (-10, 1000)
Ecoli_core.optimize()

iML1515.objective = 'EX_4abut_e'
iML1515.optimize()
gemstool.io.solution2txt(iML1515.optimize(), iML1515, 'test2.txt')

Rint_initial = cobra.io.load_json_model('R_intLreu_draft_3_refined_20-01-21.json')
Bhyd_initial = cobra.io.load_json_model('B_hydLreu_draft_3_refined_20-01-21.json')
Fpra_initial = cobra.io.load_json_model('F_praLreu_draft_3_refined_20-01-21.json')

# %% Rint
model = Rint_initial.copy()
model.optimize()

model.reactions.get_by_id('ICDHyr').bounds = (-1000, 1000)
model.reactions.get_by_id('ACONT').bounds = (-1000, 1000)
model.reactions.get_by_id('O2t').bounds = (-1000, 1000)
model.reactions.get_by_id('EX_o2_e').bounds = (-1000, 1000)
model.optimize()

for i in model.medium.keys():
    if model.medium[i] < 1:
        model.reactions.get_by_id(i).lower_bound = -1

model.optimize()
biomass = cobra.Reaction('biomass')
model.add_reaction(biomass)
model.reactions.get_by_id('biomass').reaction = Ecoli_core.reactions.get_by_id('BIOMASS_Ecoli_core_w_GAM').reaction
model.objective = 'biomass'
model.optimize()
# uptake fru
model.reactions.get_by_id('EX_glc__D_e').bounds = (0, 1000)
model.reactions.get_by_id('EX_fru_e').bounds = (-10, 1000)
model.optimize()

# uptake ac
model.reactions.get_by_id('EX_fru_e').bounds = (0, 1000)
model.reactions.get_by_id('EX_ac_e').bounds = (-10, 1000)
model.optimize()

# prduce for
model.objective = 'EX_for_e'
model.optimize()

# prduce but
model.objective = 'EX_4abut_e'
model.optimize()

cobra.io.save_json_model(model, 'Rint_growth.json')

# %% Fpra
model = Fpra_initial.copy()
model.optimize()

model.reactions.get_by_id('ICDHyr').bounds = (-1000, 1000)
model.reactions.get_by_id('ACONT').bounds = (-1000, 1000)
model.reactions.get_by_id('O2t').bounds = (-1000, 1000)
model.reactions.get_by_id('EX_o2_e').bounds = (-1000, 1000)
for i in model.medium.keys():
    if model.medium[i] < 1:
        model.reactions.get_by_id(i).lower_bound = -1
model.optimize()

biomass = cobra.Reaction('biomass')
model.add_reaction(biomass)
model.reactions.get_by_id('biomass').reaction = Ecoli_core.reactions.get_by_id('BIOMASS_Ecoli_core_w_GAM').reaction
model.objective = 'biomass'
model.optimize()

# uptake fru
model.reactions.get_by_id('EX_glc__D_e').bounds = (0, 1000)
model.reactions.get_by_id('EX_fru_e').bounds = (-10, 1000)
model.optimize()

# uptake ac
model.reactions.get_by_id('EX_fru_e').bounds = (0, 1000)
model.reactions.get_by_id('EX_ac_e').bounds = (-10, 1000)

model.optimize()

# prduce for
model.objective = 'EX_for_e'
model.optimize()

# prduce but
model.objective = 'EX_4abut_e'
model.optimize()

cobra.io.save_json_model(model, 'Fpra_growth.json')

# %% BH
model = Bhyd_initial.copy()
model.optimize()
model.reactions.get_by_id('ICDHyr').bounds = (-1000, 1000)
model.reactions.get_by_id('ACONT').bounds = (-1000, 1000)
model.reactions.get_by_id('O2t').bounds = (-1000, 1000)
model.reactions.get_by_id('EX_o2_e').bounds = (-1000, 1000)
for i in model.medium.keys():
    if model.medium[i] < 1:
        model.reactions.get_by_id(i).lower_bound = -1
model.optimize()

biomass = cobra.Reaction('biomass')
model.add_reaction(biomass)
model.reactions.get_by_id('biomass').reaction = Ecoli_core.reactions.get_by_id('BIOMASS_Ecoli_core_w_GAM').reaction
model.objective = 'biomass'
model.optimize()

# uptake fru
model.reactions.get_by_id('EX_glc__D_e').bounds = (0, 1000)
model.reactions.get_by_id('EX_fru_e').bounds = (-10, 1000)
model.optimize()

# uptake for
model.reactions.get_by_id('EX_fru_e').bounds = (0, 1000)
model.reactions.get_by_id('EX_ac_e').bounds = (0, 1000)
model.reactions.get_by_id('EX_for_e').bounds = (-10, 1000)

# model.reactions.get_by_id('ATPM').bounds = (0,1000)
model.objective = 'biomass'
model.optimize()

# prduce ac
model.objective = 'EX_ac_e'
model.optimize()

# prduce lac
model.objective = 'EX_lac__L_e'
model.optimize()

cobra.io.save_json_model(model, 'Bhyd_growth.json')

# gemstool.io.gem2txt(model,'Branch_work.txt')
# gemstool.io.solution2txt(model.optimize(),model,'test2.txt')
