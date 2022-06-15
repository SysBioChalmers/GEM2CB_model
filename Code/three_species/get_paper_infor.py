#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 3/19/22

"""get_paper_infor.py
:description : script
:param : 
:returns: 
:rtype: 
"""



import os

import cobra
import numpy as np
os.chdir('../../Data/three_species/')


name_list = ['B_hyd', 'F_pra', 'R_int']  #

# model_Bhyd = cobra.io.load_json_model('GEM_from_templates_' + name_list[0] + '_refined.json')
# model_Fpra = cobra.io.load_json_model('GEM_from_templates_' + name_list[1] + '_refined.json')
# model_Rint = cobra.io.load_json_model('GEM_from_templates_' + name_list[2] + '_refined.json')
#
# model_Bhyd = cobra.io.load_json_model('GEM_from_templates_' + name_list[0] + '_draft.json')
# model_Fpra = cobra.io.load_json_model('GEM_from_templates_' + name_list[1] + '_draft.json')
# model_Rint = cobra.io.load_json_model('GEM_from_templates_' + name_list[2] + '_draft.json')

model_Bhyd = cobra.io.read_sbml_model('GEM_from_carveme_' + name_list[0] + '_LB.xml')
model_Fpra = cobra.io.read_sbml_model('GEM_from_carveme_' + name_list[1] + '_LB.xml')
model_Rint = cobra.io.read_sbml_model('GEM_from_carveme_' + name_list[2] + '_LB.xml')


#%%
print('model_Bhyd reactions:',len(model_Bhyd.reactions))
print('model_Bhyd metabolites:',len(model_Bhyd.metabolites))
print('model_Fpra reactions:',len(model_Fpra.reactions))
print('model_Fpra metabolites:',len(model_Fpra.metabolites))
print('model_Rint reactions:',len(model_Rint.reactions))
print('model_Rint metabolites:',len(model_Rint.metabolites))

for i in name_list:
    Smz_ = np.genfromtxt(i + '_Smz_constrain.csv', delimiter=',')
    print(i,' shape: ',Smz_.shape)


