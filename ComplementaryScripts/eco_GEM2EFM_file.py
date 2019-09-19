#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-09-16

"""eco_GEM2EFM_file.py
:description : script to convert model to emtools files
:param : 
:returns: 
:rtype: 
"""

import cobra
import os
import My_def


# def GEM2CB(model,outfile,sort=False):
#
#     file = open(outfile,'w')
#
#     file.write('-ENZREV\n\n')
#
#     file.write('-ENZIRREV\n\n')
#
#     file.write('-METEXT\n\n')
#
#     METEXT = []
#     for met in model.metabolites:
#         if '_e' in met.id:
#             METEXT.append(met.id)
#
#     file.write(' '.join(METEXT))
#
#     file.write('\n-CAT\n')
#
#     for rea in model.reactions:
#         file.write(rea.id +' : ')
#         equ = rea.reaction
#         equ = equ.replace('-->','=>')
#         equ = equ.replace('<=>','=')
#         equ = equ.replace('<--','<=')
#         file.write(equ +' .\n')
#     file.close()



os.chdir('../ComplementaryData/')




# iML1515 = cobra.io.read_sbml_model('iML1515.xml')
#
# My_def.io_file.GEM2CB(iML1515,'iML1515_em.txt',sort=False)
#
# os.system('cp iML1515_em.txt ~/Documents/MATLAB/aumic-master/GEM2CB/iML1515_em.txt')


e_coli_core = cobra.io.read_sbml_model('e_coli_core.xml')

My_def.io_file.GEM2CB(e_coli_core,'e_coli_core_em.txt',sort=False)

os.system('cp e_coli_core_em.txt ~/Documents/MATLAB/aumic-master/GEM2CB/e_coli_core_em.txt')

