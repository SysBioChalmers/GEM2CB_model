#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 9/27/19

"""step0_process_initial_z.py
:description : script
:param : 
:returns: 
:rtype: 
"""
import scipy.io as sio
import cobra
import os

os.chdir('../ComplementaryData/')

# ecoli_reduced = sio.loadmat('initial_data/ecoli_reduced_model.mat')
# ecoli_reduced_Z = ecoli_reduced['Z']   # rea x pat
# ecoli_reduced_S = ecoli_reduced['S']   # met x rea
# ecoli_reduced_reas_id = ecoli_reduced['reas_id']
# ecoli_reduced_mets_id = ecoli_reduced['mets_id']
#
# yeast_reduced = sio.loadmat('initial_data/yeast_reduced_model.mat')

