#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2/17/20

"""get_GEMs_from_carveme_archive.py
:description : script
:param : 
:returns: 
:rtype: 
"""
import os

import gemstool

os.chdir('../../ComplementaryData/three_species/initial_data')

name_list = ['R_int', 'F_pra', 'B_hyd']
ref_list = ['GCF_000156535.1', 'GCF_000162015.1_ASM16201v1', 'GCF_000157975.1']
for i in name_list:
    gemstool.seq_ana.gbk2faa(i + '.gbff', i + '.faa', 'locus_tag', True)  # sequence:RefSeq: locus_tag

    '''carve --refseq GCF_000005845.2 -o ecoli_k12_mg1655.xml
    carve genome.faa -g M9,LB
    '''
    print('-----  draft models %s  -----' % i)
    # faa_file = i + '.gkbb'
    comd_str = 'carve ' + i + '.faa' + ' -g LB -o ../' + i + 'LB.xml'
    os.system(comd_str)

'''
ref_list = ['GCF_000156535.1','GCF_000162015.1_ASM16201v1','GCF_000157975.1']
for seq_id in ref_list:
    # carve --refseq GCF_000005845.2 -o ecoli_k12_mg1655.xml
    print('-----  draft models %s  -----'%i)
    # faa_file = i + '.gkbb'
    comd_str = 'carve --refseq ' + seq_id + ' -o ../'+ name_list[ref_list.index(seq_id)]+'.xml'
    os.system(comd_str)
'''
