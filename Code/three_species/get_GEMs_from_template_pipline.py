#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by lhao at 2019-02-27


"""pipeline.py
:description : script to renconstract GEM for L reuteri
:param : templates models and seqs (standardizated) iNF517,iML1515,iBT721
:returns: draft model : Lreu_draft_3_refined
:rtype:
"""

import os

import cobra
import gemstool
import pandas as pd

import My_def

# from Bio import SeqIO

os.chdir('../../Data/three_species/')

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)


def get_model_from_tp(qseq_file, sseq_file, tp_model):
    # <step blast(used dimond) >
    blast_result_s_in_q, blast_result_q_in_s = gemstool.seq_ana.blastp_pairwise(qseq_file, sseq_file, out_dir='blast/')
    # <step select blast results: cut off bitscore=100, ppos=45 best(BBHs)>
    blast_result_df = gemstool.seq_ana.select_blast(blast_result_s_in_q, blast_result_q_in_s, best_match=True,
                                                    evalue=10 ** -10, pident=0, length=0,
                                                    bitscore=100, ppos=45, qcovs=0)
    # blast_result_df.to_csv('result.csv')

    os.system('rm ' + blast_result_s_in_q)
    os.system('rm ' + blast_result_q_in_s)
    # <step get draft model from templates>
    # note: only extracted the matched reactions
    model_i_tp = My_def.get_draft_from_template(tp_model, blast_result_df, remove_missing_genes=False)
    # model_i_tp = gemstool.get_draft_from_template(tp_model, blast_result_df,remove_missing_genes = False)

    return model_i_tp


def add_basic_rea(model_i, main_template):
    reaset = set([i.id for i in model_i.reactions])
    metset = set([i.id for i in model_i.metabolites])
    for rea in main_template.reactions:
        if ('EX_' in rea.id):
            # case exchange reactions
            if rea.id not in reaset:
                rea.notes['from'] = [main_template.id, 'exchange']
                model_i.add_reaction(rea)
                reaset.add(rea.id)
        elif ('_c' in rea.reaction) and ('_e' in rea.reaction):
            # case transport
            if rea.id not in reaset:
                rea.notes['from'] = [main_template.id, 'transport']
                model_i.add_reaction(rea)
                reaset.add(rea.id)

        elif ('_LRE' in rea.id) or ('LRE_c' in rea.reaction):
            # case biomass
            # print(rea)
            if rea.functional:
                if rea.id not in reaset:
                    rea.notes['from'] = [main_template.id, 'biomass']
                    model_i.add_reaction(rea)
                    reaset.add(rea.id)
        elif rea.id in ['ATPM'] and rea.id not in reaset:
            # case ATPM
            rea.notes['from'] = [main_template.id, 'atp']
            main_template.add_reaction(rea)
            reaset.add(rea.id)
    return main_template


# %% <general step: load data: seqs, templates models >

print('----- loading data -----')
os.chdir('../../Data/three_species/')

name_list = ['B_hyd', 'F_pra', 'R_int']

# for i in name_list:
#     gbk_file = i + '.gbff'
#     faa_file = i + '.faa'
#     gemstool.seq_ana.gbk2faa(gbk_file,faa_file,locus_tag = 'locus_tag')

# counts = []
# for i in name_list:
#
#     faa_file = i + '.faa'
#     n = 0
#     for seq in SeqIO.parse(faa_file, "fasta"):
#         n = n+1
#         continue
#     counts.append(n)

print('===== Template models and seqs processed, Done =====')

# %%

for sp_name in name_list:
    os.chdir('initial_data/templates/')

    boj_model_seq = '../' + sp_name + '.faa'

    iNF517_seq = 'iNF517.faa'
    iBT721_seq = 'iBT721.faa'
    iML1515_seq = 'iML1515.faa'
    Lreuteri_530_seq = 'Lreuteri_530.faa'

    iNF517 = cobra.io.load_json_model('iNF517_standlized.json')
    iBT721 = cobra.io.load_json_model('iBT721_standlized.json')
    iML1515 = cobra.io.load_json_model('iML1515_standlized.json')
    Lreuteri_530 = cobra.io.load_json_model('Lreuteri_530_standlized.json')

    # %% <general step: blast and get draft models. steps details could be find in in get_Lreu_from_tp function (above)>

    print('----- blasting and geting draft models  -----')
    model_i_from_iNF517 = get_model_from_tp(boj_model_seq, iNF517_seq, iNF517)
    model_i_from_iBT721 = get_model_from_tp(boj_model_seq, iBT721_seq, iBT721)
    model_i_from_iML1515 = get_model_from_tp(boj_model_seq, iML1515_seq, iML1515)
    model_i_from_Lreuteri_530 = get_model_from_tp(boj_model_seq, Lreuteri_530_seq, Lreuteri_530)

    # %% <special step: process model_i_from_iML1515 (have periplasm(_p) compartment)>

    model_i_from_iML1515 = gemstool.model_refine.remove_compartment(model_i_from_iML1515, compartment='_p')

    # %% <general step: add reaction 'from' notes and save models>

    model_i_from_iNF517.id = 'model_i_from_iNF517'
    model_i_from_iBT721.id = 'model_i_from_iBT721'
    model_i_from_iML1515.id = 'model_i_from_iML1515'
    model_i_from_Lreuteri_530.id = 'model_i_from_iML1515'

    gemstool.merge_model.note_model_from(model_i_from_iNF517, ['iNF517', 'BBH'])
    gemstool.merge_model.note_model_from(model_i_from_iBT721, ['iBT721', 'BBH'])
    gemstool.merge_model.note_model_from(model_i_from_iML1515, ['iML1515', 'BBH'])
    gemstool.merge_model.note_model_from(model_i_from_Lreuteri_530, ['Lreuteri_530', 'BBH'])

    # cobra.io.save_json_model(model_i_from_iNF517,'model_i_from_iNF517.json',sort=True)
    # cobra.io.save_json_model(model_i_from_iBT721,'model_i_from_iBT721.json',sort=True)
    # cobra.io.save_json_model(model_i_from_iML1515,'model_i_from_iML1515.json',sort=True)
    # cobra.io.save_json_model(model_i_from_Lreuteri_530,'model_i_from_Lreuteri_530.json',sort=True)

    # %% <general step: select a main model, iNF517>

    os.chdir('../../')
    print('----- Mergeing draft models  -----')
    model_i_draft_1 = model_i_from_Lreuteri_530.copy()

    # %% <general step: merge models, add unique reactions from other templates, iML1515 and iBT721>

    #   option 1
    #   by cobra function
    # filed
    # Lreu = Lreu.merge(model_i_from_iML1515,inplace=False)
    # Lreu = Lreu.merge(model_i_from_iBT721,inplace=False)

    #   option 2
    # by def function

    model_i_draft_2, report_df_from_iML1515 = gemstool.merge_model.merge_draftmodels(model_i_draft_1,
                                                                                     model_i_from_iML1515)
    model_i_draft_2, report_df_from_1BT721 = gemstool.merge_model.merge_draftmodels(model_i_draft_2,
                                                                                    model_i_from_iNF517)
    model_i_draft_2, report_df_from_Lreuteri_530 = gemstool.merge_model.merge_draftmodels(model_i_draft_2,
                                                                                          model_i_from_iBT721)

    print('\033[0;34;48m')

    # %% <general step: add exchange transport and biomass reactions >
    # code not general, depends on templates model ex tran bio reaction ids

    print('----- add exchange transport and biomass reactions  -----')
    model_i_draft_2 = add_basic_rea(model_i_draft_2, Lreuteri_530)
    # cobra.io.save_json_model(model_i_draft_2,'model_i_draft_2.json',sort='True')

    # %% <general step: check draft model FBA results >

    print('----- FBA result  -----')
    model_i_draft_3 = model_i_draft_2.copy()
    Lreuteri_530.solver = 'cplex'
    model_i_draft_3.solver = 'cplex'
    # Bug: solver status is 'infeasible'

    Lreuteri_530.objective = "BIOMASS"
    print('Lreuteri_530 Biomass:', Lreuteri_530.optimize())

    model_i_draft_3.objective = "BIOMASS"
    print('model_i_draft_3 Biomass:', model_i_draft_3.optimize())

    # %% <special step: reless release constraint avoide infeasible fab results (remove midum update rates )>

    for rea in model_i_draft_3.reactions:
        if 'EX' in rea.id:
            if rea.lower_bound <= 0 and rea.upper_bound <= 0:
                rea.upper_bound = 0.0
            elif rea.lower_bound >= 0 and rea.upper_bound >= 0:
                rea.lower_bound = 0.0

    for rea in Lreuteri_530.reactions:
        if 'EX' in rea.id:
            if rea.lower_bound <= 0 and rea.upper_bound <= 0:
                rea.upper_bound = 0.0
            elif rea.lower_bound >= 0 and rea.upper_bound >= 0:
                rea.lower_bound = 0.0

    # %% <general step: gap fill >

    print('----- Gap filling  -----')
    # note, if still failed , need change other way to gap fill, such as biomass partly gapfill, check templates FVA even FBA results.
    # solution_biomass = cobra.flux_analysis.gapfill(model_i_draft_2, iNF517)
    # filed
    # note: need to set integer_threshold=1e-10
    print('gap filling')
    # gap_biomass_mets = ['LIP_LRE_c','LTAtotal_LRE_c','adeadocbl_c','btn_c','pydx5p_c']
    gap_biomass_mets = [i.id for i, j in model_i_draft_3.reactions.get_by_id('BIOMASS').metabolites.items() if j < 0]
    gap_biomass_mets = set(gap_biomass_mets) - {'atp_c', 'coa_c', 'h2o_c', 'nad_c', }
    # solution_biomass = cobra.flux_analysis.gapfill(model_i_draft_3, Lreuteri_530)
    if model_i_draft_3.optimize().objective_value == 0.0:

        try:

            solution_biomass_f = cobra.flux_analysis.gapfilling.GapFiller(model_i_draft_3, Lreuteri_530,
                                                                          demand_reactions=False,
                                                                          integer_threshold=1e-10)
            solution_biomass = solution_biomass_f.fill(iterations=1)[0]
            biomass_gaps_set = set([i.id for i in solution_biomass])
            print('biomass gaps number:', len(biomass_gaps_set))

            reaset = set([i.id for i in model_i_draft_3.reactions])
            metset = set([i.id for i in model_i_draft_3.metabolites])

            for i in biomass_gaps_set:
                rea = Lreuteri_530.reactions.get_by_id(i)
                rea.notes['from'] = [Lreuteri_530.id, 'gap']
                # Bug !!!
                model_i_draft_3.add_reaction(rea)
        except:
            print('gapfill failed')

            # gap_partly_set = set(['PYDAMt', 'ALATA_Lr', 'AOBUTDs' ,'BTNt2i','FA161tr','FA182tr','FA183tr','AACPS183'])
            gaps_set = set([])

            # gap fill unstable!!! here are the gap if the function run well
            try:
                model_template = Lreuteri_530.copy()
                model_gap = model_i_draft_3.copy()
                model_template.objective = "BIOMASS"
                # print('Lreuteri_530:', model_template.optimize())
                model_gap.objective = "BIOMASS"
                # print('model_i_draft_3_refined:',model_gap.optimize())

                if model_gap.optimize().objective_value < 1e-10:
                    solution_part_biomass = cobra.flux_analysis.gapfill(model_gap, model_template,
                                                                        demand_reactions=False)
                    gaps_set = gaps_set | set([i.id for i in solution_part_biomass[0]])
            except:
                try:
                    for gap_biomass_met in gap_biomass_mets:
                        model_template = Lreuteri_530.copy()
                        model_gap = model_i_draft_3.copy()

                        rea_temp = cobra.Reaction('object')
                        model_template.add_reactions([rea_temp])
                        model_gap.add_reactions([rea_temp])

                        rea_str = gap_biomass_met + ' --> '
                        model_template.reactions.get_by_id('object').reaction = rea_str
                        model_gap.reactions.get_by_id('object').reaction = rea_str

                        # print('biomass apart optimize: ',gap_biomass_met)
                        model_template.objective = "object"
                        # print('Lreuteri_530:', model_template.optimize())
                        model_gap.objective = "object"
                        # print('model_i_draft_3_refined:',model_gap.optimize())

                        if model_gap.optimize().objective_value < 1e-10:
                            solution_part_biomass = cobra.flux_analysis.gapfill(model_gap, model_template,
                                                                                demand_reactions=False)
                            gaps_set = gaps_set | set([i.id for i in solution_part_biomass[0]])
                except:
                    pass

    if model_i_draft_3.optimize().objective_value == 0.0:
        print('----- Gap fill failed !!! -----')
        for i in Lreuteri_530.optimize().fluxes.index:
            if model_i_draft_3.reactions.get_by_id(i).reaction != Lreuteri_530.reactions.get_by_id(i).reaction:
                print(model_i_draft_3.reactions.get_by_id(i).reaction, Lreuteri_530.reactions.get_by_id(i).reaction)
            model_i_draft_3.reactions.get_by_id(i).bounds = Lreuteri_530.reactions.get_by_id(i).bounds

    model_i_draft_3.objective = "BIOMASS"
    print('model_i_draft_3 Biomass:', model_i_draft_3.optimize())
    model_i_draft_3.id = sp_name
    cobra.io.save_json_model(model_i_draft_3, 'GEM_from_templates_' + sp_name + '_draft.json', sort='True')

    print('=====  Done =====')
