#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 9/25/19

"""step_0_prepare_ecoli_reduced_GEM.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import cobra

#make a GEM
'''
model:

-METINT
G6P T3P PEP PYR AcCoA ATP ADP NADH NAD CoA

-METEXT
GLU SUC FOR ACT LAC ETH B CO2 H2

R1 : GLU + PEP => G6P + PYR .
R2 : G6P + ATP => 2 T3P + ADP .
R3 : G6P + 6 NAD => T3P + 6 NADH .
R4 : T3P + NAD + ADP => PEP + NADH + ATP .
R5 : PEP + ADP => PYR + ATP .
R6 : PEP + CO2 + 2 NADH => SUC + 2 NAD .
R7 : PYR + CoA => AcCoA + FOR .
R8 : PYR + NADH => LAC + NAD .
R9 : AcCoA + ADP => ACT + CoA + ATP .
R10 : AcCoA + 2 NADH => ETH + CoA + 2 NAD .
R11 : FOR => CO2 + H2 .
R12 : 6.775 G6P + 82.2 ATP + 4.065 NADH => B + 82.2 ADP + 4.065 NAD .

reaids = ['R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12']
yield_rea_ids = ['EX_ac_e','EX_for_e','EX_etoh_e','EX_lac__D_e','EX_succ_e',]
yield_rea_ids = ['R1','R12','R9','R7','R10','R8','R6']
'''


reaids = ['R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12']
reaids_all = reaids +['EX_GLU','EX_SUC','EX_FOR','EX_ACT','EX_LAC','EX_ETH','EX_B','EX_CO2','EX_H2']

ecoli_reduced_model = cobra.Model('ecoli_reduced_model')
for reaid in reaids_all :
    locals()[reaid] = cobra.Reaction(reaid)
    ecoli_reduced_model.add_reaction(locals()[reaid])

euations = {
    'R1' : 'GLU + PEP --> G6P + PYR',
    'R2' : 'G6P + ATP --> 2 T3P + ADP',
    'R3' : 'G6P + 6 NAD --> T3P + 6 NADH',
    'R4' : 'T3P + NAD + ADP --> PEP + NADH + ATP',
    'R5' : 'PEP + ADP --> PYR + ATP',
    'R6' : 'PEP + CO2 + 2 NADH --> SUC + 2 NAD',
    'R7' : 'PYR + CoA --> AcCoA + FOR',
    'R8' : 'PYR + NADH --> LAC + NAD',
    'R9' : 'AcCoA + ADP --> ACT + CoA + ATP',
    'R10' : 'AcCoA + 2 NADH --> ETH + CoA + 2 NAD',
    'R11' : 'FOR --> CO2 + H2',
    'R12' : '6.775 G6P + 82.2 ATP + 4.065 NADH --> B + 82.2 ADP + 4.065 NAD',
    'EX_GLU' : 'GLU <-> ',
    'EX_SUC' : 'SUC <--> ',
    'EX_FOR' : 'FOR <--> ',
    'EX_ACT' : 'ACT <--> ',
    'EX_LAC' : 'LAC <--> ',
    'EX_ETH' : 'ETH <--> ',
    'EX_B' : 'B <--> ',
    'EX_CO2' : 'CO2 <--> ',
    'EX_H2' : 'H2 <--> '}

for reaid, euqation in euations.items():
    ecoli_reduced_model.reactions.get_by_id(reaid).reaction = euqation


ecoli_reduced_model.reactions.get_by_id('R1').bounds = (0.0,10.0)
ecoli_reduced_model.objective = 'R12'
ecoli_reduced_model.optimize()
for met in ecoli_reduced_model.metabolites:
    met.compartment = 'c'

cobra.io.write_sbml_model(ecoli_reduced_model,'../ComplementaryData/case1_ecoli_reduced/ecoli_reduced_model.xml')
cobra.io.save_json_model(ecoli_reduced_model,'../ComplementaryData/case1_ecoli_reduced/ecoli_reduced_model.json')

# ecoli_reduced_model2 = cobra.io.load_json_model('../ComplementaryData/case1_ecoli_reduced/ecoli_reduced_model_R.json')
#
# ecoli_reduced_model3 = cobra.io.read_sbml_model('../ComplementaryData/case1_ecoli_reduced/ecoli_reduced_model_R.xml')

# %%

