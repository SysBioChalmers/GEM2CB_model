#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-09-12

"""Step1_iML1515_test.py
:description : script
:param : 
:returns: 
:rtype: 
"""
import cobra
import matplotlib.pyplot as plt
import GEM2pathways
import ConvexHull_yield


e_coli_core = cobra.io.read_sbml_model('../ComplementaryData/e_coli_core.xml')
e_coli_core.reactions.get_by_id('EX_o2_e').bounds = (0.0,1000.0)

biomass_rea_id = 'BIOMASS_Ecoli_core_w_GAM'
carbon_rea_ids = ['EX_glc__D_e']
yield_rea_ids = ['EX_ac_e','EX_for_e','EX_etoh_e','EX_lac__D_e','EX_succ_e',]

step_of_biomass = 10

for carbon_rea_id in carbon_rea_ids:        #set carbon lb as -10
    e_coli_core.reactions.get_by_id(carbon_rea_id).bounds = (-10.0,0.0)
model = e_coli_core


#%% < Step2 get yield dataFrame(matrix/ points for MYA/Next function)>   FVA
yield_normalized_df = GEM2pathways.get_fba_yield_df(model, biomass_rea_id,carbon_rea_ids,yield_rea_ids,step_of_biomass)



#%% <plot initial yield>: TODO not divided by carbon, so not yield

def trim_axs(axs, N):
    """little helper to massage the axs list to have correct length..."""
    axs = axs.flat
    for ax in axs[N:]:
        ax.remove()
    return axs[:N]

cols = 3
rows = len(yield_rea_ids) // cols + 1
figsize = (10, 8)
fig, axs = plt.subplots(cols, rows,figsize=figsize)
axs = trim_axs(axs, len(yield_rea_ids))

for ax , exmet_reaction in zip(axs,yield_rea_ids):
    temp_columns_max = [column for column in yield_normalized_df.columns if (exmet_reaction in column) and ('max' in column)]
    temp_columns_min = [column for column in yield_normalized_df.columns if (exmet_reaction in column) and ('min' in column)]
    x = yield_normalized_df.loc[biomass_rea_id,temp_columns_max].values
    y_max = yield_normalized_df.loc[exmet_reaction,temp_columns_max]
    y_min = yield_normalized_df.loc[exmet_reaction,temp_columns_min]

    ax.plot(x, y_max,'x-',color = 'black',alpha=0.5)
    ax.plot(x, y_min,'x-',color = 'black',alpha=0.5)
    # temp_df.plot.area(x='biomass', y=['maximum','minimum'],label=['max', 'max'],color=['r', 'w'],color=['tab:blue', 'b'],stacked=False);
    ax.set_ylabel(exmet_reaction)
ax.set_xlabel('biomass')
fig.show()

# %% <SETP3 yield analysis >:

all_points = yield_normalized_df.values.T
all_points = all_points[:,1:]       #exclude the carbon colume

experiment_datas = [ ]      # TODO experiment data!!!
qhull_options = 'Qt QJ Pp Qw Qx'      #'QG0'mean expect point 0(index)
cutoff_persent = 0.99

indexes, hulls , weightss , estimated_datas, in_hulls = ConvexHull_yield.pipeline_mya(all_points, experiment_datas , qhull_options = 'Qt QJ Pp Qw Qx', cutoff_persent = 0.99)






