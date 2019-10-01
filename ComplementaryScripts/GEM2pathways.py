#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 9/24/19

"""GEM2pathways.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import cobra
import My_def
import matplotlib.pyplot as plt
import numpy as np
import statistics
from cobra.flux_analysis import production_envelope
# from cobra.flux_analysis import flux_variability_analysis
import pandas as pd
from cobra.core import Configuration
from pandas import DataFrame
from optlang.symbolics import Zero
from numpy import zeros
import multiprocessing
from warnings import warn
from cobra.flux_analysis.parsimonious import add_pfba
from cobra.util import solver as sutil
from cobra.flux_analysis.loopless import loopless_fva_iter
import pandas as pd



CONFIGURATION = Configuration()

def _init_worker(model, loopless, sense):
    """Initialize a global model object for multiprocessing."""
    global _model
    global _loopless
    _model = model
    _model.solver.objective.direction = sense
    _loopless = loopless


def _fva_step(reaction_id):
    global _model
    global _loopless
    rxn = _model.reactions.get_by_id(reaction_id)
    # The previous objective assignment already triggers a reset
    # so directly update coefs here to not trigger redundant resets
    # in the history manager which can take longer than the actual
    # FVA for small models
    _model.solver.objective.set_linear_coefficients(
        {rxn.forward_variable: 1, rxn.reverse_variable: -1})
    _model.slim_optimize()
    sutil.check_solver_status(_model.solver.status)
    if _loopless:
        value = loopless_fva_iter(_model, rxn)
        #flux =
    else:
        value = _model.solver.objective.value
        pfba_solution = cobra.flux_analysis.pfba(_model)
        fluxes = pfba_solution.fluxes

    _model.solver.objective.set_linear_coefficients(
        {rxn.forward_variable: 0, rxn.reverse_variable: 0})
    return reaction_id, value , fluxes


def flux_variability_analysis2(model, reaction_list=None, loopless=False,
                              fraction_of_optimum=1.0, pfba_factor=None,
                              processes=None):
    """
    most code of this function From cobrapy !!!!!!!

    Determine the minimum and maximum possible flux value for each reaction.

    Parameters
    ----------
    model : cobra.Model
        The model for which to run the analysis. It will *not* be modified.
    reaction_list : list of cobra.Reaction or str, optional
        The reactions for which to obtain min/max fluxes. If None will use
        all reactions in the model (default).
    loopless : boolean, optional
        Whether to return only loopless solutions. This is significantly
        slower. Please also refer to the notes.
    fraction_of_optimum : float, optional
        Must be <= 1.0. Requires that the objective value is at least the
        fraction times maximum objective value. A value of 0.85 for instance
        means that the objective has to be at least at 85% percent of its
        maximum.
    pfba_factor : float, optional
        Add an additional constraint to the model that requires the total sum
        of absolute fluxes must not be larger than this value times the
        smallest possible sum of absolute fluxes, i.e., by setting the value
        to 1.1 the total sum of absolute fluxes must not be more than
        10% larger than the pFBA solution. Since the pFBA solution is the
        one that optimally minimizes the total flux sum, the ``pfba_factor``
        should, if set, be larger than one. Setting this value may lead to
        more realistic predictions of the effective flux bounds.
    processes : int, optional
        The number of parallel processes to run. If not explicitly passed,
        will be set from the global configuration singleton.

    Returns
    -------
    pandas.DataFrame
        A data frame with reaction identifiers as the index and two columns:
        - maximum: indicating the highest possible flux
        - minimum: indicating the lowest possible flux

    Notes
    -----
    This implements the fast version as described in [1]_. Please note that
    the flux distribution containing all minimal/maximal fluxes does not have
    to be a feasible solution for the model. Fluxes are minimized/maximized
    individually and a single minimal flux might require all others to be
    suboptimal.

    Using the loopless option will lead to a significant increase in
    computation time (about a factor of 100 for large models). However, the
    algorithm used here (see [2]_) is still more than 1000x faster than the
    "naive" version using ``add_loopless(model)``. Also note that if you have
    included constraints that force a loop (for instance by setting all fluxes
    in a loop to be non-zero) this loop will be included in the solution.

    References
    ----------
    .. [1] Computationally efficient flux variability analysis.
       Gudmundsson S, Thiele I.
       BMC Bioinformatics. 2010 Sep 29;11:489.
       doi: 10.1186/1471-2105-11-489, PMID: 20920235

    .. [2] CycleFreeFlux: efficient removal of thermodynamically infeasible
       loops from flux distributions.
       Desouki AA, Jarre F, Gelius-Dietrich G, Lercher MJ.
       Bioinformatics. 2015 Jul 1;31(13):2159-65.
       doi: 10.1093/bioinformatics/btv096.
    """
    if reaction_list is None:
        reaction_ids = [r.id for r in model.reactions]
    else:
        reaction_ids = [r.id
                        for r in model.reactions.get_by_any(reaction_list)]

    if processes is None:
        processes = CONFIGURATION.processes
    num_reactions = len(reaction_ids)
    processes = min(processes, num_reactions)

    fva_result = DataFrame({
        "minimum": zeros(num_reactions, dtype=float),
        "maximum": zeros(num_reactions, dtype=float)
    }, index=reaction_ids)


    fva_fluxes_result = DataFrame([])


    prob = model.problem
    with model:
        # Safety check before setting up FVA.
        model.slim_optimize(error_value=None,
                            message="There is no optimal solution for the "
                                    "chosen objective!")
        # Add the previous objective as a variable to the model then set it to
        # zero. This also uses the fraction to create the lower/upper bound for
        # the old objective.
        # TODO: Use utility function here (fix_objective_as_constraint)?
        if model.solver.objective.direction == "max":
            fva_old_objective = prob.Variable(
                "fva_old_objective",
                lb=fraction_of_optimum * model.solver.objective.value)
        else:
            fva_old_objective = prob.Variable(
                "fva_old_objective",
                ub=fraction_of_optimum * model.solver.objective.value)
        fva_old_obj_constraint = prob.Constraint(
            model.solver.objective.expression - fva_old_objective, lb=0, ub=0,
            name="fva_old_objective_constraint")
        model.add_cons_vars([fva_old_objective, fva_old_obj_constraint])

        if pfba_factor is not None:
            if pfba_factor < 1.:
                warn("The 'pfba_factor' should be larger or equal to 1.",
                     UserWarning)
            with model:
                add_pfba(model, fraction_of_optimum=0)
                ub = model.slim_optimize(error_value=None)
                flux_sum = prob.Variable("flux_sum", ub=pfba_factor * ub)
                flux_sum_constraint = prob.Constraint(
                    model.solver.objective.expression - flux_sum, lb=0, ub=0,
                    name="flux_sum_constraint")
            model.add_cons_vars([flux_sum, flux_sum_constraint])

        model.objective = Zero  # This will trigger the reset as well
        for what in ("minimum", "maximum"):
            if processes > 1:
                # We create and destroy a new pool here in order to set the
                # objective direction for all reactions. This creates a
                # slight overhead but seems the most clean.
                chunk_size = len(reaction_ids) // processes
                pool = multiprocessing.Pool(
                    processes,
                    initializer=_init_worker,
                    initargs=(model, loopless, what[:3])
                )
                for rxn_id, value ,fluxes in pool.imap_unordered(_fva_step,
                                                         reaction_ids,
                                                         chunksize=chunk_size):
                    fva_result.at[rxn_id, what] = value
                    fva_fluxes_result[rxn_id + what] = fluxes
                pool.close()
                pool.join()
            else:
                _init_worker(model, loopless, what[:3])
                for rxn_id, value, fluxes in map(_fva_step, reaction_ids):
                    fva_result.at[rxn_id, what] = value
                    fva_fluxes_result[rxn_id + what] = fluxes

    return fva_result[["minimum", "maximum"]],fva_fluxes_result

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


#fva_value_result,fva_fluxes_result = fva_self_def(model, reaction_list=yield_rea_ids, pfba=False,)

def get_fba_flux_df(model, biomass_rea_id,carbon_rea_ids,yield_rea_ids,step_of_biomass,):

    #%% <Step1 check data>
    # initial fba check
    model.objective = biomass_rea_id
    solution = model.optimize()
    print('model.optimize():\t', solution.objective_value)

    #%% < Step2 get yield dataFrame(matrix/ points for MYA/Next function)>   FVA

    #pfba_values_all = pd.DataFrame()
    pfba_fluxes_all = pd.DataFrame()
    for persent_i in np.arange(0,1+1/step_of_biomass,1/step_of_biomass):

        #set biomass bounds

        model.reactions.get_by_id(biomass_rea_id).bounds = (solution.objective_value*persent_i,solution.objective_value*persent_i)
        #pfba_values,pfba_fluxes = flux_variability_analysis2(model,yield_rea_ids)     #FVA
        pfba_values,pfba_fluxes = fva_self_def(model, reaction_list=yield_rea_ids, pfba=False,)

        # if np.isnan(np.sum(pfba_fluxes)) and persent_i == 0:        # somtimen biomass = 0  no solutions
        #     persent_i = 1e-6
        #     model.reactions.get_by_id(biomass_rea_id).bounds = (solution.objective_value*persent_i,solution.objective_value*persent_i)
        #     pfba_values,pfba_fluxes = flux_variability_analysis2(model,yield_rea_ids)     #FVA

        #pfba_values['biomass'] = solution.objective_value*persent_i
        #pfba_values_all = pfba_values_all.append(pfba_values)

        pfba_fluxes.columns = [column + '_' + str(round(persent_i,2)) for column in pfba_fluxes.columns ]
        pfba_fluxes_all = pd.concat([pfba_fluxes_all, pfba_fluxes], axis=1)

    yield_df = pfba_fluxes_all.loc[ carbon_rea_ids+[biomass_rea_id]+ yield_rea_ids, : ]     #select reactions
    yield_normalized_df = yield_df.copy()

    yield_df_values  = yield_normalized_df.values      #normalize and sort
    yield_df_values = yield_df_values/abs(yield_df_values[0,:])

    yield_normalized_df.loc[:,:] = yield_df_values
    yield_normalized_df = yield_normalized_df.sort_values(by = [biomass_rea_id]+yield_rea_ids,axis = 1,)
    return yield_normalized_df



# %%
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

# yield_value,fluxes = get_yield_m(model2,production_reaid,carbon_rea_id,max_or_min = 'max')

def get_fba_yield_df(model2, biomass_rea_id,carbon_rea_id,yield_rea_id,step_of_biomass,):


    # <Step1 check data and find the biomass yield max and min>
    # initial fba check
    model = model2.copy()
    model.objective = biomass_rea_id
    solution = model.optimize()
    print('model.optimize():\t', solution.objective_value)

    biomass_yield_value_max,_ = get_yield_m(model,biomass_rea_id,carbon_rea_id,max_or_min = 'max')
    biomass_yield_value_min,_ = get_yield_m(model,biomass_rea_id,carbon_rea_id,max_or_min = 'min')


    # < Step2 get yield dataFrame(matrix/ points for MYA/Next function)>

    fluxes_all = pd.DataFrame()

    for persent_i in np.arange(0,1+1/step_of_biomass,1/step_of_biomass):

        yield_k = (biomass_yield_value_max + biomass_yield_value_min) * persent_i

        #set biomass yield
        model = model2.copy()
        yield_point_flux = model.problem.Constraint(
            model.reactions.get_by_id(biomass_rea_id).flux_expression + yield_k * model.reactions.get_by_id(carbon_rea_id).flux_expression,
            lb=0,
            ub=0)
        model.add_cons_vars(yield_point_flux)

        for rea_id in yield_rea_ids:
            _,fluxes_max = get_yield_m(model,rea_id,carbon_rea_id,max_or_min = 'max')
            _,fluxes_min = get_yield_m(model,rea_id,carbon_rea_id,max_or_min = 'min')
            # fluxes_all = pd.concat([fluxes_all, fluxes_max], axis=1)
            # fluxes_all = pd.concat([fluxes_all, fluxes_min], axis=1)
            fluxes_all[rea_id + 'max'+str(persent_i)] = fluxes_max
            fluxes_all[rea_id + 'min'+str(persent_i)] = fluxes_min

    yield_df = fluxes_all.loc[ [carbon_rea_id]+[biomass_rea_id]+ yield_rea_ids, : ]     #select reactions
    yield_normalized_df = yield_df.copy()
    # yield_normalized_df = yield_normalized_df.T.drop_duplicates().T

    yield_df_values  = yield_normalized_df.values      #normalize and sort
    yield_df_values = yield_df_values/abs(yield_df_values[0,:])

    yield_normalized_df.loc[:,:] = yield_df_values
    yield_normalized_df = yield_normalized_df.sort_values(by = [biomass_rea_id]+yield_rea_ids,axis = 1,)
    return yield_normalized_df





# %%
if __name__ == '__main__':

    #%% <Step1 set data> model, biomass_rea_id,carbon_rea_ids,yield_rea_ids,step_of_biomass

    e_coli_core = cobra.io.read_sbml_model('../ComplementaryData/e_coli_core.xml')
    #e_coli_core.reactions.get_by_id('EX_o2_e').bounds = (0.0,1000.0)

    biomass_rea_id = 'BIOMASS_Ecoli_core_w_GAM'
    carbon_rea_ids = ['EX_glc__D_e']
    yield_rea_ids = ['EX_ac_e','EX_for_e','EX_etoh_e','EX_lac__D_e','EX_succ_e',]

    step_of_biomass = 10

    for carbon_rea_id in carbon_rea_ids:        #set carbon lb as -10
        e_coli_core.reactions.get_by_id(carbon_rea_id).bounds = (-10.0,0.01)
    model = e_coli_core


    #%% < Step2 get yield dataFrame(matrix/ points for MYA/Next function)>   FVA
    # yield_normalized_df = get_fba_flux_df(model, biomass_rea_id,carbon_rea_ids,yield_rea_ids,step_of_biomass)
    yield_normalized_df = get_fba_yield_df(model, biomass_rea_id,carbon_rea_ids[0],yield_rea_ids,step_of_biomass)



    #%% <Step3 plot initial yield>: TODO not divided by carbon, so not yield

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
        ax.set_ylabel(exmet_reaction+'/glc')
    ax.set_xlabel('biomass/glc')
    fig.show()
