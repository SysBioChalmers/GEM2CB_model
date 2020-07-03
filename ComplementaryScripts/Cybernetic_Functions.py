#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 11/7/19

"""Cybernetic_Functions.py
:description : script to build cybernetic model
:param : 
:returns: 
:rtype: 
"""
import copy
import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.integrate import odeint


# %% <def functions>

class Cybernetic_Model(dict):

    def __init__(self, name):
        # dict.__init__(self)

        self['name'] = self.name = name
        # self = dict.fromkeys(['Smz', 'x0', 'kmx', 'K', 'ke', 'alpha', 'beta', 'n_carbon'], 0)
        self['Smz'] = None
        self['mets_name'] = None
        self['initial_mets'] = None
        self['initial_enzymes'] = None
        self['x0'] = None
        self['kmax'] = None
        self['K'] = None
        self['ke'] = None
        self['alpha'] = None
        self['beta'] = None
        self['n_carbon'] = None
        self['sub_index'] = None
        self['Biomass_index'] = None
        self['experiment_data_df'] = pd.DataFrame()
        self['weight_of_method_t_f'] = 0
        self['weights_of_mets'] = 1
        self['weights_of_time'] = 1
        self['xtol'] = 0.8

    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError(r"'Model' object has no attribute '%s'" % key)

    def __setattr__(self, key, value):
        if key in self.keys():
            self[key] = value
        else:
            raise AttributeError(r"'Model' object has no attribute '%s'" % key)

    def check_type(self):
        for key in ['kmax', 'ke', 'alpha', 'beta', 'n_carbon']:
            if type(self[key]) != np.ndarray:
                try:
                    self[key] = np.array(self[key])
                except:
                    print("Warning: type(self[%s]) != np.ndarray, check it " % key)

    def check_shape(self):
        (n_mets, n_path) = self['Smz'].shape
        for key in ['kmax', 'K', 'ke', 'alpha', 'beta', 'n_carbon']:
            if self[key].shape != (n_path,):
                print('self[%s].shape != (%d, check it )' % key, n_path)
        if self['initial_x0'].shape != n_mets + n_mets:
            print('initial_x0 shape !=(%d, check it )' % (n_mets + n_mets,))

    def weight_of_method(self, exp_y, model_y, res_sq):

        if self.weight_of_method_t_f == 0:  # minimize variance
            # weight_of_method = 1  # no weight_of_method
            res_sq_weight_of_method = res_sq

        elif self.weight_of_method_t_f == 1:  # from AUMIC
            weight_of_method = 1 / np.mean(res_sq, axis=0) ** 2 / res_sq.shape[0]
            weight_of_method = weight_of_method / sum(weight_of_method)
            res_sq_weight_of_method = res_sq * weight_of_method

        elif self.weight_of_method_t_f == 2:  # minimize relative variance
            res_sq_weight_of_method = ((exp_y - model_y) / (np.median(exp_y, axis=0))) ** 2

        elif self.weight_of_method_t_f == 2.1:  # minimize relative variance
            res_sq_weight_of_method = ((exp_y - model_y) / (np.mean(exp_y, axis=0))) ** 2


        elif self.weight_of_method_t_f == 2.5:  # minimize relative variance
            # weight_of_method = 1 / (exp_y ** 2 + 1e-20)
            exp_y_ = np.where(exp_y <= 0, 0.2, exp_y)  # avoid division by zero

            res_sq_weight_of_method = ((res_sq) / exp_y_) ** 2


        elif self.weight_of_method_t_f == 3:  # model_y and exp_y , Normalization by col note : error 0/0

            model_y_norm = (model_y - model_y.min(axis=0) + 1e-20) / (
                    model_y.max(axis=0) - model_y.min(axis=0) + 1e-20)
            exp_y_norm = (exp_y - exp_y.min(axis=0) + 1e-20) / (
                    exp_y.max(axis=0) - exp_y.min(axis=0) + 1e-20)

            weight_of_method = (exp_y_norm - model_y_norm) ** 2
            res_sq_weight_of_method = res_sq * weight_of_method

        elif self.weight_of_method_t_f == 4:  # model_y and exp_y , Normalization by col note : error 0/0

            model_y_norm = (model_y - model_y.min(axis=0) + 1e-20) / (
                    model_y.max(axis=0) - model_y.min(axis=0) + 1e-20)
            exp_y_norm = (exp_y - exp_y.min(axis=0) + 1e-20) / (
                    exp_y.max(axis=0) - exp_y.min(axis=0) + 1e-20)

            weight_of_method = (exp_y_norm - model_y_norm) ** 2

            weight_of_method = 1 / np.mean(weight_of_method, axis=0) ** 2 / weight_of_method.shape[0]
            weight_of_method = weight_of_method / sum(weight_of_method)
            res_sq_weight_of_method = res_sq * weight_of_method


        elif self.weight_of_method_t_f == 5:  # tol for each pathways
            exp_y_ = np.where(exp_y <= 0, 1, exp_y)  # avoid division by zero
            relative_variance = ((res_sq) / exp_y_) ** 2
            res_sq_weight_of_method = np.where(relative_variance <= (self.xtol) ** 2, 0, 1)
            res_sq_weight_of_method = np.array([np.max(res_sq_weight_of_method, axis=0)])
            res_sq_weight_of_method = np.concatenate((res_sq_weight_of_method,
                                                      np.array([np.mean(res_sq_weight_of_method, axis=0)]))
                                                     , axis=1)


        elif self.weight_of_method_t_f == 6:  # tol for each experiment point
            exp_y_ = np.where(exp_y <= 0, 1, exp_y)  # avoid division by zero
            relative_variance = ((res_sq) / exp_y_) ** 2
            res_sq_weight_of_method = np.where(relative_variance <= (self.xtol) ** 2, 0, 1)

        elif self.weight_of_method_t_f == -1:  # res_sq Normalization all
            weight_of_method = (res_sq - res_sq.min() + 1e-20) / (res_sq.max() - res_sq.min() + 1e-20) / (
                    res_sq + 1e-20)  # relative error squire

        elif self.weight_of_method_t_f == -2:  # res_sq Normalization by col note : error 0/0
            weight_of_method = (res_sq - res_sq.min(axis=0) + 1e-20) / (
                    res_sq.max(axis=0) - res_sq.min(axis=0) + 1e-20) / (res_sq + 1e-20)  # relative error squire

        elif self.weight_of_method_t_f == -1:  # res_sq Normalization all
            weight_of_method = (res_sq - res_sq.min() + 1e-20) / (res_sq.max() - res_sq.min() + 1e-20) / (
                    res_sq + 1e-20)  # relative error squire

        #
        # elif self.weight_of_method_t_f == -4:        #   Normalization by col
        #     weight_of_method = (res_sq - res_sq.min(axis=0) )/(res_sq.max(axis=0) -  res_sq.min(axis=0) ) / res_sq# relative error squire
        #
        # elif self.weight_of_method_t_f == -5:        #   Normalization all
        #     weight_of_method = (res_sq - res_sq.min() )/(res_sq.max() -  res_sq.min() ) / res_sq# relative error squire

        else:
            weight_of_method = 1

            res_sq_weight_of_method = res_sq * weight_of_method

        return res_sq_weight_of_method

    def infmations(self):
        print('''
        Matrix_Smz: 
        InitialConditions_x0
        K
        kmax
        EnzymeRate_ke
        EnzymeSynthesis_alpha
        EnzymeDegradation_beta
        n_carbon
        ''')

    def expend_model(self, model2, name='connected'):
        '''connect two model
        '''
        model1 = self
        model = copy.deepcopy(model1)
        model.name = name
        if model1.Smz.shape[0] != model2.Smz.shape[0]:
            print('Metabolites munber not equal!!! please check')
        model.Smz = np.concatenate((model1.Smz, model2.Smz), axis=1)
        model.mets_name = model1.mets_name
        model.initial_enzymes = np.concatenate((model1.initial_enzymes, model2.initial_enzymes), )
        model.x0 = np.concatenate((model.initial_mets, model.initial_enzymes))
        model.kmax = np.concatenate((model1.kmax, model2.kmax), )
        model.K = np.concatenate((model1.K[0:-1], model2.K), )
        model.ke = np.concatenate((model1.ke, model2.ke), )
        model.alpha = np.concatenate((model1.alpha, model2.alpha), )
        model.beta = np.concatenate((model1.beta, model2.beta), )
        model.n_carbon = np.concatenate((model1.n_carbon, model2.n_carbon), )
        model.sub_index = np.concatenate((model1.sub_index, model2.sub_index), )
        if len(model1.Biomass_index) != 1 and len(model2.Biomass_index) != 1:
            model.Biomass_index = np.concatenate((model1.Biomass_index, model2.Biomass_index), )
        return model


def dxdy(x, t, CB_model):
    # Smz=Smz, kmax=kmax, K=K, ke=ke, alpha=alpha, beta=beta
    # TODO the result i different from matlab.... :(
    Smz = CB_model.Smz
    kmax = CB_model.kmax
    K = CB_model.K
    ke = CB_model.ke
    alpha = CB_model.alpha
    beta = CB_model.beta
    Biomass_index = CB_model.Biomass_index
    sub_index = CB_model.sub_index
    n_carbon = CB_model.n_carbon

    (n_mets, n_path) = Smz.shape

    rM, rE, rG = CB_model.rate_def(x)

    # rM, rE, rG = __main__.rate_def2(x, CB_model).rM , rate_def2(x, CB_model).rE , rate_def2(x, CB_model).rG

    try:
        u, v = CB_model.cybernetic_var_def(rM)

    except:
        # print('stand cybernetic_vars')
        # cybernetic_var = abs(Smz[sub_index, :] * rM * n_carbon)
        cybernetic_var = rM * n_carbon
        # cybernetic_var[cybernetic_var==0] = 1
        cybernetic_var = np.array(cybernetic_var)
        cybernetic_var[cybernetic_var < 0] = 0
        if sum(cybernetic_var) > 0:
            u = cybernetic_var / sum(cybernetic_var)
            v = cybernetic_var / np.max(abs(cybernetic_var))
        else:
            u = v = np.zeros(n_path)

    Growth_rate = rG * v
    mu = sum(Growth_rate)

    # V = np.eye(n_path) * v
    # dy_dx_mets = Smz @ V @ (rM * x[Biomass_index]) ##Smz @ V @ (rM * x[Biomass_index]) == Smz @ (v * rM * x[Biomass_index])
    dy_dx_mets = Smz @ (v * rM * x[Biomass_index])

    # dy_dx_enzyme = alpha + rE * u - (beta + mu) * x[n_mets:];
    # Smz[Biomass_index, :] == np.array([Smz1[Biomass_index[i], i] for i in np.arange(0,n_path)])
    if type(Biomass_index) != int and len(Biomass_index) == n_path:
        mumax = kmax * np.array([Smz[Biomass_index[i], i] for i in np.arange(0, n_path)])
    else:
        mumax = kmax * Smz[Biomass_index, :]
    dy_dx_enzyme = (mumax + beta) / (alpha + ke) * (alpha + rE * u) - (beta + mu) * x[n_mets:];

    dxdt_vector = list(dy_dx_mets) + list(dy_dx_enzyme)
    # if abs(t - 0.0) < 0.01:
    #     print('t', t)
    #     print('rM', rM)
    #     print('(mumax + beta) / (alpha + ke) * (alpha + rE * u)', (mumax + beta) / (alpha + ke) * (alpha + rE * u))
    #     print('(beta + mu) * x[n_mets:]', (beta + mu) * x[n_mets:])
    #     print('x', x)

    return dxdt_vector


def cb_model_simulate(CB_model, tspan, draw=True):
    try:
        CB_model.x0 = np.concatenate((CB_model.initial_mets, CB_model.initial_enzymes))
    except:
        pass
    initial_x0 = CB_model.x0
    sol = odeint(dxdy, initial_x0, tspan, args=(CB_model,))

    if draw:
        fig, ax = plt.subplots()
        ax.plot(tspan, sol[:, np.arange(CB_model.Smz.shape[0])])

        if CB_model.experiment_data_df.shape != (0, 0):
            experiment_data_df = CB_model.experiment_data_df
            plt.gca().set_prop_cycle(None)
            ax.plot(experiment_data_df.iloc[:, 0], experiment_data_df[CB_model.mets_name], 'o')

        try:
            ax.legend(CB_model.mets_name)
        except:
            pass
            ax.legend(("Simulation",))
        ax.set_xlabel("Time (h)", fontsize=20)
        ax.set_ylabel("Abundance (mM)", fontsize=20)
        ax.set_title(CB_model.name)
        fig.show()
    return sol


def update_paras_func(x_paras, _CB_models, para_to_fit, retern_model_or_paras='model', full_output=False):
    '''
    from para_to_fit get x_paras from model, return x_paras
    or apply x_paras to model, return model
    '''

    num_kmax = len(para_to_fit['kmax'])
    num_K = len(para_to_fit['K'])
    if full_output:
        print('kmax = np.array(', list(x_paras)[0:num_kmax], ')')
        print('K = np.array(', list(x_paras)[num_kmax:], ')')
    if type(_CB_models) != list:
        _CB_models = [_CB_models]

    if retern_model_or_paras == 'model':
        _CB_models = [copy.deepcopy(i) for i in _CB_models]

        if num_kmax > 0:
            for i in range(0, len(_CB_models)):
                index = para_to_fit['kmax'][:, i]
                valuae = x_paras[0:num_kmax][~np.isnan(index)]
                index = index[~np.isnan(index)]
                index = index.astype(int)
                _CB_models[i]['kmax'][index] = valuae

        if num_K > 0:
            for i in range(0, len(_CB_models)):
                index = para_to_fit['K'][:, i]
                valuae = x_paras[num_kmax:num_kmax + num_K][~np.isnan(index)]
                index = index[~np.isnan(index)]
                index = index.astype(int)
                _CB_models[i]['K'][index] = valuae

        return _CB_models

    else:
        x_paras = []  # [0] * (num_kmax + num_K)
        if num_kmax > 0:
            for kmax_i in para_to_fit['kmax']:
                for i in range(0, len(_CB_models)):
                    if not np.isnan(kmax_i[i]):
                        x_paras.append(_CB_models[i]['kmax'][int(kmax_i[i])])
                        break

        if num_K > 0:
            for K_i in para_to_fit['K']:
                for i in range(0, len(_CB_models)):
                    if not np.isnan(K_i[i]):
                        x_paras.append(_CB_models[i]['K'][int(K_i[i])])
                        break

        # paras_initial = x_paras
        return x_paras
    # TODO return other para_to_fit keys


def residuals_func(x_paras, _CB_models, para_to_fit, exp_ys, tspan, time_points_indexs, metas_indexs,
                   draw=False, full_output=False):
    ''':arg get residuals

    '''

    # _CB_model = update_paras_func(x_paras, _CB_model, para_to_fit, retern_model_or_paras='model')
    _CB_models = update_paras_func(x_paras, _CB_models, para_to_fit, retern_model_or_paras='model',
                                   full_output=full_output)

    sol_temps = []
    fs = np.array([])
    for i in range(0, len(_CB_models)):
        _CB_model = _CB_models[i]
        time_points_index = time_points_indexs[i]
        metas_index = metas_indexs[i]

        sol_temp = cb_model_simulate(_CB_model, tspan, draw=False)
        sol_temps.append(sol_temp)

        model_y = sol_temp[time_points_index, :][:, metas_index]

        exp_y = exp_ys[i]

        res = exp_y - model_y
        res_sq = res ** 2

        res_sq_weight_of_method = _CB_model.weight_of_method(exp_y, model_y, res_sq)
        # res_sq_weight_of_method = res_sq * weight_of_method

        weights_of_mets = _CB_model.weights_of_mets

        f = sum(res_sq_weight_of_method * weights_of_mets)
        f = np.array(f)
        if i == 0:
            fs = f
        else:
            fs = np.concatenate((fs, f))

        # try:
        #
        # except:
        #     fs =

        if full_output:
            print('res:', res)
            print('res_sq:', res_sq)
            print('weight_of_method:', weight_of_method)
            print('weights_of_mets:', weights_of_mets)
            print('f:', f)
            print('fs:', fs)
        # print('fs', fs)
        # f = sum(abs(exp_y - model_y)) * weights
        # print('f', f)
        # print('x_paras',x_paras)

    if draw:
        global fitting_fig, fitting_axs, model_lines
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        for i_model in range(0, len(_CB_models)):
            _CB_model = _CB_models[i_model]
            fitting_ax = fitting_axs[i_model]
            sol_temp = sol_temps[i_model]
            # plt.gca().set_prop_cycle(None)
            for i in range(0, _CB_model.Smz.shape[0]):
                try:
                    model_lines[i_model][i][0].set_alpha(0.3)
                except:
                    pass
                # model_lines[i][0].remove()
                model_line_i = fitting_ax.plot(tspan, sol_temp[:, i], color=colors[i])
                model_lines[i_model][i] = model_line_i
        plt.pause(1e-30)
        # fitting_fig.show()
    # return 0
    # return fs
    result = abs(sum(fs))

    return result


def parameters_fitting(CB_models, experiment_data_dfs, para_to_fit, tspan, draw=False, full_output=False,
                       method='Nelder-Mead', jac=None, hess=None, hessp=None, bounds=None, constraints=(), tol=None,
                       callback=None, options={'xtol': 0.1, 'ftol': 0.1, 'maxiter': 1000, 'disp': True}):
    options_def = {'xtol': 0.1, 'ftol': 0.1, 'maxiter': 1000, 'disp': True}
    for i in options_def.keys():
        if i not in options.keys():
            options[i] = options_def[i]

    if type(CB_models) != list:
        print('convert to model list')
        CB_models = [CB_models]
    if type(experiment_data_dfs) != list:
        experiment_data_dfs = [experiment_data_dfs]
    _CB_models = []
    exp_ys = []
    sol_initials = []
    time_points_indexs = []
    metas_indexs = []
    for index_i in range(0, len(CB_models)):
        _CB_model = copy.deepcopy(CB_models[index_i])
        _CB_models.append(_CB_model)
        experiment_data_df = experiment_data_dfs[index_i]
        time_points = experiment_data_df.iloc[:, 0]
        time_points_index = []

        for time_point in time_points:
            index = np.argmin(abs(tspan - time_point))
            time_points_index.append(index)
        time_points_indexs.append(time_points_index)
        metas_index = []
        for exp_met_name in experiment_data_df.columns:
            if exp_met_name in _CB_model.mets_name:
                index = _CB_model.mets_name.index(exp_met_name)
                metas_index.append(index)
        metas_indexs.append(metas_index)
        exp_y = experiment_data_df[_CB_model.mets_name].values
        exp_ys.append(exp_y)

        sol_initial = cb_model_simulate(_CB_model, tspan, draw=False)
        sol_initials.append(sol_initial)
    paras_initial = update_paras_func([], _CB_models, para_to_fit, retern_model_or_paras='paras')
    # paras_initials.append(paras_initial)

    # model_y = sol[timepoints_index,:][:,metas_index]
    print('paras Fitting')
    if draw:
        global fitting_fig, fitting_axs, model_lines

        print('drawing')
        # plt.set_cmap('jet')
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        plt.gca().set_prop_cycle(None)
        # matplotlib.use("qt5agg")
        matplotlib.use("macosx")
        # import matplotlib.pyplot as plt

        # matplotlib.use("agg")
        plt.ion()
        fitting_fig = plt.figure()
        fitting_axs = []
        model_lines = []
        for i_model in range(0, len(CB_models)):
            fitting_ax = fitting_fig.add_subplot(len(CB_models), 1, i_model + 1)
            fitting_ax.set_xlabel("Time (h)")
            fitting_ax.set_ylabel("Concentration (mM)")
            exp_ponts = []
            model_i_lines = []
            for i in range(0, len(_CB_model.mets_name)):
                met_i_name = _CB_model.mets_name[i]
                exp_pont_i = fitting_ax.plot(experiment_data_dfs[i_model].iloc[:, 0], exp_ys[i_model][:, i], 'o',
                                             color=colors[i], label=met_i_name)
                exp_ponts.append(exp_pont_i)
                # plt.gca().set_prop_cycle(None)
                model_line_i = fitting_ax.plot(tspan, sol_initials[i_model][:, i], color=colors[i], alpha=0.5)
                model_line_i = fitting_ax.plot(tspan, sol_initials[i_model][:, i], color=colors[i], alpha=0.5)
                model_i_lines.append(model_line_i)
            model_lines.append(model_i_lines)

            fitting_ax.legend()
            fitting_ax.set_title(CB_models[i_model].name)
            fitting_axs.append(fitting_ax)
            # fitting_fig.show()
        plt.pause(0.01)
        # plt.show()
        # plt.show()
        # fitting_fig.clf()
        # plt.close()

    # minimum = scipy.optimize.fmin(residuals_func, paras_initial, args=(
    #     _CB_models, para_to_fit, exp_ys, tspan, time_points_indexs, metas_indexs, draw, full_output),
    #                               xtol=0.01, ftol=0.01, maxiter=maxiter, full_output=True)
    # print(minimum[0])
    # print(minimum[1])

    if bounds == None and method in ['SLSQP', 'L-BFGS-B', 'TNC']:
        bounds = [(0.1, None)] * len(paras_initial)

    minimum = scipy.optimize.minimize(residuals_func, paras_initial, args=(
        _CB_models, para_to_fit, exp_ys, tspan, time_points_indexs, metas_indexs, draw, full_output),
                                      method=method, jac=jac, hess=hess, hessp=hessp, bounds=bounds,
                                      constraints=constraints, tol=tol,
                                      callback=callback, options=options)
    print(list(minimum.x))
    print(minimum.fun)

    # if draw:
    #     # import matplotlib.pyplot as plt
    #     plt.ioff()
    #     plt.show()

    return minimum


if __name__ == '__main__':
    # %% < def input:>
    os.chdir('../ComplementaryData/')

    tStart = 0.0
    tStop = 10
    tStep = 0.1
    tspan = np.linspace(tStart, tStop, int(((tStop - tStart) / tStep) + 1))

    # matrix Z: reactions x pathways; and Smz metabolites x pathways , S metabolites x reactions : Smz = Sm @ Z
    print('\n---------- Loading FBA modes ... ---------- ')
    # Z = np.genfromtxt('Case1_ecoli_reduced/FBA_em_reduced.csv', delimiter=',')
    # Z = Z.T
    # Smz = Z
    # Smz[6, :] = Smz[6, :] - Smz[10, :]
    # Smz = Smz[[0, 5, 6, 7, 8, 9, 11], :]
    # Smz[0, :] = -Smz[0, :]
    Smz = np.array([[-1., -1., -1., -1., -1.,
                     0.],
                    [0.01123911, 0.02796746, 0.02796746, 0.01960329, 0.02294896,
                     0.],
                    [0., 0.8673643, 0.8673643, 0.01980717, 0.89116456,
                     0.],
                    [1., 1.62104086, 0., 0.01980717, 1.15071905,
                     -1.],
                    [1., 0.75367656, 0.75367656, 0., 0.25955448,
                     0.],
                    [0., 0., 0., 1.70458824, 0.,
                     0.],
                    [0.57722474, 0., 0., 0., 0.53832256,
                     0.]])

    # metabolites and pathways number
    (n_mets, n_path) = Smz.shape

    # experiment data
    print('\n---------- Loading Experiment Data ... ---------- ')
    experiment_data_df = pd.read_csv('Case1_ecoli_reduced/ecoli_reduce_experiment_data.txt', delimiter='\t',
                                     header=0)
    metabObj = ['glc', 'biomass', 'ac', 'for', 'etoh', 'lac', 'succ', ]

    # initial metabolites at tome 0, t0: initial_mets.shape = (n_mets,)
    initial_mets = experiment_data_df[metabObj].values[0, :]
    # initial_mets = initial_mets[[0, 6, 4, 2, 5, 3, 1]]

    # initial:Enzyme: initial_enzyme.shape = (n_path,)
    initial_enzyme = np.array([0.9, 0.9, 0.9, 0.9, 0.9, 1])

    # initial data x0 :initial_x0.shape  = (n_mets + n_path,)
    initial_x0 = np.concatenate((initial_mets, initial_enzyme))

    # Enzyme Rate Parameters: alpha,beta,ke : de/dt =  alpha + rE(ke) * u - (beta + mu) * e
    alpha = np.array([0.04] * n_path)
    beta = np.array([0.05] * n_path)
    ke = np.array([0.620342] * n_path)  # or 0.5

    # Metabolites rate Parameters kmax , Ki : dm/dt =  Smz @ V @ rM(kmax,K) * c
    # TODO the equations not the same , should defined by user

    # kmax : n_path

    kmax = np.array([
        10,
        35.17,
        4.1210e-09,
        2.31,
        4.12e-7,
        21.08
    ])

    # K : n_path
    K = np.array([
        2.5235,
        24.39,
        33.50,
        1,
        1,
        1.042,
    ])

    # carbon number for each pathways
    n_carbon = np.array([6, 6, 6, 6, 6, 0])

    # Construct the model:
    ecoli_reduced_cb = Cybernetic_Model('CB model for Ecoli reduced matrix ')
    ecoli_reduced_cb.Smz = Smz
    ecoli_reduced_cb.x0 = initial_x0
    ecoli_reduced_cb.kmax = kmax
    ecoli_reduced_cb.K = K
    ecoli_reduced_cb.ke = ke
    ecoli_reduced_cb.alpha = alpha
    ecoli_reduced_cb.beta = beta
    ecoli_reduced_cb.n_carbon = n_carbon
    ecoli_reduced_cb.sub_index = 0
    ecoli_reduced_cb.Biomass_index = 1
    CB_model = copy.deepcopy(ecoli_reduced_cb)
    CB_model.mets_name = metabObj


    # CB_model = copy.deepcopy(ecoli_reduced_cb)
    # CB_model.kmax[[1,2,3]] = [11,22,33]
    # CB_model.K[[1,2]] = [77,55]

    # if num_kmax > 0:
    #     index = para_to_fit['kmax']
    #     CB_model.kmax[index] = [11,22,33]
    # num_K = len(para_to_fit['K'])
    # if num_K > 0:
    #     index = para_to_fit['K']
    #     CB_model.K[index] = [44,55]

    def rate_def(self, x):
        # print('def rate')
        Smz = self.Smz
        kmax = self.kmax
        K = self.K
        ke = self.ke
        alpha = self.alpha
        beta = self.beta
        Biomass_index = self.Biomass_index
        sub_index = self.sub_index

        (n_mets, n_path) = Smz.shape

        r_kin_basic = [
            kmax[0] * x[sub_index] / (K[0] + x[sub_index]),
            kmax[1] * x[sub_index] / (K[1] + x[sub_index]),
            kmax[2] * x[sub_index] / (K[2] + x[sub_index]),
            kmax[3] * x[sub_index] / (K[3] + x[sub_index]),
            kmax[4] * x[sub_index] / (K[4] + x[sub_index]),
            kmax[5] * (x[3] ** 2) / ((K[5] ** 2) + (x[3] ** 2)), ]

        rM = r_kin_basic * x[n_mets:]

        rE = ke * r_kin_basic / kmax

        rG = Smz[Biomass_index, :] * rM[:]  # 11: biomass index

        return rM, rE, rG


    setattr(Cybernetic_Model, 'rate_def', rate_def)


    def cybernetic_var_def(self, rM):
        # print('def cybernetic_var')
        (n_mets, n_path) = self.Smz.shape
        # cybernetic_var = abs(Smz[sub_index, :] * rM * n_carbon)
        cybernetic_var = rM * self.n_carbon
        # cybernetic_var[cybernetic_var==0] = 1
        cybernetic_var[cybernetic_var < 0] = 0
        if sum(cybernetic_var) > 0:
            u = cybernetic_var / sum(cybernetic_var)
            v = cybernetic_var / np.max(abs(cybernetic_var))
        else:
            u = v = np.zeros(n_path)

        if CB_model.name in ['CB model for Ecoli reduced matrix ', 'CB model for Ecoli iML1515 matrix ']:
            u[-1] = 1.0
            v[-1] = 1.0
        return u, v


    setattr(Cybernetic_Model, 'cybernetic_var_def', cybernetic_var_def)
    # CB_model.cybernetic_var_def = cybernetic_var_def

    CB_model.experiment_data_df = experiment_data_df
    sol = cb_model_simulate(CB_model, tspan, draw=True)

    # %%
    para_to_fit = {'kmax': np.array([[0], [1], [4], [5]]), 'K': np.array([[1], [2], [4], [5]])}

    # Branch_work update_paras_func
    # paras_initials = update_paras_func([], [CB_model], para_to_fit, retern_model_or_paras='paras')
    # print(paras_initials)
    # print('kmax',kmax,'K',K)
    # [model_2] = update_paras_func(np.arange(0,8), [CB_model], para_to_fit, retern_model_or_paras='model')
    #
    # print(model_2.kmax)
    # print(model_2.K)
    CB_model.weight_of_method_t_f = False
    minimum = parameters_fitting([CB_model], [experiment_data_df], para_to_fit, tspan, draw=True,
                                 full_output=False, options={'xtol': 0.01, 'ftol': 0.01, 'maxiter': 10, 'disp': True})
