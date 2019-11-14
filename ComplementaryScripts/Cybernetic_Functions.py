#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 11/7/19

"""Cybernetic_Functions.py
:description : script to build cybernetic model
:param : 
:returns: 
:rtype: 
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import odeint
import __main__



# %% <def functions>

class Cybernetic_Model(dict):

    def __init__(self, name):
        # dict.__init__(self)

        self['name'] = self.name = name
        # self = dict.fromkeys(['Smz', 'x0', 'kmx', 'K', 'ke', 'alpha', 'beta', 'n_carbon'], 0)
        self['Smz'] = None
        self['x0'] = None
        self['kmax'] = None
        self['K'] = None
        self['ke'] = None
        self['alpha'] = None
        self['beta'] = None
        self['n_carbon'] = None
        self['sub_index'] = None
        self['Biomass_index'] = None

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


def dxdy(x, t, CB_model):
    # Smz=Smz, kmax=kmax, K=K, ke=ke, alpha=alpha, beta=beta
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

    rM, rE, rG = __main__.rate_def(x, CB_model)

    # rM, rE, rG = __main__.rate_def2(x, CB_model).rM , rate_def2(x, CB_model).rE , rate_def2(x, CB_model).rG

    # cybernetic_var = abs(Smz[sub_index, :] * rM * n_carbon)
    cybernetic_var = rM * n_carbon
    # cybernetic_var[cybernetic_var==0] = 1
    cybernetic_var[cybernetic_var < 0] = 0
    if sum(cybernetic_var) > 0:
        u = cybernetic_var / sum(cybernetic_var)
        v = cybernetic_var / np.max(abs(cybernetic_var))
    else:
        u = v = np.zeros(n_path)

    u[u < 0.0001] = 1.0
    v[v < 0.0001] = 1.0


    V = np.eye(n_path) * v

    Growth_rate = rG * v
    mu = sum(Growth_rate)

    dy_dx_mets = Smz @ V @ rM * x[Biomass_index]
    # dy_dx_enzyme = alpha + rE * u - (beta + mu) * x[n_mets:];

    mumax = Smz[Biomass_index, :].T * kmax
    dy_dx_enzyme = (mumax + beta) / (alpha + ke) * (alpha + rE * u) - (beta + mu) * x[n_mets:];

    dxdt_vector = list(dy_dx_mets) + list(dy_dx_enzyme)
    if abs(t - 5.0) < 0.1:
        print('t', t)
        print('rM', rM)
        print('u', u)
        print('v', v)
        print('x', x)


    return dxdt_vector


def cb_model_simulate(CB_model, tspan, draw=True):
    initial_x0 = CB_model.x0
    sol = odeint(dxdy, initial_x0, tspan, args=(CB_model,))

    if draw:
        for key in range(0, CB_model.Smz.shape[0]):
            model_line = plt.plot(tspan, sol[:, key], color="k", linewidth=2)
        plt.legend((model_line[0],), ("HCM FBA",), fontsize=18)
        plt.xlabel("Time (hr)", fontsize=20)
        plt.ylabel("Abundance (mM)", fontsize=20)
        plt.show()
    return sol


if __name__ == '__main__':
    # %% < def input:>
    os.chdir('../ComplementaryData/')
    case = 2
    if case == 1:
        # DefineTime
        tStart = 0.0
        tStop = 15
        tStep = 0.1
        tspan = np.linspace(tStart, tStop, (tStop - tStart) / tStep)

        # matrix Z: reactions x pathways; and Smz metabolites x pathways , S metabolites x reactions : Smz = Sm @ Z
        print('\n---------- Loading FBA modes ... ---------- ')
        Z = np.genfromtxt('Case1_ecoli_reduced/FBA_em_reduced.csv', delimiter=',')
        Z = Z.T
        Smz = Z
        Smz[6, :] = Smz[6, :] - Smz[10, :]
        Smz = Smz[[0, 5, 6, 7, 8, 9, 11], :]
        Smz[0, :] = -Smz[0, :]

        # metabolites and pathways number
        (n_mets, n_path) = Smz.shape

        # experiment data
        print('\n---------- Loading Experiment Data ... ---------- ')
        experiment_data_df = pd.read_csv('Case1_ecoli_reduced/ecoli_reduce_experiment_data.txt', delimiter='\t',
                                         header=0)
        metabObj = ['glc', 'succ', 'for', 'lac', 'ac', 'etoh', 'biomass', ]

        # initial metabolites at tome 0, t0: initial_mets.shape = (n_mets,)
        initial_mets = experiment_data_df[metabObj].values[0, :]

        # initial:Enzyme: initial_enzyme.shape = (n_path,)
        initial_enzyme = np.array([0.9] * n_path)

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
            0.346572,
            0.0124164,
            0.032764,
            0.049244,
            0.203273,
            0.249942,
            7.90653,
        ])

        # K : n_path
        K = np.array([
            1.19695,  # 1
            10.3871,  # 2
            0.312798,  # 3
            5.42755,  # 4
            4.40557,  # 5
            11.9541,  # 6
            10.5902,  # 7 Formate
        ])

        # carbon number for each pathways
        n_carbon = np.array([6, 6, 6, 6, 6, 6, 0])

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
        ecoli_reduced_cb['sub_index'] = 0
        ecoli_reduced_cb['Biomass_index'] = 6
        CB_model = ecoli_reduced_cb


        def rate_def(x, CB_model):
            Smz = CB_model.Smz
            kmax = CB_model.kmax
            K = CB_model.K
            ke = CB_model.ke
            alpha = CB_model.alpha
            beta = CB_model.beta
            Biomass_index = CB_model.Biomass_index
            sub_index = CB_model['sub_index']

            (n_mets, n_path) = Smz.shape

            rM = [
                kmax[0] * x[0 + n_mets] * x[sub_index] / (K[0] + x[sub_index]),
                kmax[1] * x[1 + n_mets] * x[sub_index] / (K[1] + x[sub_index]),
                kmax[2] * x[2 + n_mets] * x[sub_index] / (K[2] + x[sub_index]),
                kmax[3] * x[3 + n_mets] * x[sub_index] / (K[3] + x[sub_index]),
                kmax[4] * x[4 + n_mets] * x[sub_index] / (K[4] + x[sub_index]),
                kmax[5] * x[5 + n_mets] * x[sub_index] / (K[5] + x[sub_index]),
                kmax[6] * (x[2] ** 2) / ((K[6] ** 2) + (x[2] ** 2)),
            ]

            rE = ke * rM / kmax / np.array(list(x[n_mets:-1]) + [1])

            rG = Smz[Biomass_index, :] * rM[:]  # 11: biomass index

            return rM, rE, rG
        #

    # case 2 DefineTime

    if case == 2:
        tStart = 0.0
        tStop = 15
        tStep = 0.1
        tspan = np.linspace(tStart, tStop, (tStop - tStart) / tStep)

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
        metabObj = ['glc', 'succ', 'for', 'lac', 'ac', 'etoh', 'biomass', ]

        # initial metabolites at tome 0, t0: initial_mets.shape = (n_mets,)
        initial_mets = experiment_data_df[metabObj].values[0, :]

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
            11.03,
            19.18,
            4.1210e-05,
            2.18,
            8.4682e-05,
            19.8
        ])

        K = np.array([
            0.002,
            6.38,
            5.27,
            0.002,
            0.002,
            15,
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
        ecoli_reduced_cb['sub_index'] = 0
        ecoli_reduced_cb['Biomass_index'] = 1
        CB_model = ecoli_reduced_cb


        def rate_def(x, CB_model):
            Smz = CB_model.Smz
            kmax = CB_model.kmax
            K = CB_model.K
            ke = CB_model.ke
            alpha = CB_model.alpha
            beta = CB_model.beta
            Biomass_index = CB_model.Biomass_index
            sub_index = CB_model['sub_index']

            (n_mets, n_path) = Smz.shape

            rM = [
                kmax[0] * x[0 + n_mets] * x[sub_index] / (K[0] + x[sub_index]),
                kmax[1] * x[1 + n_mets] * x[sub_index] / (K[1] + x[sub_index]),
                kmax[2] * x[2 + n_mets] * x[sub_index] / (K[2] + x[sub_index]),
                kmax[3] * x[3 + n_mets] * x[sub_index] / (K[3] + x[sub_index]),
                kmax[4] * x[4 + n_mets] * x[sub_index] / (K[4] + x[sub_index]),
                kmax[5] * (x[3] ** 2) / ((K[5] ** 2) + (x[3] ** 2)),
            ]

            rE = ke * rM / kmax / np.array(list(x[n_mets:-1]) + [1])

            rG = Smz[Biomass_index, :] * rM[:]  # 11: biomass index

            return rM, rE, rG
        #
    sol = cb_model_simulate(CB_model, tspan, draw=True)
