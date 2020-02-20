#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 10/29/19

"""nothing.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import multiprocessing

def add(x, y = 1, z = 2):
	v = x+y+z
	list = [x,y,z]
	return v,list

# Get all worker processes
cores = multiprocessing.cpu_count()

# Start all worker processes
# pool = multiprocessing.Pool(processes=cores)
# x1 = list(range(3))
# y1 = list(range(3))
# z1 = list(range(3))
#
# tasks = [(x,y) for x in x1 for y in y1 ]
#
# a = []
# b = []
# for i in pool.starmap(add,tasks):
# 	a.append(i[0])
# 	b.append(i[1])
# print(a,b)
# pool.close()
# pool.join()

#
# pool = multiprocessing.Pool(processes=cores)
# tasks = [(y,x,z) for x in x1 for y in y1 for z in z1]
# a = []
# b = []
# for i in pool.starmap(add,tasks):
# 	a.append(i[0])
# 	b.append(i[1])
# print(a,b)

import cobra.test
model = cobra.test.create_test_model("textbook")
# flux_variability_analysis(model, model.reactions[:10])

import cvxpy as cp

# %%

import numpy as np
from cvxopt import solvers, matrix, spmatrix
from scipy import sparse


def scipy_sparse_to_spmatrix(A):
    coo = A.tocoo()
    SP = spmatrix(coo.data, coo.row.tolist(), coo.col.tolist())
    return SP


def spmatrix_sparse_to_scipy(A):
    data = np.array(A.V).squeeze()
    rows = np.array(A.I).squeeze()
    cols = np.array(A.J).squeeze()
    return sparse.coo_matrix((data, (rows, cols)))


def sparse_None_vstack(A1, A2):
    if A1 is None:
        return A2
    else:
        return sparse.vstack([A1, A2])


def numpy_None_vstack(A1, A2):
    if A1 is None:
        return A2
    else:
        return np.vstack([A1, A2])


def numpy_None_concatenate(A1, A2):
    if A1 is None:
        return A2
    else:
        return np.concatenate([A1, A2])


def get_shape(A):
    if isinstance(C, spmatrix):
        return C.size
    else:
        return C.shape


def numpy_to_cvxopt_matrix(A):
    if A is None:
        return A
    if sparse.issparse(A):
        if isinstance(A, sparse.spmatrix):
            return scipy_sparse_to_spmatrix(A)
        else:
            return A
    else:
        if isinstance(A, np.ndarray):
            if A.ndim == 1:
                return matrix(A, (A.shape[0], 1), 'd')
            else:
                return matrix(A, A.shape, 'd')
        else:
            return A


def cvxopt_to_numpy_matrix(A):
    if A is None:
        return A
    if isinstance(A, spmatrix):
        return spmatrix_sparse_to_scipy(A)
    elif isinstance(A, matrix):
        return np.array(A).squeeze()
    else:
        return np.array(A).squeeze()


def lsqlin(C, d, reg=0, A=None, b=None, Aeq=None, beq=None, \
           lb=None, ub=None, x0=None, opts=None):
    '''
        Solve linear constrained l2-regularized least squares. Can
        handle both dense and sparse matrices. Matlab's lsqlin
        equivalent. It is actually wrapper around CVXOPT QP solver.

            min_x ||C*x  - d||^2_2 + reg * ||x||^2_2
            s.t.  A * x <= b
                  Aeq * x = beq
                  lb <= x <= ub

        Input arguments:
            C   is m x n dense or sparse matrix
            d   is n x 1 dense matrix
            reg is regularization parameter
            A   is p x n dense or sparse matrix
            b   is p x 1 dense matrix
            Aeq is q x n dense or sparse matrix
            beq is q x 1 dense matrix
            lb  is n x 1 matrix or scalar
            ub  is n x 1 matrix or scalar

        Output arguments:
            Return dictionary, the output of CVXOPT QP.

        Dont pass matlab-like empty lists to avoid setting parameters,
        just use None:
            lsqlin(C, d, 0.05, None, None, Aeq, beq) #Correct
            lsqlin(C, d, 0.05, [], [], Aeq, beq) #Wrong!
    '''
    sparse_case = False
    if sparse.issparse(A):  # detects both np and cxopt sparse
        sparse_case = True
        # We need A to be scipy sparse, as I couldn't find how
        # CVXOPT spmatrix can be vstacked
        if isinstance(A, spmatrix):
            A = spmatrix_sparse_to_scipy(A)

    C = numpy_to_cvxopt_matrix(C)
    d = numpy_to_cvxopt_matrix(d)
    Q = C.T * C
    q = - d.T * C
    nvars = C.size[1]

    if reg > 0:
        if sparse_case:
            I = scipy_sparse_to_spmatrix(sparse.eye(nvars, nvars, \
                                                    format='coo'))
        else:
            I = matrix(np.eye(nvars), (nvars, nvars), 'd')
        Q = Q + reg * I

    lb = cvxopt_to_numpy_matrix(lb)
    ub = cvxopt_to_numpy_matrix(ub)
    b = cvxopt_to_numpy_matrix(b)

    if lb is not None:  # Modify 'A' and 'b' to add lb inequalities
        if lb.size == 1:
            lb = np.repeat(lb, nvars)

        if sparse_case:
            lb_A = -sparse.eye(nvars, nvars, format='coo')
            A = sparse_None_vstack(A, lb_A)
        else:
            lb_A = -np.eye(nvars)
            A = numpy_None_vstack(A, lb_A)
        b = numpy_None_concatenate(b, -lb)
    if ub is not None:  # Modify 'A' and 'b' to add ub inequalities
        if ub.size == 1:
            ub = np.repeat(ub, nvars)
        if sparse_case:
            ub_A = sparse.eye(nvars, nvars, format='coo')
            A = sparse_None_vstack(A, ub_A)
        else:
            ub_A = np.eye(nvars)
            A = numpy_None_vstack(A, ub_A)
        b = numpy_None_concatenate(b, ub)

    # Convert data to CVXOPT format
    A = numpy_to_cvxopt_matrix(A)
    Aeq = numpy_to_cvxopt_matrix(Aeq)
    b = numpy_to_cvxopt_matrix(b)
    beq = numpy_to_cvxopt_matrix(beq)

    # Set up options
    if opts is not None:
        for k, v in opts.items():
            solvers.options[k] = v

    # Run CVXOPT.SQP solver
    sol = solvers.qp(Q, q.T, A, b, Aeq, beq, None, x0)

    return sol


# C = points_2d.T
# %%
C = np.array([[1.47568331, 0.126809, 1.72088113, 1.86881024, 1.96879662,
               1.71787763, 1.57111798, 1.02675484, 0.35520492, 0.79717899,
               0.2678625, 0.0617791, 1.87828341, 0.60261213, 0.59106767,
               0.66587256, 0.93413637, 1.29639681, 0.05045636, 1.68441322,
               1.11806509, 1.7081999, 0.69575839, 0.8920533, 0.10847897,
               0.35421507, 1.32561612, 0.66165799, 1.79697228, 0.2363104],
              [1.97683586, 1.0799642, 1.41383484, 1.99898324, 0.57569869,
               0.82904508, 0.92967988, 1.52791416, 1.63640808, 0.20044308,
               0.35623391, 0.71926983, 0.11340938, 1.04377135, 0.67169795,
               0.35133806, 0.41789335, 1.81030712, 1.35078235, 0.9369364,
               1.82426495, 0.20802315, 1.49109215, 1.47253491, 1.12372285,
               0.3683882, 1.1944227, 0.59987398, 0.26824587, 0.42520307]])
d = [0.5, 0.5]
# set the pramater of lsqlin
# ret = lsqlin.lsqlin(C_in, d_in, 0, None, None, Aeq, beq, \
#         lb=None, ub=None, x0=None, opts=None
(mC, nC) = C.shape
lb = np.zeros(nC)
ub = np.ones(nC)
A = None
b = None
x0 = None

# normalize C and d
C0 = C  # initial C and d
d0 = np.array(d)

Aeq = np.ones((1, nC))
beq = np.ones((1,))
ret = lsqlin(C0, d0, 0, None, None, Aeq, beq, \
             lb, ub, None, {'show_progress': False})
weights = ret['x'].T
x1 = weights

x = cp.Variable(nC)
objective = cp.Minimize(cp.sum_squares(C0 * x - d0))
constraints = [x <= ub,
               x >= lb,
               Aeq @ x == beq,
               ]
prob = cp.Problem(objective, constraints)
result = prob.solve(solver=cp.CPLEX)
x2 = x.value

a = np.where(x2 > 1e-6)[0]

if len(a) > mC + 1:
    print('inhull')
    Aeq = np.vstack([Aeq, C0])
    beq = np.vstack([beq, d0.reshape(2, 1)])
    beq = beq.reshape(3)
    C = np.eye(nC)  # np.array([mC]*nC)
    d = np.ones(nC)
    x = cp.Variable(nC)
    objective = cp.Minimize(cp.norm1(x))
    constraints = [x <= ub,
                   x >= lb,
                   Aeq @ x <= beq * 1.0,
                   Aeq @ x >= beq * 1.0
                   ]
    prob = cp.Problem(objective, constraints)
    result = prob.solve(solver=cp.CPLEX)
    x2 = x.value

else:
    print('not in_hull')
#
# temp_dt = np.tile(d0, (nC, 1))
# C = C0 * (1 / temp_dt.T)
# d = np.ones((1, mC))[0]
#
#
# ret = lsqlin(C, d, 0, None, None, Aeq, beq, \
#              lb, ub, None, {'show_progress': False})
# weights = ret['x'].T
#
# x = cp.Variable(nC)
# objective = cp.Minimize(cp.sum_squares(C*x - d ))
# constraints = [x <= ub,
#                x >= lb,
#                Aeq @ x == beq,
#                ]
# prob = cp.Problem(objective, constraints)
# result = prob.solve(solver=cp.CPLEX)
# x2 = x.value
#
#
#
#
# ret = lsqlin(C, d, 0, None, None, Aeq, beq, \
#              lb, ub, None, {'show_progress': False})
# weights = ret['x'].T
#
#
#
#
# %%in op1
#
# Aeq =  np.vstack([Aeq, C0])
# beq = np.vstack([beq, d0.reshape(2,1)])
# beq = beq.reshape(3)
# C = np.array([mC]*nC)
# d = np.ones(1)
# x = cp.Variable(nC)
# objective = cp.Minimize(cp.sum_squares(C*x - d ))
# constraints = [x <= ub,
#                x >= lb,
#             Aeq@x <= beq+0.01,
#              Aeq@x >= beq-0.01,
#                ]
# prob = cp.Problem(objective, constraints)
# result = prob.solve(solver=cp.CPLEX)
# x2 = x.value
#
# %% opt2
# Aeq =  np.vstack([Aeq, C0])
# beq = np.vstack([beq, d0.reshape(2,1)])
# beq = beq.reshape(3)
# C = np.eye(nC) #np.array([mC]*nC)
# d = np.ones(nC)
# x = cp.Variable(nC)
# objective = cp.Minimize(cp.sum_squares(-C*x + d ))
# constraints = [x <= ub,
#                x >= lb,
#                 Aeq@x <= beq*1.1,
#                Aeq@x >= beq*0.9
#                ]
# prob = cp.Problem(objective, constraints)
# result = prob.solve(solver=cp.CPLEX)
# x2 = x.value


# %% opt3


if False:
    temp_dt = np.tile(d0, (nC, 1))
    C = C0 * (1 / temp_dt.T)
    d = np.ones((1, mC))[0]

    in_hull = True
    if in_hull:
        #
        # lb = np.ones(nC)
        # ub = np.ones(nC)*np.inf

        C_in = np.eye(nC)  # one C and zero d

        d_in = np.ones((1, nC))[0]
        d_in = np.ones((1, nC))[0] * 0.5
        # nored
        Aeq = np.ones((1, nC))
        # beq = np.ones((1, 1))
        Aeq = np.vstack([Aeq, C])
        beq = np.ones((Aeq.shape[0]))

        ret = lsqlin(C_in, d_in, 0, None, None, Aeq, beq, \
                     lb, ub, None, {'show_progress': False})
        # print (ret['x'].T)
        weights = ret['x'].T

        x = cp.Variable(nC)
        objective = cp.Minimize(cp.sum_squares(C_in * x - d_in))
        # c = np.array([0.5]*nC)
        # objective = cp.Minimize(c@x)

        constraints = [x <= ub,
                       x >= lb,
                       Aeq @ x == beq,
                       ]
        prob = cp.Problem(objective, constraints)
        result = prob.solve()
        x2 = x.value

    if not in_hull:
        # Bug: why there are too many pathways !!!!
        # lb = np.ones(nC)
        # ub = np.zeros(nC)
        Aeq = np.ones((1, nC))
        beq = np.ones((1, 1))
        ret = lsqlin(C, d, 0, None, None, Aeq, beq, \
                     lb, ub, None, {'show_progress': False})
        # ret = lsqlin.lsqnonneg(C, d, {'show_progress': False})

        # print (ret['x'].T)
        weights = ret['x'].T

        x = cp.Variable(nC)
        objective = cp.Minimize(cp.sum_squares(C * x - d))
        constraints = [x <= ub,
                       x >= lb,
                       Aeq @ x == beq,
                       ]
        prob = cp.Problem(objective, constraints)
        result = prob.solve()
        x2 = x.value

print(np.around(x1, 4))
# print(C @ x1.T)
print(d0)
print(C0 @ x1.T)
print(sum(x1))

print('----------------------------')

print(np.around(x2, 4))
# print(C @ x2.T)
print(d0)
print(C0 @ x2.T)
print(sum(x2))
print(sum((C0 @ x.value - d0) ** 2))

# %%erro

# %% ode example:

y0, t0 = [1.0j, 2.0], 0


def f(t, y, arg1):
    return [1j * arg1 * y[0] + y[1], -arg1 * y[1] ** 2]


def jac(t, y, arg1):
    return [[1j * arg1, 1], [0, -arg1 * 2 * y[1]]]


r = ode(f, jac).set_integrator('zvode', method='bdf')
r.set_initial_value(y0, t0).set_f_params(2.0).set_jac_params(2.0)
t1 = 10
dt = 1
while r.successful() and r.t < t1:
    print(r.t + dt, r.integrate(r.t + dt))


# %% < examples 2 >


# Define a function which calculates the derivative
def dy_dx(y, x):
    return x - y


xs = np.linspace(0, 5, 100)
y0 = 1.0  # the initial condition
ys = odeint(dy_dx, y0, xs)
# ys = np.array(ys).flatten()
# plt.rcParams.update({'font.size': 14})  # increase the font size
plt.xlabel("x")
plt.ylabel("y")
plt.plot(xs, ys);
plt.show()

# %%
import numpy as np
from scipy.integrate import odeint
import scipy

Ab = np.array([[-0.25, 0, 0],
               [0.25, -0.2, 0],
               [0, 0.2, -0.1]])
Ab = [1, 2, 3]


def deriv(A, t, Ab):
    c = [0, 0, 0]
    c[0] = Ab[0] + A[0]
    c[1] = Ab[1] + A[1]
    c[2] = Ab[2] + A[2]
    return c
    # return np.dot(Ab, A)


time = np.linspace(0, 25, 101)
A0 = [10, 20, 30]

MA = odeint(deriv, A0, time, args=(Ab,))

MA2 = scipy.integrate.ode(deriv).set_integrator('vode', method='bdf', order=15)
MA2.set_f_params(*args)

# %%
from scipy.integrate import ode
import numpy as np

y0, t0 = [1.0j, 2.0], 0


def f(t, y, arg1):
    return [1j * arg1 * y[0] + y[1], -arg1 * y[1] ** 2]


def jac(t, y, arg1):
    return [[1j * arg1, 1], [0, -arg1 * 2 * y[1]]]


r = ode(f, jac).set_integrator('zvode', method='bdf')
r.set_initial_value(y0, t0).set_f_params(2.0).set_jac_params(2.0)
t1 = 10
dt = 1
while r.successful() and r.t < t1:
    print(r.t + dt, r.integrate(r.t + dt))


def pend(t, y, b, c):
    theta, omega = y
    dydt = [omega, -b * omega - c * np.sin(theta)]
    return dydt


b = 0.25
c = 5.0
y0 = [np.pi - 0.1, 0.0]
t = np.linspace(0, 10, 101)
from scipy.integrate import odeint

sol = odeint(pend, y0, t, args=(b, c), tfirst=True)

import matplotlib.pyplot as plt

plt.plot(t, sol[:, 0], 'b', label='theta(t)')
plt.plot(t, sol[:, 1], 'g', label='omega(t)')
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()

# %%
from scipy.integrate import odeint
from scipy.integrate import ode
import numpy as np


def dy_dt(t, y, a, b, c):
    dydt = [y[1], 1000 * (1 - y[0] ** 2) * y[1] - y[0]]
    return dydt


# y(2); 1000*(1-y(1)^2)*y(2)-y(1)
a = 1
b = 0
c = 1000.0
y0 = [2, 0]
t = np.linspace(0, 3000, 30000)

sol = odeint(dy_dt, y0, t, args=(a, b, c), tfirst=True)

import matplotlib.pyplot as plt

plt.plot(t, sol[:, 0], 'b', label='1')
plt.plot(t, sol[:, 1], 'g', label='2')
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()

import julia
from julia.api import Julia

jl = Julia(compiled_modules=False)
from diffeqpy import de


def f(u, p, t):
    return -u


u0 = 0.5
tspan = (0., 1.)
prob = de.ODEProblem(f, u0, tspan)
sol = de.solve(prob)

r = ode(dy_dt, jac=None).set_integrator('vode', method='bdf', order=15)
r.set_f_params(a, b, c)
r.set_initial_value(y0, 0)

t1 = 1000
dt = 1
while r.successful() and r.t < t1:
    plt.plot(r.t + dt, r.integrate(r.t + dt)[0], 'b', label='1')
    plt.plot(r.t + dt, r.integrate(r.t + dt)[1], 'g', label='2')
    # print(r.t+dt, r.integrate(r.t+dt))
plt.show()

# %%
from scipy.interpolate import UnivariateSpline


def f(y, t, a1, a2, a3, b1, b2, b3, c1, c2, c3):
    dy1_dt = a1(t) * y[1] - b1(t) * y[2] - c1(t) * y[0]
    dy2_dt = a2(t) * y[0] - b2(t) * y[0] - c2(t) * y[2]
    dy3_dt = a3(t) * y[1] - b3(t) * y[2] - c3(t) * y[0]
    return list([dy1_dt, dy2_dt, dy3_dt])


y0 = [1, 0.5, 0.02]
t = np.linspace(0, 10.0, 1000)
a1 = b1 = c1 = UnivariateSpline(t, -t / 10)
a2 = b2 = c2 = UnivariateSpline(t, t * 2)
a3 = b3 = c3 = UnivariateSpline(t, np.sin(t / 20))

args = (a1, a2, a3, b1, b2, b3, c1, c2, c3,)
soln = odeint(f, y0, t, args)

y1 = soln[:, 0]
y2 = soln[:, 1]
y3 = soln[:, 2]


# %%

class Cybernetic_Model(dict):

    def __init__(self, name):
        dict.__init__(self)
        self.name = name

    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError(r"'Model' object has no attribute '%s'" % key)

    def __setattr__(self, key, value):
        self[key] = value

    def check_type(self):
        for key in ['kmax', 'ke', 'alpha', 'beta', 'n_carbon']:
            if type(self[key]) != np.ndarray:
                try:
                    self[key] = np.array(self[key])
                except:
                    print("Warning: type(self[%s]) != np.ndarray, check it " % key)

    def check_shape(self):
        (n_mets, n_path) = self['Matrix_Smz'].shape
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


ecoli_reduced_cb = Cybernetic_Model('CB model for Ecoli reduced matrix ')

# %% <fminsearch test:>
from scipy import optimize
import random


# print(random.randint(0,9))
def f(x, y):
    k = [x[0] * 3, (x[1] ** 2) * y]
    return random.randint(0, 9)


minimum = optimize.fmin(f, [1000, 1], args=(10,), full_output=True)
minimum[0]

# banana = lambda x: 100*(x[1]-x[0]**2)**2+(1-x[0])**2
# xopt = scipy.optimize.fmin(func=banana, x0=[-1.2,1])

# %%
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("TkAgg")

plt.ion()  # 开启interactive mode 成功的关键函数
fig, ax = plt.subplots(1)
plt.figure(1)

plt.plot([1, 2, 3], [4, 5, 5], '-r')

for i in range(200):
    fig.clf(ax)  # 清空画布上的所有内容
    t_now = i * 0.1
    t.append(t_now)  # 模拟数据增量流入，保存历史数据
    m.append(sin(t_now))  # 模拟数据增量流入，保存历史数据
    ax.plot(t, m, '-r')
    plt.pause(0.01)

# %%
import matplotlib.pyplot as plt
import numpy as np
import time
from math import *

plt.ion()  # 开启interactive mode 成功的关键函数
plt.figure(1)
t = [0]
t_now = 0
m = [sin(t_now)]

for i in range(2000):
    # plt.clf() # 清空画布上的所有内容。此处不能调用此函数，不然之前画出的点，将会被清空。
    t_now = i * 0.1
    """
    由于第次只画一个点，所以此处有两种方式，第一种plot函数中的样式选
    为点'.'、'o'、'*'都可以，就是不能为线段'-'。因为一条线段需要两
    个点才能确定。第二种方法是scatter函数，也即画点。
    """
    plt.plot(t_now, sin(t_now), '.')  # 第次对画布添加一个点，覆盖式的。
    # plt.scatter(t_now, sin(t_now))

    plt.draw()  # 注意此函数需要调用
    time.sleep(0.01)

# %%
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("TkAgg")
import numpy as np
import time
from math import *

plt.ion()  # 开启interactive mode 成功的关键函数
plt.figure(1)
t = np.linspace(0, 20, 100)

for i in range(20):
    # plt.clf() # 清空画布上的所有内容。此处不能调用此函数，不然之前画出的轨迹，将会被清空。
    y = np.sin(t * i / 10.0)
    plt.plot(t, y)  # 一条轨迹
    plt.draw()  # 注意此函数需要调用
    time.sleep(1)

# %%
import matplotlib

gui_env = ['TKAgg', 'GTKAgg', 'Qt4Agg', 'WXAgg']
import matplotlib.pyplot as plt

for env in gui_env:
    try:
        matplotlib.use(i)
        fig, ax = plt.subplots()
        y1 = []
        line = ax.plot(y1, label='test')
        for i in range(50):
            y1.append(i)
            # line[0].remove()
            line = ax.plot(y1, label='test')
            ax.legend()
            # plt.pause(1e-9)
            plt.show()
        print(env)
    except:
        pass

    #
# 其中y1是数据的Y值，只要不停地更y1的数组内容，就可以0.3S刷新一次
# %%
import matplotlib

gui_env = ['TKAgg', 'GTKAgg', 'Qt4Agg', 'WXAgg']
for gui in gui_env:
    try:
        print("testing", gui)
        matplotlib.use(gui, warn=False, force=True)
        from matplotlib import pyplot as plt

        break
    except:
        continue
print("Using:", matplotlib.get_backend())
