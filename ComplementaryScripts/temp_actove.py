#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 9/20/19

"""temp_actove.py
:description : script
:param : 
:returns: 
:rtype: 
"""


# %% <ConvexHull cover_experiment>

from scipy.optimize import lsq_linear
import lsqlin


# Input arguments:
#             C   is m x n dense or sparse matrix
#             d   is n x 1 dense matrix
#             reg is regularization parameter
#             A   is p x n dense or sparse matrix
#             b   is p x 1 dense matrix
#             Aeq is q x n dense or sparse matrix
#             beq is q x 1 dense matrix
#             lb  is n x 1 matrix or scalar
#             ub  is n x 1 matrix or scalar


# out side of the hull

# C, d, reg=0, A=None, b=None, Aeq=None, beq=None, \
#         lb=None, ub=None, x0=None, opts=None
#%%
import numpy as np

C = hull_cutoff_points.T

d = np.array([0.6, 4.2])
d = np.array([1.6, 1.9])
lb = np.zeros(C.shape[1])
ub = np.ones(C.shape[1])
Aeq = np.ones((1,C.shape[1]))
beq = np.ones((1,1))
A = None;
b = None;
x0 = None


#not normalliseed
# res = lsq_linear(C, d)
# print(res.x)

ret = lsqlin.lsqlin(C, d, 0, None, None, Aeq, beq, \
         lb, ub, None, {'show_progress': False})
print (ret['x'].T)

#%% normalize C and d

C0 = C
d0 = d
temp_dt = np.tile(d0,(C.shape[1],1))
C = C*(1/temp_dt.T)
d = np.ones((1,C.shape[0]))[0]

# res = lsq_linear(C, d)
# print(res.x)

ret = lsqlin.lsqlin(C, d, 0, None, None, Aeq, beq, \
         lb, ub, None, {'show_progress': False})
print (ret['x'].T)



# %% <inside of the hull>

C = np.eye(C.shape[1])
d = np.zeros((1,C.shape[1]))[0]
Aeq = np.ones((1,C.shape[1]))
beq = np.ones((1,1))
Aeq = np.vstack([Aeq, C0])
d1 = d0.reshape(Aeq.shape[0]-1,1)
beq = np.vstack([beq, d1])

# res = lsq_linear(C, d)
# print(res.x)

ret = lsqlin.lsqlin(C, d, 0, None, None, Aeq, beq, \
          lb, ub, None, {'show_progress': False})
print (ret['x'].T)