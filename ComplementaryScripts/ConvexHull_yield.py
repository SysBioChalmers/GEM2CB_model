#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-09-17

"""ConvexHull_yield.py
:description : script
:param : 
:returns: 
:rtype: 
"""
# %% ConvexHull base

#hull = ConvexHull(points)
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D #<-- Note the capitalization!
import math
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import numpy as np
from itertools import combinations
import os


os.chdir('../ComplementaryData/')


points = np.array([[0.04471453, 0.60130876],[0.88212444, 0.38443542],
       [0.8032823 , 0.77469327],
       [0.20236158, 0.0311318 ],
       [0.34418886, 0.73012773],
       [0.26171132, 0.49639897],
       [0.57424653, 0.67457535],
       [0.5455598 , 0.24871044],
       [0.52914157, 0.36050058],
       [0.32934503, 0.02044649],
       [0.15627265, 0.68306348],
       [0.87374454, 0.16020225],
       [0.65344122, 0.92569922],
       [0.44822616, 0.22729576],
       [0.36277612, 0.17726654],
       [0.92333561, 0.49174639],
       [0.14659749, 0.48762022],

       [0.95701288, 0.83669329],
       [0.78608686, 0.93377683],
       [0.45363602, 0.3757689 ],
       [0.72192405, 0.66718681],
       [0.44475373, 0.97439096],
       [0.10256349, 0.26267084],
       [0.86974979, 0.85784314],
       [0.65280078, 0.89913258],
       [0.9637922 , 0.90677449],
       [0.32237504, 0.45332723],
       [0.71291365, 0.70448372],
       [0.50491881, 0.70261496],
       [0.66722671, 0.85257059]])

points  = np.array([[0,0],[2.1,1.5],[2,2],[2,0],[0,2],[1.5,1.5]])

# import scipy.io as sio
# A = sio.loadmat('../../../MATLAB/aumic-master/MYA_test_2/A_2d_rand.mat')
# oct_a = A['A']
points =np.genfromtxt("../../../MATLAB/aumic-master/MYA_test_2/A_2d_rand.txt", delimiter=",")

#points = np.genfromtxt("model_ymas_ymas.csv", delimiter=",")
#points = np.random.rand(30, 4)
#points = np.around(points, decimals=2)

# %% <ConvexHull base>
#from pyhull.convex_hull import ConvexHull
qhull_options = ''      #'QG0'mean expect point 0(index)

hull_all = ConvexHull(points,qhull_options)
ndim = hull_all.ndim
hull_all_index = hull_all.vertices
hull_all_points = points[hull_all_index,:]
volume = hull_all.volume

# hull_all_index = set()
# for i in hull.vertices:
#     hull_all_index = hull_all_index | set(i)



# %% <ConvexHull 90%> <option 1 >

persent = 0.99

cutoff_v = persent*hull_all.volume


#%% <option 1 > base points : max in each dim and max distant one: max in each dim

def get_hull_cutoff(hull_all,persent = 0.99,qhull_options = '',options = 1):

    cutoff_v = persent*hull_all.volume

    if cutoff_v>=1:
        hull_cutoff = hull_all
        return hull_cutoff

    ndim = hull_all.ndim
    hull_all_index = hull_all.vertices
    hull_all_points = points[hull_all_index,:]

    if options == 1 :       # find ndim ymax and one max distant as base comb; add points one by one
        hull_cutoff_index_set = set()

        for i in range(0,ndim):     #get yn max
            temp_index = hull_all_points.T[i,:].argmax()
            index = hull_all_index[temp_index]
            hull_cutoff_index_set.add(index)

        while len(hull_cutoff_index_set) < ndim+1:      # get max distant one
            print('base point is not enough')
            add_list = set(hull_all_index) - hull_cutoff_index_set
            distant = 0
            for i in add_list:
                points[list(hull_cutoff_index_set)]
                temp_distant = get_distance(points[list(hull_cutoff_index_set)],points[[i]])
                if temp_distant > distant:
                    distant = temp_distant
                    max_dist_i = i

            hull_cutoff_index_set.add(max_dist_i)
            add_list.remove(max_dist_i)

    pass

hull_cutoff_index_set = set()
for i in range(0,ndim):
    temp_index = hull_all_points.T[i,:].argmax()
    index = hull_all_index[temp_index]
    hull_cutoff_index_set.add(index)





def get_distance(V, p):

    #V = np.eye(3)
    #p = np.array([[1,2,3]])
    p = np.mat(p)
    V = np.mat(V)
    v0=V[:,-1]
    v=V[:,0:V.shape[1]-1]
    (mv, nv)=v.shape
    A = np.zeros((nv,nv))
    b = np.zeros((nv,1))
    for k in range(0,nv):

        for i in range(0,nv):
            A[k,i]=0
            for j in range(0,mv):
                A[k,i]=A[k,i]+(v[j,k]-v0[j])*(v[j,i]-v0[j])
        b[k,0]=0
        for j in range(0,mv):
            b[k,0]=b[k,0]+(v[j,k]-v0[j])*(p[0,j]-v0[j])


    x=np.linalg.inv(A)@b
    sum=np.zeros((mv,1))
    for i in range(0,nv):
        sum=sum + float(x[i])*(v[:,i]-v0)

    dist=np.linalg.norm(v0-p+sum)
    return dist


hull_cutoff_points = points[list(hull_cutoff_index_set),:]

current_v = 0
while current_v < cutoff_v:     #add point one by one
    hull_temp = ConvexHull(points[list(hull_cutoff_index_set),:],qhull_options='Qt QJ Pp Qw Qx')
    temp_current_v = hull_temp.volume

    for i in add_list:
        #hull_temp =hull_temp.add_points(points[i],False)
        temp_selected_set = hull_cutoff_index_set | set([i])
        hull_temp = ConvexHull(points[list(temp_selected_set),:],qhull_options='Qt QJ Pp Qw Qx')
        if hull_temp.volume > temp_current_v:
            #print(i)
            temp_current_v = hull_temp.volume
            temp_i = i
    hull_cutoff_index_set.add(temp_i)
    #add_list = set(hull_all_index) - hull_cutoff_index_set
    add_list.remove(temp_i)
    current_v = temp_current_v

hull_cutoff_index = np.array(hull_cutoff_index_set)
hull_cutoff_points = points[list(hull_cutoff_index_set), :]



# %% <ConvexHull 90%> <option 2 >



current_v = 0
combs = list(combinations(hull_all_index, ndim+1))
for comb in combs:  # frist selection
    #base_points = np.array([points[i] for i in comb])
    base_points = points[list(comb),:]
    try:
        hull_temp = ConvexHull(base_points,qhull_options=['Qt','QJ','Pp','Qw','Qx'])
        # print(base_points)
        # print(hull_temp.volume)
    except:
        print(base_points)
        print('1 dim, line!')
    if hull_temp.volume > current_v:
        current_v = hull_temp.volume
        hull_cutoff_index_comb = comb

hull_cutoff_index_set = set(hull_cutoff_index_comb)
add_list = set(hull_all_index) - hull_cutoff_index_set

hull_cutoff_points = points[list(hull_cutoff_index_set),:]

while current_v < cutoff_v:
    hull_temp = ConvexHull(points[list(hull_cutoff_index_set),:],qhull_options='Qt QJ Pp Qw Qx')
    temp_current_v = hull_temp.volume

    for i in add_list:
        #hull_temp =hull_temp.add_points(points[i],False)
        temp_selected_set = hull_cutoff_index_set | set([i])
        hull_temp = ConvexHull(points[list(temp_selected_set),:],qhull_options='Qt QJ Pp Qw Qx')
        if hull_temp.volume > temp_current_v:
            #print(i)
            temp_current_v = hull_temp.volume
            temp_i = i
    hull_cutoff_index_set.add(temp_i)
    #add_list = set(hull_all_index) - hull_cutoff_index_set
    add_list.remove(temp_i)
    current_v = temp_current_v

hull_cutoff_index = np.array(hull_cutoff_index_set)
hull_cutoff_points = points[list(hull_cutoff_index_set), :]



# %% <ConvexHull 90%> <option 3 >



current_v = 0

for i in range(ndim+1,len(hull_all_index)+1):
    combs = list(combinations(hull_all_index, i))
    for comb in combs:  # frist selection
        base_points = points[list(comb),:]
        try:
            hull_temp = ConvexHull(base_points,qhull_options=['Qt','QJ','Pp','Qw','Qx'])
        except:
            print(base_points)
            print('1 dim, line!')
        if hull_temp.volume > cutoff_v and hull_temp.volume > current_v:
            current_v = hull_temp.volume
            hull_cutoff_index_comb = comb
    if current_v >= cutoff_v:
        break

hull_cutoff_index = np.array(hull_cutoff_index_comb)
hull_cutoff_points = points[list(hull_cutoff_index_comb), :]


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














# %%plot

if hull.ndim ==2:
    fig = plt.figure()
    plt.plot(points[:,0], points[:,1], 'o')
    for simplex in hull.simplices:
         plt.plot(points[simplex, 0], points[simplex, 1], '-',color = 'grey')

    hull_cutoff = ConvexHull(hull_cutoff_points,qhull_options='Qt QJ Pp Qw Qx')
    for simplex in hull_cutoff.simplices:
         plt.plot(hull_cutoff_points[simplex, 0], hull_cutoff_points[simplex, 1], 'r-.')

    # hull_cutoff2 = ConvexHull(hull_cutoff_points2,qhull_options='Qt QJ Pp Qw Qx')
    # for simplex in hull_cutoff2.simplices:
    #      plt.plot(hull_cutoff_points2[simplex, 0], hull_cutoff_points2[simplex, 1], 'g--')


elif hull.ndim ==3:
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot(points[:,0], points[:,1], points[:,2], 'o')
    for simplex in hull.simplices:
        ax.plot(points[simplex, 0], points[simplex, 1], points[simplex, 2], '-',color = 'grey')

    hull_cutoff = ConvexHull(hull_cutoff_points,qhull_options='Qt QJ Pp Qw Qx')
    for simplex in hull_cutoff.simplices:
         plt.plot(hull_cutoff_points[simplex, 0], hull_cutoff_points[simplex, 1],hull_cutoff_points[simplex, 2], 'r-.')

    # hull_cutoff2 = ConvexHull(hull_cutoff_points2,qhull_options='Qt QJ Pp Qw Qx')
    # for simplex in hull_cutoff2.simplices:
    #      plt.plot(hull_cutoff_points2[simplex, 0], hull_cutoff_points2[simplex, 1],hull_cutoff_points2[simplex, 2], 'g--')

plt.show()




