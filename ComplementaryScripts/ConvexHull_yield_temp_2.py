#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-09-17

"""ConvexHull_yield.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import numpy as np
from itertools import combinations
import os
import lsqlin
from mpl_toolkits.mplot3d import axes3d, Axes3D #<-- Note the capitalization!
import math
import scipy.io as sio


def get_hull_cutoff(all_points,hull_all_index,ndim,cutoff_v,qhull_options = '',options = 1):

    # return hull_cutoff_index
    hull_all_points = all_points[hull_all_index,:]

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
                temp_points = all_points[list(hull_cutoff_index_set)]
                temp_distant = get_distance(temp_points,all_points[[i]])
                if temp_distant > distant:
                    distant = temp_distant
                    max_dist_i = i

            hull_cutoff_index_set.add(max_dist_i)
            add_list.remove(max_dist_i)

        hull_cutoff_points = all_points[list(hull_cutoff_index_set),:]
        hull_cutoff = ConvexHull(hull_cutoff_points,qhull_options)
        current_v = hull_cutoff.volume

        while current_v < cutoff_v:     #add point one by one (max volume)
            # hull_temp = ConvexHull(points[list(hull_cutoff_index_set),:],qhull_options='Qt QJ Pp Qw Qx')
            # temp_current_v = hull_temp.volume

            for i in add_list:
                #hull_temp =hull_temp.add_points(points[i],False)
                temp_selected_set = hull_cutoff_index_set | set([i])
                hull_temp = ConvexHull(points[list(temp_selected_set),:],qhull_options)

                if hull_temp.volume > current_v:
                    current_v = hull_temp.volume
                    temp_i = i

                if current_v >= cutoff_v:
                    hull_cutoff = hull_temp
                    hull_cutoff_index_set = temp_selected_set
                    #print(hull_cutoff_index_set)
                    return list(hull_cutoff_index_set)

            hull_cutoff_index_set.add(temp_i)       #add the v max point and continue
            add_list.remove(temp_i)

    elif options == 2 :     #always select the bigest point

        current_v = 0

        combs = list(combinations(hull_all_index, ndim+1))
        for comb in combs:      # frist selection max v base points
            #base_points = np.array([points[i] for i in comb])
            base_points = all_points[list(comb),:]
            try:
                hull_temp = ConvexHull(base_points,qhull_options)
                # print(base_points)
                # print(hull_temp.volume)
            except:
                print(base_points)
                print('point in edge')
            if hull_temp.volume > current_v:
                current_v = hull_temp.volume
                hull_cutoff_index_comb = comb

        hull_cutoff_index_set = set(hull_cutoff_index_comb)
        add_list = set(hull_all_index) - hull_cutoff_index_set

        while current_v < cutoff_v:

            for i in add_list:
                #hull_temp =hull_temp.add_points(points[i],False)
                temp_selected_set = hull_cutoff_index_set | set([i])
                hull_temp = ConvexHull(points[list(temp_selected_set),:],qhull_options)

                if hull_temp.volume > current_v:
                    current_v = hull_temp.volume
                    temp_i = i

                # if current_v >= cutoff_v:     # find return
                #     hull_cutoff = hull_temp
                #     hull_cutoff_index_set = temp_selected_set
                #     #print(hull_cutoff_index_set)
                #     return list(hull_cutoff_index_set)

            hull_cutoff_index_set.add(temp_i)       #add the v max point and continue
            add_list.remove(temp_i)
        return list(hull_cutoff_index_set)

    elif options == 3 :         # all comb for ndim+1 points to all len(points) points
        current_v = 0

        for i in range(ndim+1,len(hull_all_index)+1):
            combs = list(combinations(hull_all_index, i))
            for comb in combs:  # frist selection
                base_points = all_points[list(comb),:]
                try:
                    hull_temp = ConvexHull(base_points,qhull_options)
                except:
                    print(base_points)
                    print('1 dim, line!')
                if hull_temp.volume > cutoff_v:
                    return list(set(comb))

    elif options == 4 :         # remove points one by one
        hull_cutoff_index_set = set(hull_all_index)

        while True :
            current_v = 0
            for i in hull_cutoff_index_set:

                temp_points_index = list(hull_cutoff_index_set - set([i]))
                hull_temp = ConvexHull(all_points[temp_points_index,:],qhull_options)

                if hull_temp.volume > current_v :
                    current_v = hull_temp.volume
                    temp_i = i

            if current_v < cutoff_v :
                return list(hull_cutoff_index_set)

            else:
                hull_cutoff_index_set.remove(temp_i)

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

def hull_plot(points,hulls = [],line_type = [],line_color = []):
    #hulls = [hull_1,]
    ndim = points.shape[1]
    if len(line_type) == 0:
        line_type = ['-','-.','--','-.','--']

    if len(line_color) == 0:
        line_color = ['grey','red','green','blue','yellow']

    if ndim ==2:
        # fig = plt.figure()
        plt.plot(points[:,0], points[:,1], 'o')
        for index in range(0,len(hulls)):
            hull = hulls[index]
            for simplex in hull.simplices:
                 plt.plot(hull.points[simplex, 0], hull.points[simplex, 1], line_type[index],color = line_color[index])
        plt.show()


    elif ndim ==3:
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot(points[:,0], points[:,1], points[:,2], 'o')
        for index in range(0,len(hulls)):
            hull = hulls[index]
            for simplex in hull.simplices:
                ax.plot(points[simplex, 0], points[simplex, 1], points[simplex, 2], line_type[index],line_color[index])
        plt.show()

def point_in_hull(point, hull, tolerance=1e-12):
    return all(
        (np.dot(eq[:-1], point) + eq[-1] <= tolerance)
        for eq in hull.equations)

def get_hull_active(C,d,hull_cutoff, normalize = True):

    in_hull = point_in_hull(d, hull_cutoff, tolerance=1e-12)

    # C, d, reg=0, A=None, b=None, Aeq=None, beq=None, \
    #         lb=None, ub=None, x0=None, opts=None
    (mC,nC) = C.shape
    lb = np.zeros(nC)
    ub = np.ones(nC)
    A = None
    b = None
    x0 = None

    # normalize C and d
    C0 = C      # initial C and d
    d0 = d

    temp_dt = np.tile(d0,(nC,1))
    C = C0*(1/temp_dt.T)
    d = np.ones((1,mC))[0]

    if in_hull:

        C_in = np.eye(nC)       #one C and zero d
        d_in = np.zeros((1,nC))[0]

        if normalize:
            # nored
            Aeq = np.ones((1, nC))
            #beq = np.ones((1, 1))
            Aeq = np.vstack([Aeq, C])
            beq = np.ones((Aeq.shape[0],1))

            ret = lsqlin.lsqlin(C_in, d_in, 0, None, None, Aeq, beq, \
                      lb, ub, None, {'show_progress': False})
            print (ret['x'].T)
            weights = ret['x'].T

        else:
            # not nored the same as nored
            Aeq0 = np.ones((1, nC))
            beq0 = np.ones((1, 1))
            Aeq0 = np.vstack([Aeq0, C0])
            beq0 = np.vstack([beq0, d0.reshape(Aeq0.shape[0]-1,1)])

            ret0 = lsqlin.lsqlin(C_in, d_in, 0, None, None, Aeq0, beq0, \
                 lb, ub, None, {'show_progress': False})
            print (ret0['x'].T)
            weights = ret0['x'].T

    if  not in_hull:
        Aeq = np.ones((1,nC))
        beq = np.ones((1,1))

        if normalize:
            ret = lsqlin.lsqlin(C, d, 0, None, None, Aeq, beq, \
                 lb, ub, None, {'show_progress': False})
            print (ret['x'].T)
            weights = ret['x'].T

        else:
            ret0 = lsqlin.lsqlin(C0, d0, 0, None, None, Aeq, beq, \
                 lb, ub, None, {'show_progress': False})
            print (ret0['x'].T)
            weights = ret0['x'].T
    return list(weights)


os.chdir('../ComplementaryData/')

# %% load data

print('load data')
#points = np.random.rand(30, 2)
points  = np.array([[0,0],[2.1,1.5],[2,2],[2,0],[0,2],[1.5,1.5]])

# A = sio.loadmat('../../../MATLAB/aumic-master/MYA_test_2/A_2d_rand.mat')
# oct_a = A['A']
# points =np.genfromtxt("../../../MATLAB/aumic-master/MYA_test_2/A_2d_rand.txt", delimiter=",")

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



# %% <ConvexHull cutoff>

persent = 0.90
cutoff_v = persent*hull_all.volume

hull_cutoff_index_1 = get_hull_cutoff(points,hull_all_index,ndim,cutoff_v,qhull_options = '',options = 1)
print(hull_cutoff_index_1)

hull_cutoff_index_2 = get_hull_cutoff(points,hull_all_index,ndim,cutoff_v,qhull_options = '',options = 2)
print(hull_cutoff_index_2)

hull_cutoff_index_3 = get_hull_cutoff(points,hull_all_index,ndim,cutoff_v,qhull_options = '',options = 3)
print(hull_cutoff_index_3)

hull_cutoff_index_4 = get_hull_cutoff(points,hull_all_index,ndim,cutoff_v,qhull_options = '',options = 4)
print(hull_cutoff_index_4)

hull_cutoff_1 = ConvexHull(points[hull_cutoff_index_1,:],qhull_options)
hull_cutoff_2 = ConvexHull(points[hull_cutoff_index_2,:],qhull_options)
hull_cutoff_3 = ConvexHull(points[hull_cutoff_index_3,:],qhull_options)
hull_cutoff_4 = ConvexHull(points[hull_cutoff_index_4,:],qhull_options)

hulls = [hull_all,hull_cutoff_1,hull_cutoff_2,hull_cutoff_3,hull_cutoff_4]
hull_plot(points, hulls,line_type = [],line_color = [])

# %% <ConvexHull active>



d_t = np.array([0.5,1.9])
d_f = np.array([1,2.5])
hull_cutoff = hull_cutoff_4

C = points[hull_cutoff_index_4,:].T
d = d_t


get_hull_active(C,d,hull_cutoff, normalize = True)



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





