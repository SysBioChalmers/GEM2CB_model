#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-09-17

"""ConvexHull_yield.py
:description : script to caculate the edge pathways, the fondmental method is ConvexHull
:param : all matrix(points) and experiment data cutoff persent.
:returns: selected points(pathways)
:rtype: 
"""

import os
import warnings
from itertools import combinations

import cvxpy as cp
import lsqlin
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from scipy.spatial import ConvexHull
from scipy.spatial import Delaunay


#%%

def get_distance(V, p):
    '''
    function to get the distance of some points and the point
    for get_hull_cutoff function
    :param V: some points set
    :param p: the point
    :return: the distance
    '''

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

def get_hull_cutoff(all_points, hull_all_index, cutoff_v, qhull_options='Qt QJ Pp Qw Qx', method=1):
    '''
    function to get the get_hull_cutoff
    get the minimum number of point and the hull.volume >= cutoff_v

    :param all_points: points set
    :param hull_all_index: the index of all points (most case it will could complete hull set)
    :param cutoff_v: cutoff hull volume , keep new hull.volume >= cutoff_v
    :param qhull_options: ConvexHull function parameters see http://www.qhull.org/
    :param method: our method is 1 efficient and defult the reference method is 2 but all method should get the result
    :return: the indexes of points, get the minimum number of point and the hull.volume >= cutoff_v
    '''

    if len(hull_all_index) == 0:
        hull_all = ConvexHull(all_points, qhull_options=qhull_options)
        hull_all_index = hull_all.vertices

    hull_all_points = all_points[hull_all_index,:]
    ndim = hull_all_points.shape[1]

    if method == 1:  # remove all_points one by one
        hull_cutoff_index_set = set(hull_all_index)  # all points

        while True:
            current_v = 0
            for i in hull_cutoff_index_set:
                temp_points_index = list(hull_cutoff_index_set - set([i]))  # remove all_points one by one
                try:
                    hull_temp = ConvexHull(all_points[temp_points_index, :], qhull_options=qhull_options)

                    if hull_temp.volume > current_v:  # find the max v (get the point index i that affects the volume minimum)
                        current_v = hull_temp.volume
                        temp_i = i
                except:
                    continue

            if current_v < cutoff_v:  # untill the current_v < cutoff_v
                return list(hull_cutoff_index_set)

            else:
                try:
                    hull_cutoff_index_set.remove(temp_i)
                except:
                    return list(hull_cutoff_index_set)

    elif method == 2:  # find ndim ymax points and one max distant point, as base comb;and then add points one by one (the same as reference)

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
        hull_cutoff = ConvexHull(hull_cutoff_points,qhull_options = qhull_options)
        current_v = hull_cutoff.volume

        while current_v < cutoff_v:     #add point one by one (max volume)
            # hull_temp = ConvexHull(all_points[list(hull_cutoff_index_set),:],qhull_options='Qt QJ Pp Qw Qx')
            # temp_current_v = hull_temp.volume

            for i in add_list:
                #hull_temp =hull_temp.add_points(all_points[i],False)
                temp_selected_set = hull_cutoff_index_set | set([i])
                hull_temp = ConvexHull(all_points[list(temp_selected_set),:],qhull_options = qhull_options)

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

    elif method == 3:  # always select the bigest point

        current_v = 0

        combs = list(combinations(hull_all_index, ndim+1))
        for comb in combs:      # frist selection max v base all_points
            #base_points = np.array([all_points[i] for i in comb])
            base_points = all_points[list(comb),:]
            try:
                hull_temp = ConvexHull(base_points,qhull_options = qhull_options)
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
                #hull_temp =hull_temp.add_points(all_points[i],False)
                temp_selected_set = hull_cutoff_index_set | set([i])
                hull_temp = ConvexHull(all_points[list(temp_selected_set),:],qhull_options = qhull_options)

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

    elif method == 4:  # all comb for ndim+1 points to all len(all_points) points
        current_v = 0

        for i in range(ndim+1,len(hull_all_index)+1):
            combs = list(combinations(hull_all_index, i))
            for comb in combs:  # frist selection
                base_points = all_points[list(comb),:]
                try:
                    hull_temp = ConvexHull(base_points,qhull_options = qhull_options)
                except:
                    print(base_points)
                    print('1 dim, line!')
                if hull_temp.volume > cutoff_v:
                    return list(set(comb))


def point_in_hull(point, hull, tolerance=1e-6):
    # Judge whether the point is inside or outside of the hull
    # TODO: more test
    return all(
        (np.dot(eq[:-1], point) + eq[-1] <= tolerance)
        for eq in hull.equations)


def point_in_hull_passed(point, hull, tolerance=1e-12):
    # Judge whether the point is inside or outside of the hull
    p = point
    hull_points = hull.points
    if not isinstance(hull_points, Delaunay):
        hull_ = Delaunay(hull_points)
    return hull_.find_simplex(p) >= 0


def get_hull_active(all_points_,d_,hull_cutoff_index ,qhull_options='Qt QJ Pp Qw Qx', normalize = True):
    '''
    TODO find points without combnation, from the Constrained linear least squares
    math:
    if experiment data in the hull :
        min the sum of weight^2
        min 1/2 ||h||^2_2
        Zyh - ym  = 0,
        h>0,

    elif the data not in the hull :
        min the distant the error
        min 1/2 ||Zyh - ym||^2_2
    reference : https://doi.org/10.1002/bit.22062
    :param all_points_: points set use how many dimations depend on the experiment data ['']
    :param d_: experimen data, the data could be imcomplete use '' replace [1,2,3,'',5]
    :param hull_cutoff_index: the indexes of points to decide which points will be use
    :param qhull_options: ConvexHull function parameters see http://www.qhull.org/
    :param normalize:
    :return:
    '''

    if len(hull_cutoff_index) == 0:  # data check
        hull_cutoff_index = np.arange(0,all_points_.shape[0])
    if len(d_) != all_points_.shape[1]:
        warnings.warn('the experiment data dimation is different with points dimation:\t')
        print(len(d_), all_points_.shape[1])
        return hull_cutoff_index, [], [], []

    index_experiment = []  # experiment could be not complated
    index_empty = []
    d = []
    for i in range(0,len(d_)):
        if d_[i] !='':
            d.append(d_[i])
            index_experiment.append(i)
        else:
            index_empty.append(i)

    all_points = all_points_[:,index_experiment]

    hull_cutoff = ConvexHull(all_points[hull_cutoff_index, :], qhull_options=qhull_options)

    # in_hull = point_in_hull(d, hull_cutoff, tolerance=1e-12)  # in or out of the hull
    # print('point in the hull:\t', in_hull)
    in_hull = False

    C = all_points[hull_cutoff_index,:].T
    d = np.array(d)
    C0 = C
    d0 = d
    # set the pramater of lsqlin
    # ret = lsqlin.lsqlin(C_in, d_in, 0, None, None, Aeq, beq, \
    #         lb=None, ub=None, x0=None, opts=None
    (mC,nC) = C.shape

    lb = np.zeros(nC)
    ub = np.ones(nC)
    # A = None
    # b = None
    # x0 = None
    Aeq = np.ones((1, nC))
    beq = np.ones((1,))
    x = cp.Variable(nC)
    objective = cp.Minimize(cp.sum_squares(C0 * x - d0))
    constraints = [x <= ub,
                   x >= lb,
                   Aeq @ x == beq,
                   ]
    prob = cp.Problem(objective, constraints)
    result = prob.solve(solver=cp.CPLEX)
    weights = x.value
    # print(weights)

    if len(np.where(weights > 1e-6)[0]) > mC + 1 and sum((C0 @ x.value - d0) ** 2) < 0.001:
        in_hull = True
        Aeq = np.vstack([Aeq, C0])
        beq = np.vstack([beq, d0.reshape(mC, 1)])
        beq = beq.reshape(mC + 1)
        # C = np.eye(nC) #np.array([mC]*nC)
        # d = np.ones(nC)
        x = cp.Variable(nC)
        objective = cp.Minimize(cp.norm1(x))
        constraints = [x <= ub,
                       x >= lb,
                       Aeq @ x == beq * 1.0,
                       ]  # Aeq @ x >= beq * 1.0
        prob = cp.Problem(objective, constraints)
        result = prob.solve(solver=cp.CPLEX)
        # print(prob.status)
        if prob.status == 'optimal':
            weights = x.value
        else:
            constraints = [x <= ub,
                           x >= lb,
                           Aeq @ x <= beq * 1.05,
                           Aeq @ x >= beq * 0.95]
            prob = cp.Problem(objective, constraints)
            result = prob.solve(solver=cp.CPLEX)
            if prob.status == 'optimal':
                print('give 5% error')
                weights = x.value
            else:
                in_hull = False

    # normalize C and d
    # C0 = C      # initial C and d
    # d0 = np.array(d)

    # temp_dt = np.tile(d0,(nC,1))
    # C = C0*(1/temp_dt.T)
    # d = np.ones((1,mC))[0]
    #
    # if in_hull:
    #
    #     C_in = np.eye(nC)       #one C and zero d
    #     d_in = np.zeros((1,nC))[0]
    #
    #     if normalize:
    #         # nored
    #         Aeq = np.ones((1, nC))
    #         #beq = np.ones((1, 1))
    #         Aeq = np.vstack([Aeq, C])
    #         beq = np.ones((Aeq.shape[0]))
    #
    #         # ret = lsqlin.lsqlin(C_in, d_in, 0, None, None, Aeq, beq, \
    #         #           lb, ub, None, {'show_progress': False})
    #
    #         x = cp.Variable(nC)
    #         objective = cp.Minimize(cp.sum_squares(C_in*x - d_in ))
    #         constraints = [x <= ub, x >= lb, Aeq @ x == beq]
    #         prob = cp.Problem(objective, constraints)
    #         result = prob.solve()
    #
    #         #print (ret['x'].T)
    #         weights = x.value  #ret['x'].T
    #         weights = combinations_points(weights, mC, nC, C_in, d_in, Aeq, beq, in_hull)
    #
    #     else:
    #         # not nored the same as nored
    #         Aeq0 = np.ones((1, nC))
    #         beq0 = np.ones((1))
    #         Aeq0 = np.vstack([Aeq0, C0])
    #         beq0 = np.vstack([beq0, d0.reshape(Aeq0.shape[0]-1,1)])
    #
    #         # ret0 = lsqlin.lsqlin(C_in, d_in, 0, None, None, Aeq0, beq0, \
    #         #      lb, ub, None, {'show_progress': False})
    #
    #         x = cp.Variable(nC)
    #         objective = cp.Minimize(cp.sum_squares(C_in*x - d_in ))
    #         constraints = [x <= ub, x >= lb, Aeq0 @ x == beq0]
    #         prob = cp.Problem(objective, constraints)
    #         result = prob.solve()
    #         #print (ret0['x'].T)
    #         weights = x.value  #ret['x'].T
    #         weights = combinations_points(weights, mC, nC, C_in, d_in, Aeq0, beq0, in_hull)
    #
    # if  not in_hull:
    #     # Bug: why there are too many pathways !!!!
    #     Aeq = np.ones((1,nC))
    #     beq = np.ones((Aeq.shape[0]))
    #
    #     if normalize:
    #         # ret = lsqlin.lsqlin(C, d, 0, None, None, Aeq, beq, \
    #         #      lb, ub, None, {'show_progress': False})
    #         #print (ret['x'].T)
    #         x = cp.Variable(nC)
    #         objective = cp.Minimize(cp.sum_squares(C*x - d ))
    #         constraints = [x <= ub, x >= lb, Aeq @ x == beq]
    #         prob = cp.Problem(objective, constraints)
    #         result = prob.solve()
    #
    #         weights = x.value  #ret['x'].T
    #         weights = combinations_points(weights, mC, nC, C, d, Aeq, beq, in_hull)
    #
    #     else:
    #         # ret0 = lsqlin.lsqlin(C0, d0, 0, None, None, Aeq, beq, \
    #         #      lb, ub, None, {'show_progress': False})
    #         #print (ret0['x'].T)
    #
    #         x = cp.Variable(nC)
    #         objective = cp.Minimize(cp.sum_squares(C0*x - d0 ))
    #         constraints = [x <= ub, x >= lb, Aeq @ x == beq]
    #         prob = cp.Problem(objective, constraints)
    #         result = prob.solve()
    #
    #         weights = x.value  #ret['x'].T
    #         weights = combinations_points(weights, mC, nC, C0, d0, Aeq, beq, in_hull)

    weights = list(weights)
    weights = np.around(weights, decimals=6)

    hull_active_index = np.array(hull_cutoff_index)[[i for i in range(0, len(weights)) if weights[i] > 1e-6]]
    estimated_data = C0@weights

    estimated_data_ = list(estimated_data[:])
    for index in index_empty:
        estimated_data_.insert(index, '')

    print('Point in or not in hull: ', in_hull)
    print('weights =\t', weights)
    print('estimated_data vs experiment_data: \t')
    print(np.around(estimated_data_, decimals=6))
    print(np.around(d_, decimals=6))

    return hull_active_index,weights,estimated_data,in_hull


def combinations_points(weights, mC, nC, C, d, Aeq, beq, in_hull):
    '''
    scrips to reduce the points
    :param weights:
    :param mC:
    :param nC:
    :param C:
    :param d:
    :param Aeq:
    :param beq:
    :param in_hull:
    :return:
    '''
    # (mC,nC) = C.shape
    nweights = 0  # check how many points
    for i in weights:
        if i > 1e-6:
            nweights = nweights + 1
    if in_hull:
        dim = mC + 1
    else:
        dim = mC
    # print(nweights,'..............',mC)

    if nweights > dim:  # if the number of points > dim+1 ,reduce
        print('Looking for points by combinations...')
        lb_ = np.zeros(dim)
        ub_ = np.ones(dim)
        min_dis_ave_2 = np.inf
        comb_list = list(combinations(np.arange(nC), dim))
        if len(comb_list) > 1e9:
            print('Too many combinations....try reduce the points!!!')
            return weights

        for i in comb_list:
            Aeq_ = Aeq[:, i]
            beq_ = np.ones((Aeq_.shape[0], 1))
            C_part = C[:, i]
            ret = lsqlin.lsqlin(C_part, d, 0, None, None, Aeq_, beq_, \
                                lb_, ub_, None, {'show_progress': False})

            if ret['gap'] < 1e-6:
                # TODO check the cutoff and sum((C_part @ ret['x'] - d.reshape(, 1)) ** 2) < 0.1 0.1 ???
                # print(sum((C_part @ ret['x'] - d.reshape(3, 1)) ** 2))
                dis_ave = np.array(ret['x']) - np.array([1 / (dim)])
                dis_ave_2 = dis_ave ** 2
                if sum(dis_ave_2) < min_dis_ave_2:
                    min_dis_ave_2 = sum(dis_ave_2)
                    min_i = i
                    weights_i = ret['x']
        weights_ = np.zeros(nC)
        weights_[list(min_i)] = list(weights_i)
    else:
        weights_ = weights
    return weights_


def hull_plot(all_points,hulls = [],labels = [], markers = [],colors = [],alphas = []):
    #hulls = [hull_1,]
    #color='green', marker='o', linestyle='dashed'

    if len(labels) == 0:
        labels = ['all_points','hull_all','hull_cutoff','hull_active','data_5','data_6']
    if len(markers) == 0:
        markers = ['o','-','-.','--','-.','--']
    if len(colors) == 0:
        colors = ['tab:blue','grey','red','green','blue','yellow']
    if len(alphas) == 0:
        alphas = [0.8,0.8,0.8,0.8,0.8,0.8]

    ndim = all_points.shape[1]
    if ndim ==2:
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        lines_0 = ax1.plot(all_points[:,0], all_points[:,1], 'o',label = labels[0],color = colors[0],alpha = alphas[0])
        #markersize=8

        for index in range(0,len(hulls)):
            hull = hulls[index]

            if str(type(hull))=="<class 'numpy.ndarray'>":
                lines_n = ax1.plot(hull[:,0], hull[:,1], "D",label = labels[index+1], color = colors[index+1], \
                                   alpha = alphas[index+1])

            else:
                for simplex in hull.simplices:
                    lines_n = ax1.plot(hull.points[simplex, 0], hull.points[simplex, 1], markers[index+1], \
                                       color = colors[index+1],alpha = alphas[index+1])
                lines_n[0].set_label(labels[index+1])

    elif ndim ==3:
        fig = plt.figure()
        ax1 = fig.gca(projection='3d')
        ax1.plot(all_points[:,0], all_points[:,1], all_points[:,2], 'o',label = labels[0],color = colors[0],\
                 alpha = alphas[0])

        for index in range(0,len(hulls)):
            hull = hulls[index]

            if str(type(hull))=="<class 'numpy.ndarray'>":
                ax1.plot(hull[:,0], hull[:,1],  hull[:,2], "D", label = labels[index+1], color = colors[index+1], \
                         alpha = alphas[index+1])
            else:
                for simplex in hull.simplices:
                    lines_n = ax1.plot(hull.points[simplex, 0], hull.points[simplex, 1], hull.points[simplex, 2], \
                                       markers[index],color = colors[index],alpha = alphas[index])
                    lines_n[0].set_label(labels[index+1])
    return fig,ax1


def pipeline_mya(all_points, experiment_datas=[], qhull_options='Qt QJ Pp Qw Qx', cutoff_persent=0.99, method=1,
                 normalize=True):
    # indexes, hulls , weightss , estimated_datas, in_hulls = pipeline_mya(all_points, experiment_datas = [], qhull_options = 'Qt QJ Pp Qw Qx', cutoff_persent = 0.99)

    #  < Setp1 ConvexHull base >
    print('\n---------- ConvexHull base ... ---------- ')

    hull_all = ConvexHull(all_points, qhull_options=qhull_options)
    ndim = hull_all.ndim
    hull_all_index = hull_all.vertices
    print('hull_all_index = \t', list(hull_all_index), '\n len:\t', len(list(hull_all_index)))

    #  < Setp2 ConvexHull cutoff >
    print('\n---------- ConvexHull cutoff ... ---------- ')

    if cutoff_persent >= 1:
        hull_cutoff_index = hull_all_index
        # hull_cutoff = hull_all
        print('cutoff_persent >= 1')

    else:
        cutoff_v = cutoff_persent * hull_all.volume
        hull_cutoff_index = get_hull_cutoff(all_points, hull_all_index, cutoff_v, qhull_options=qhull_options,
                                            method=method)
        # hull_cutoff = ConvexHull(all_points[hull_cutoff_index, :], qhull_options=qhull_options)

    print('hull_cutoff_index = \t', hull_cutoff_index, '\n len:\t', len(list(hull_cutoff_index)))

    #  < Setp3 ConvexHull active >
    print('\n---------- ConvexHull active ... ---------- ')

    hull_active_indexes = []
    weightss = []
    estimated_datas = []
    in_hulls = []
    hull_cutoff_actives = []

    if len(experiment_datas) == 0:
        print('No eperiment data')
        hull_active_indexes = [hull_cutoff_index]

    else:  #
        for experiment_data in experiment_datas:
            print('\n\t\t\t\t---------- Experimrnt data i ... ---------- ')
            hull_cutoff_index_ = hull_cutoff_index[:]

            hull_active_index, weights, estimated_data, in_hull = get_hull_active(all_points, experiment_data, \
                                                                                  hull_cutoff_index_,
                                                                                  qhull_options=qhull_options,
                                                                                  normalize=normalize)
            print('hull_active_index = \t', list(hull_active_index), '\n \t\t len:\t', len(list(hull_active_index)),
                  '\n')

            # if in_hull:
            #     hull_cutoff_active = ConvexHull(all_points[hull_active_index, :], qhull_options)
            # else:
            #     hull_cutoff_active = all_points[hull_active_index, :]

            hull_active_indexes.append(hull_active_index)
            weightss.append(weights)
            estimated_datas.append(estimated_data)
            in_hulls.append(in_hull)
            # hull_cutoff_actives.append(hull_cutoff_active)

    hulls = ['return removed, please caculate it ']
    if len(experiment_datas) <= 1:
        indexes = [hull_all_index, np.array(hull_cutoff_index), hull_active_indexes[0]]
        # hulls = [hull_all, hull_cutoff, hull_cutoff_actives]
    else:
        indexes = [hull_all_index, np.array(hull_cutoff_index), hull_active_indexes]
        # hulls = [hull_all, hull_cutoff, hull_cutoff_actives]

    return indexes, weightss, estimated_datas, in_hulls

if __name__ == '__main__':
    # %% <set data and pramaters>   all_points, experiment_data ,qhull_options = 'Qt QJ Pp Qw Qx', cutoff_persent = 0.99
    os.chdir('../ComplementaryData/')
    print('loading data ...')
    points_all = sio.loadmat('../../../MATLAB/aumic-master/MYA_test_2/points_all.mat')
    points_2d = points_all['points_2d']
    points_3d = points_all['points_3d']
    points_4d = points_all['points_4d']
    points_sq = points_all['points_sq']         #points_sq  = np.array([[0,0],[2.1,1.5],[2,2],[2,0],[0,2],[1.5,1.5]])
    points_glc_33 = points_all['points_glc_33']

    # all_points = points_sq
    # d_t = np.array([0.5, 1.9])
    # d_f = np.array([1, 2.5])
    # experiment_datas = [d_t, d_f, d_t]
    # experiment_data = experiment_datas[0]

    all_points = points_4d
    experiment_data = [0.6, 0.6, 0.6, 0.6]

    # all_points = points_glc_33
    # dataValue = np.array([0.0169, 1.8878, 0.0556])
    # experiment_datas = [dataValue]
    # experiment_data = experiment_datas[0]

    # %%
    # all_points =  np.random.rand(200,6)#points_2d  # try different data points
    # experiment_datas = []
    # experiment_data = [0.26112345, 0.68739509, 0.85738014, 0.16804945, 0.40893545,
    #     0.3448415 ]

    cutoff_persent = 0.95  # 0.99
    qhull_options = 'QJ A0.999'

    #  <Setp1 ConvexHull base>

    print('\nConvexHull base ...')
    hull_all = ConvexHull(all_points, qhull_options=qhull_options)
    ndim = hull_all.ndim
    hull_all_index = hull_all.vertices
    print('hull_all_index = ', list(hull_all_index), '\n len:\t', len(list(hull_all_index)))

    # hull_all_points = all_points[hull_all_index,:]

    #  <Setp2 ConvexHull cutoff>
    print('ConvexHull cutoff ...')
    cutoff_v = cutoff_persent * hull_all.volume
    hull_cutoff_index = get_hull_cutoff(all_points, hull_all_index, cutoff_v, qhull_options=qhull_options, method=1)
    print('hull_cutoff_index = ', hull_cutoff_index, '\n len:\t', len(list(hull_cutoff_index)))

    hull_cutoff = ConvexHull(all_points[hull_cutoff_index, :], qhull_options=qhull_options)

    # %% <Setp3 ConvexHull active>
    print('ConvexHull active ...')
    # experiment_data = [0.3, 0.68739509, 0.85738014, 0.16804945, 0.40893545,0.3448415 ]

    hull_active_index, weights, estimated_data, in_hull = get_hull_active(all_points, experiment_data,
                                                                          hull_cutoff_index,
                                                                          qhull_options=qhull_options, normalize=True)
    print(hull_active_index)

    try:
        hull_cutoff_active = ConvexHull(all_points[hull_active_index, :], qhull_options)
    except:
        hull_cutoff_active = all_points[hull_active_index, :]


    # %% < pipeline_mya >
    indexes, weightss, estimated_datas, in_hulls = pipeline_mya(all_points, experiment_datas=[experiment_data],
                                                                cutoff_persent=cutoff_persent,
                                                                qhull_options=qhull_options,
                                                                method=1, normalize=True)

    # %% plot  only for 2d
    if ndim == 2:
        hulls = [hull_all, hull_cutoff, hull_cutoff_active]

        fig, ax1 = hull_plot(all_points, hulls, labels=['all_points', 'hull_all', 'hull_cutoff', 'hull_active'])
        # ax1.plot(experiment_data[0], experiment_data[1], 'r*', label='experiment data')
        ax1.legend(loc='lower left', ncol=3, bbox_to_anchor=(0, -0.3, 1, 0), mode="expand", borderaxespad=0.)
        ax1.set_xlabel('Biomass/ Sourse')
        ax1.set_ylabel('Production_1/ Sourse')
        fig.show()
    # %% < test resuts:>
    '''
    points  experiment  estimated_data  in_hull act_index number
    points_sq:  [0.5, 1.9]  [0.5, 1.9]  True [0 2 4] 3  passed
    points_sq:  [1,3]  [1.0, 2.0]  False [2 4] 2        passed
    
    points_2d:  [0.8,0.6]  [0.8000002661226724, 0.5999998133271094]  True [ 3  8 12] 3  passed
    points_2d:  [1,3]  [0.9344051189612134, 0.9994916200977046]  False [3] 1    passed
    
    points_3d:  [0.8, 0.6, 0.6]  [0.7999998126345972, 0.599999961784013, 0.5999999865075214]  True [ 5 12 23 24] 4  passed
    points_3d:  [1,0.6,0.6]  [0.9190085753931684, 0.591700064282075, 0.6046583530568356]  False [12 23 24] 1    passed
    
    points_glc_33:  [0.0169 1.8878 0.0556]  [0.01662157720811822, 1.464659815092023, 0.055366326919507505]  False [ 1 22 23] 3  passed
    
    
    '''
