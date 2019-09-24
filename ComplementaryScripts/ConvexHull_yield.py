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
from itertools import combinations

import lsqlin
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from scipy.spatial import ConvexHull


#%%
def get_hull_cutoff(all_points,hull_all_index,cutoff_v,qhull_options = '',options = 1):

    # return hull_cutoff_index
    hull_all_points = all_points[hull_all_index,:]
    ndim = hull_all_points.shape[1]

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

    elif options == 2 :     #always select the bigest point

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

    elif options == 3 :         # all comb for ndim+1 points to all len(all_points) points
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

    elif options == 4 :         # remove all_points one by one
        hull_cutoff_index_set = set(hull_all_index)

        while True :
            current_v = 0
            for i in hull_cutoff_index_set:

                temp_points_index = list(hull_cutoff_index_set - set([i]))
                hull_temp = ConvexHull(all_points[temp_points_index,:],qhull_options = qhull_options)

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


def point_in_hull(point, hull, tolerance=1e-12):
    return all(
        (np.dot(eq[:-1], point) + eq[-1] <= tolerance)
        for eq in hull.equations)


def get_hull_active(all_points,d,hull_cutoff,hull_cutoff_index = [], normalize = True):

    if len(hull_cutoff_index)==0:
        hull_cutoff_index = np.arange(0,all_points.shape[0])

    C = all_points[hull_cutoff_index,:].T

    in_hull = point_in_hull(d, hull_cutoff, tolerance=1e-12)

    # ret = lsqlin.lsqlin(C_in, d_in, 0, None, None, Aeq, beq, \
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
            print('experiment data in hull')
            # nored
            Aeq = np.ones((1, nC))
            #beq = np.ones((1, 1))
            Aeq = np.vstack([Aeq, C])
            beq = np.ones((Aeq.shape[0],1))

            ret = lsqlin.lsqlin(C_in, d_in, 0, None, None, Aeq, beq, \
                      lb, ub, None, {'show_progress': False})
            #print (ret['x'].T)
            weights = ret['x'].T

        else:
            # not nored the same as nored
            Aeq0 = np.ones((1, nC))
            beq0 = np.ones((1, 1))
            Aeq0 = np.vstack([Aeq0, C0])
            beq0 = np.vstack([beq0, d0.reshape(Aeq0.shape[0]-1,1)])

            ret0 = lsqlin.lsqlin(C_in, d_in, 0, None, None, Aeq0, beq0, \
                 lb, ub, None, {'show_progress': False})
            #print (ret0['x'].T)
            weights = ret0['x'].T

    if  not in_hull:
        print('experiment data not in hull')
        Aeq = np.ones((1,nC))
        beq = np.ones((1,1))

        if normalize:
            ret = lsqlin.lsqlin(C, d, 0, None, None, Aeq, beq, \
                 lb, ub, None, {'show_progress': False})
            #print (ret['x'].T)
            weights = ret['x'].T

        else:
            ret0 = lsqlin.lsqlin(C0, d0, 0, None, None, Aeq, beq, \
                 lb, ub, None, {'show_progress': False})
            #print (ret0['x'].T)
            weights = ret0['x'].T

    weights = list(weights)
    weights = np.around(weights, decimals=4)


    hull_active_index = np.array(hull_cutoff_index)[[i for i in range(0,len(weights)) if weights[i] > 1e-4 ]]
    estimated_data = C0@weights

    return hull_active_index,weights,estimated_data,in_hull


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


def pipeline_mya(all_points, experiment_datas = [], qhull_options = 'Qt QJ Pp Qw Qx', cutoff_persent = 0.99):
        #indexes, hulls , weightss , estimated_datas, in_hulls = pipeline_mya(all_points, experiment_datas = [], qhull_options = 'Qt QJ Pp Qw Qx', cutoff_persent = 0.99)

        #  < Setp1 ConvexHull base >
        print('\n---------- ConvexHull base ... ---------- ')

        hull_all = ConvexHull(all_points,qhull_options = qhull_options )
        ndim = hull_all.ndim
        hull_all_index = hull_all.vertices
        print('hull_all_index = \t', hull_all_index)

        #  < Setp2 ConvexHull cutoff >
        print('\n---------- ConvexHull cutoff ... ---------- ')

        if cutoff_persent >=1 :
            hull_cutoff_index = hull_all_index
            hull_cutoff = hull_all
            print('cutoff_persent >= 1')

        else:
            cutoff_v = cutoff_persent*hull_all.volume
            hull_cutoff_index = get_hull_cutoff(all_points,hull_all_index,cutoff_v,qhull_options = '',options = 1)
            hull_cutoff = ConvexHull(all_points[hull_cutoff_index,:],qhull_options = qhull_options)

        print('hull_cutoff_index = \t', hull_cutoff_index)

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

        # elif len(experiment_datas) == 1:
        #     experiment_data = experiment_datas[0]
        #     hull_active_index,weights,estimated_data,in_hull = get_hull_active(all_points,experiment_data,hull_cutoff, \
        #                                                                        hull_cutoff_index, normalize = True)
        #     print(hull_active_index)
        #
        #     # if in_hull:
        #     #     hull_cutoff_active = ConvexHull(all_points[hull_active_index,:],qhull_options)
        #     # else:
        #     #     hull_cutoff_active = all_points[hull_active_index,:]
        else: #

            for experiment_data in experiment_datas :
                #experiment_data = experiment_datas[0]
                hull_active_index,weights,estimated_data,in_hull = get_hull_active(all_points,experiment_data,hull_cutoff, \
                                                                                   hull_cutoff_index, normalize = True)
                print('hull_active_index = \t'  ,hull_active_index )
                print('weights =\t' , weights)
                print('estimated_data vs experiment_data: \t' )
                print(estimated_data)
                print(estimated_data,'\n')

                if in_hull:
                    hull_cutoff_active = ConvexHull(all_points[hull_active_index,:],qhull_options)
                else:
                    hull_cutoff_active = all_points[hull_active_index,:]

                hull_active_indexes.append(hull_active_index)
                weightss.append(weights)
                estimated_datas.append(estimated_data)
                in_hulls.append(in_hull)
                hull_cutoff_actives.append(hull_cutoff_active)

        if len(experiment_datas) <=1:
            indexes = [hull_all_index,hull_cutoff_index,hull_active_indexes[0]]
            hulls = [hull_all,hull_cutoff,hull_cutoff_actives]
        else:
            indexes = [hull_all_index,hull_cutoff_index,hull_active_indexes]
            hulls = [hull_all,hull_cutoff,hull_cutoff_actives]

        return indexes, hulls , weightss , estimated_datas, in_hulls


#%%
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
    dataValue =np.array([0.0169,1.8878,0.0556])
    d_t = np.array([0.5,1.9])
    d_f = np.array([1,2.5])


    all_points = points_glc_33
    experiment_data_1 = dataValue
    experiment_data_2 = np.array([0.01,0.50,0.0556])

    experiment_datas = [experiment_data_1,experiment_data_2]
    qhull_options = 'Qt QJ Pp Qw Qx'      #'QG0'mean expect point 0(index)
    cutoff_persent = 0.99


    indexes, hulls , weightss , estimated_datas, in_hulls = pipeline_mya(all_points, experiment_datas , qhull_options = 'Qt QJ Pp Qw Qx', cutoff_persent = 0.99)



    # # %% <Setp1 ConvexHull base>
    #
    # print('\nConvexHull base ...')
    # hull_all = ConvexHull(all_points,qhull_options = qhull_options )
    # ndim = hull_all.ndim
    # hull_all_index = hull_all.vertices
    # print('hull_all_index = ',hull_all_index)
    #
    # #hull_all_points = all_points[hull_all_index,:]
    #
    # # %% <Setp2 ConvexHull cutoff>
    # print('ConvexHull cutoff ...')
    #
    # cutoff_v = cutoff_persent*hull_all.volume
    # hull_cutoff_index = get_hull_cutoff(all_points,hull_all_index,cutoff_v,qhull_options = '',options = 1)
    # print(hull_cutoff_index)
    # hull_cutoff = ConvexHull(all_points[hull_cutoff_index,:],qhull_options = qhull_options)
    #
    # # %% <Setp3 ConvexHull active>
    # print('ConvexHull active ...')
    #
    # hull_active_index,weights,estimated_data,in_hull = get_hull_active(all_points,experiment_data,hull_cutoff, \
    #                                                                    hull_cutoff_index, normalize = True)
    # print(hull_active_index)
    #
    # if in_hull:
    #     hull_cutoff_active = ConvexHull(all_points[hull_active_index,:],qhull_options)
    # else:
    #     hull_cutoff_active = all_points[hull_active_index,:]
    #
    # # %%
    #
    # hulls = [hull_all,hull_cutoff,hull_cutoff_active]

    # fig,ax1 = hull_plot(all_points, hulls,labels =['all_points','hull_all','hull_cutoff','hull_active'])
    # ax1.plot(experiment_data[0], experiment_data[1], 'r*',label = 'experiment data')
    # ax1.legend(loc='lower left',ncol=3,bbox_to_anchor=(0, -0.3,1,0),mode="expand", borderaxespad=0.)
    # ax1.set_xlabel('Biomass/ Sourse')
    # ax1.set_ylabel('Production_1/ Sourse')
    # fig.show()




