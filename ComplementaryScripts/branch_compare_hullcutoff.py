#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 9/22/19

"""branch_compare_hullcutoff.py
:description : script
:param : 
:returns: 
:rtype: 
"""
# %% <set data and pramaters>   all_points, qhull_options, cutoff_persent
#     os.chdir('../ComplementaryData/')
#     print('loading data ...')
#
#
#     points_all = sio.loadmat('../../../MATLAB/aumic-master/MYA_test_2/points_all.mat')
#     points_2d = points_all['points_2d']
#     points_3d = points_all['points_3d']
#     points_4d = points_all['points_4d']
#     points_sq = points_all['points_sq']
#     points_glc_33 = points_all['points_glc_33']
#     dataValue =np.array([0.0169,1.8878,0.0556])
#     d_t = np.array([0.5,1.9])
#     d_f = np.array([1,2.5])
#     #points_sq  = np.array([[0,0],[2.1,1.5],[2,2],[2,0],[0,2],[1.5,1.5]])
#
#     all_points = points_glc_33
#     experiment_data = dataValue
#
#
#     qhull_options = 'Qt QJ Pp Qw Qx'      #'QG0'mean expect point 0(index)
#     cutoff_persent = 0.99
#
#     # %% <Setp1 ConvexHull base>
#     print('ConvexHull base ...')
#
#     hull_all = ConvexHull(all_points,qhull_options = qhull_options )
#     ndim = hull_all.ndim
#     hull_all_index = hull_all.vertices
#     print(hull_all_index)
#     #hull_all_points = all_points[hull_all_index,:]
#
#     # %% <ConvexHull cutoff>
#     print('ConvexHull cutoff ...')
#
#     cutoff_v = cutoff_persent*hull_all.volume
#     hull_cutoff_index = get_hull_cutoff(all_points,hull_all_index,cutoff_v,qhull_options = '',options = 4)
#     print(hull_cutoff_index)
#     hull_cutoff = ConvexHull(all_points[hull_cutoff_index,:],qhull_options = qhull_options)
#
#
#     # %% <ConvexHull active>
#     print('ConvexHull active ...')
#
#
#
#     compare_options = False     #test
#     if compare_options == True:
#         hull_cutoff_index_1 = get_hull_cutoff(all_points,hull_all_index,cutoff_v,qhull_options = '',options = 1)
#         print(hull_cutoff_index_1)
#
#         hull_cutoff_index_2 = get_hull_cutoff(all_points,hull_all_index,cutoff_v,qhull_options = '',options = 2)
#         print(hull_cutoff_index_2)
#
#         hull_cutoff_index_3 = get_hull_cutoff(all_points,hull_all_index,cutoff_v,qhull_options = '',options = 3)
#         print(hull_cutoff_index_3)
#
#         hull_cutoff_index_4 = get_hull_cutoff(all_points,hull_all_index,cutoff_v,qhull_options = '',options = 4)
#         print(hull_cutoff_index_4)
#
#         hull_cutoff_1 = ConvexHull(all_points[hull_cutoff_index_1,:],qhull_options = qhull_options)
#         hull_cutoff_2 = ConvexHull(all_points[hull_cutoff_index_2,:],qhull_options = qhull_options)
#         hull_cutoff_3 = ConvexHull(all_points[hull_cutoff_index_3,:],qhull_options = qhull_options)
#         hull_cutoff_4 = ConvexHull(all_points[hull_cutoff_index_4,:],qhull_options = qhull_options)
#         hulls = [hull_all,hull_cutoff_1,hull_cutoff_2,hull_cutoff_3,hull_cutoff_4]
#         hull_plot(all_points, hulls,line_type = [],line_color = [])
