#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 9/19/19

"""tp.py
:description : script to test hull result of scipy, pyhull and matlab hull
:param : matlab result
:returns:   all the same !!
:rtype: 
"""


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D #<-- Note the capitalization!
import math
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import numpy as np
from itertools import combinations
import os
import scipy.io as sio
from scipy.spatial import ConvexHull
from pyhull.convex_hull import ConvexHull as ConvexHull2

os.chdir('../ComplementaryData/')
def get_hull_points_scipy(all_points, option):

    hull = ConvexHull(all_points,qhull_options = option)    #'QG0'mean expect point 0(index)
    hull_index = hull.vertices
    sub_points = all_points[hull_index,:]
    hull_index.sort()
    return hull_index,sub_points,hull.volume

def get_hull_points_pyhull(all_points, option):

    hull = ConvexHull2(all_points,option)    #'QG0'mean expect point 0(index)
    hull_index = set()
    for i in hull.vertices:
        hull_index = hull_index | set(i)

    sub_points = all_points[list(hull_index),:]
    return hull_index,sub_points,0

# SCV 精度不够！！！
# points_2d =np.genfromtxt("../../../MATLAB/aumic-master/MYA_test_2/points_2d.csv", delimiter=",")
# points_3d =np.genfromtxt("../../../MATLAB/aumic-master/MYA_test_2/points_3d.csv", delimiter=",")
# points_4d =np.genfromtxt("../../../MATLAB/aumic-master/MYA_test_2/points_4d.csv", delimiter=",")
# points_sq =np.genfromtxt("../../../MATLAB/aumic-master/MYA_test_2/points_sq.csv", delimiter=",")
# points_glc_33 =np.genfromtxt("../../../MATLAB/aumic-master/MYA_test_2/points_glc_33.csv", delimiter=",")

points_all = sio.loadmat('../../../MATLAB/aumic-master/MYA_test_2/points_all.mat')
points_2d = points_all['points_2d']
points_3d = points_all['points_3d']
points_4d = points_all['points_4d']
points_sq = points_all['points_sq']
points_glc_33 = points_all['points_glc_33']


#%%

points_hull_index_matlab = sio.loadmat('../../../MATLAB/aumic-master/MYA_test_2/points_hull_index.mat')
hull_index_2d_matlab = points_hull_index_matlab['points_index_2d'].T - 1
hull_index_3d_matlab = points_hull_index_matlab['points_index_3d'].T - 1
hull_index_4d_matlab = points_hull_index_matlab['points_index_4d'].T - 1
hull_index_sq_matlab = points_hull_index_matlab['points_index_sq'].T - 1
hull_index_glc_33_matlab = points_hull_index_matlab['points_index_glc_33'].T - 1

option='Qt QJ Pp Qw Qx'     #'Qt QJ Pp Qw Qx' QJ is important
# %%

hull_index_2d_sc,sub_points_2d_sc,hull_volume_2d = get_hull_points_scipy(points_2d, option)
hull_index_3d_sc,sub_points_3d_sc,hull_volume_3d = get_hull_points_scipy(points_3d, option)
hull_index_4d_sc,sub_points_4d_sc,hull_volume_4d = get_hull_points_scipy(points_4d, option)
hull_index_sq_sc,sub_points_sq_sc,hull_volume_sq = get_hull_points_scipy(points_sq, option)
hull_index_glc_33_sc,sub_points_glc_33_sc,hull_volume_glc_33 = get_hull_points_scipy(points_glc_33, option)

# %%

hull_index_2d_pyhull,sub_points_2d_pyhull,hull_volume_2d_py = get_hull_points_pyhull(points_2d, option)
hull_index_3d_pyhull,sub_points_3d_pyhull,hull_volume_3d_py = get_hull_points_pyhull(points_3d, option)
hull_index_4d_pyhull,sub_points_4d_pyhull,hull_volume_4d_py = get_hull_points_pyhull(points_4d, option)
hull_index_sq_pyhull,sub_points_sq_pyhull,hull_volume_sq_py = get_hull_points_pyhull(points_sq, option)
hull_index_glc_33_pyhull,sub_pointglc_33_pyhull,hull_volume_glc_33_pyhull = get_hull_points_pyhull(points_glc_33, option)

#%%

print('compare scipy with matlab')

print(set(hull_index_2d_sc) - set(hull_index_2d_matlab[0]) , set(hull_index_2d_matlab[0]) - set(hull_index_2d_sc))
print(set(hull_index_3d_sc) - set(hull_index_3d_matlab[0]) , set(hull_index_3d_matlab[0]) - set(hull_index_3d_sc))
print(set(hull_index_4d_sc) - set(hull_index_4d_matlab[0]) , set(hull_index_4d_matlab[0]) - set(hull_index_4d_sc))
print(set(hull_index_sq_sc) - set(hull_index_sq_matlab[0]) , set(hull_index_sq_matlab[0]) - set(hull_index_sq_sc))
print(set(hull_index_glc_33_sc) - set(hull_index_glc_33_matlab[0]) , set(hull_index_glc_33_matlab[0]) - set(hull_index_glc_33_sc))


print('compare scipy with pyhull')

print(set(hull_index_2d_sc) - hull_index_2d_pyhull , hull_index_2d_pyhull - set(hull_index_2d_sc))
print(set(hull_index_3d_sc) - hull_index_3d_pyhull , hull_index_3d_pyhull - set(hull_index_3d_sc))
print(set(hull_index_4d_sc) - hull_index_4d_pyhull , hull_index_4d_pyhull - set(hull_index_4d_sc))
print(set(hull_index_sq_sc) - hull_index_sq_pyhull , hull_index_sq_pyhull - set(hull_index_sq_sc))
print(set(hull_index_glc_33_sc) - hull_index_glc_33_pyhull , hull_index_glc_33_pyhull - set(hull_index_glc_33_sc))

