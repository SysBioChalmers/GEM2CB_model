#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 9/19/19

"""tp.py
:description : script
:param : 
:returns: 
:rtype: 
"""


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D #<-- Note the capitalization!
import math
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import numpy as np
from itertools import combinations
import os


os.chdir('../ComplementaryData/')
def get_hull_points_scipy(all_points, option):

    hull = ConvexHull(all_points,option)    #'QG0'mean expect point 0(index)
    hull_index = hull.vertices
    sub_points = all_points[hull_index,:]
    hull_index.sort()
    return hull_index,sub_points,hull.volume


def get_hull_points_pyhull(all_points, option):

    hull = ConvexHull(all_points,option)    #'QG0'mean expect point 0(index)
    hull_index = set()
    for i in hull.vertices:
        hull_index = hull_index | set(i)

    sub_points = all_points[list(hull_index),:]
    return hull_index,sub_points,0




# import scipy.io as sio
# points_all = sio.loadmat('../../../MATLAB/aumic-master/MYA_test_2/A_2d_rand.mat')
# oct_a = A['A']

points_2d =np.genfromtxt("../../../MATLAB/aumic-master/MYA_test_2/points_2d.csv", delimiter=",")
points_3d =np.genfromtxt("../../../MATLAB/aumic-master/MYA_test_2/points_3d.csv", delimiter=",")
points_4d =np.genfromtxt("../../../MATLAB/aumic-master/MYA_test_2/points_4d.csv", delimiter=",")
points_sq =np.genfromtxt("../../../MATLAB/aumic-master/MYA_test_2/points_sq.csv", delimiter=",")

option=''#'Qt QJ Pp Qw Qx'
# hull = ConvexHull(points_2d,qhull_options)    #'QG0'mean expect point 0(index)
# hull_all_index = hull.vertices
# hull_all_points = points_2d[hull_all_index,:]
# %%

from scipy.spatial import ConvexHull

hull_index_2d,sub_points_2d,hull_volume_2d = get_hull_points_scipy(points_2d, option)
hull_index_3d,sub_points_3d,hull_volume_3d = get_hull_points_scipy(points_3d, option)
hull_index_4d,sub_points_4d,hull_volume_4d = get_hull_points_scipy(points_4d, option)
hull_index_sq,sub_points_sq,hull_volume_sq = get_hull_points_scipy(points_sq, option)


# %%
from pyhull.convex_hull import ConvexHull
hull_index_2d_py,sub_points_2d_py,hull_volume_2d_py = get_hull_points_pyhull(points_2d, option)
hull_index_3d_py,sub_points_3d_py,hull_volume_3d_py = get_hull_points_pyhull(points_3d, option)
hull_index_4d_py,sub_points_4d_py,hull_volume_4d_py = get_hull_points_pyhull(points_4d, option)
hull_index_sq_py,sub_points_sq_py,hull_volume_sq_py = get_hull_points_pyhull(points_sq, option)


#%%

print(set(hull_index_2d) == hull_index_2d_py )

print(set(hull_index_3d) == hull_index_3d_py )

print(set(hull_index_4d) == hull_index_4d_py )

print(set(hull_index_sq) == hull_index_sq_py )



print(sub_points_2d - sub_points_2d_py )
print(sub_points_3d - sub_points_3d_py )
print(sub_points_4d - sub_points_4d_py )
print(sub_points_sq - sub_points_sq_py )

hull = ConvexHull(sub_points_2d_py,option)

#%%
import scipy.io as sio
sub_points_mat = sio.loadmat('../../../MATLAB/aumic-master/MYA_test_2/sub_points.mat')

hull_index_2d_mat = sub_points_mat['points_index_2d'].T-1
hull_index_3d_mat = sub_points_mat['points_index_3d'].T-1
hull_index_4d_mat = sub_points_mat['points_index_4d'].T-1
hull_index_sq_mat = sub_points_mat['points_index_sq'].T-1


print(set(hull_index_2d_mat[0]) == hull_index_2d_py )

print(set(hull_index_3d_mat[0]) == hull_index_3d_py )

print(set(hull_index_4d_mat[0]) == hull_index_4d_py )

print(set(hull_index_sq_mat[0]) == hull_index_sq_py )