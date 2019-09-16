#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2019-09-12

"""Step1_iML1515_test.py
:description : script
:param : 
:returns: 
:rtype: 
"""
import cobra
import My_def
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import statistics
from cobra.flux_analysis import production_envelope
from cobra.flux_analysis import flux_variability_analysis
import pandas as pd



iML1515 = cobra.io.read_sbml_model('../ComplementaryData/iML1515.xml')




iML1515.reactions.get_by_id('EX_o2_e').bounds = (0.0,1000.0)
# % Constrain the phosphotransferase system
# model = changeRxnBounds(model, 'GLCabcpp', -1000, 'l');
# model = changeRxnBounds(model, 'GLCptspp', -1000, 'l');
# model = changeRxnBounds(model, 'GLCabcpp', 1000, 'u');
# model = changeRxnBounds(model, 'GLCptspp', 1000, 'u');
# model = changeRxnBounds(model, 'GLCt2pp', 0, 'b');
iML1515.reactions.get_by_id('GLCabcpp').bounds = (-1000.0,1000.0)
iML1515.reactions.get_by_id('GLCptspp').bounds = (-1000.0,1000.0)
iML1515.reactions.get_by_id('GLCt2pp').bounds = (0.0,0.0)
# iML1515.reactions.get_by_id('EX_o2_e').bounds = (0.0,1000.0)
solution = iML1515.optimize()
print(solution)



yiled_rea_ids = ['EX_ac_e','EX_for_e','EX_succ_e','EX_lac__D_e','EX_etoh_e']
yiled_rea = [iML1515.reactions.get_by_id(i) for i in yiled_rea_ids]
# %%

b = pd.DataFrame()
for i in range(1,10):
    iML1515.reactions.get_by_id('BIOMASS_Ec_iML1515_core_75p37M').bounds = (solution.objective_value/i,solution.objective_value/i)
    a = flux_variability_analysis(iML1515,yiled_rea)
    a['biomass'] = solution.objective_value/i
    b = b.append(a)

c = b.loc[ 'EX_etoh_e' , : ].plot(kind='line', x='biomass', y='maximum')


# c.plot(
#      kind='line', x='biomass', y='maximum')

plt.show()

# %%

import matplotlib.pyplot as plt
import math
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import numpy as np


# %% ConvexHull base

#hull = ConvexHull(points)
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

points  = np.array([[1,1],[2.01,1.5],[2,2],[2,1],[1,2],[1.5,1.5]])

hull = ConvexHull(points,qhull_options=['Qt','QJ','Pp','Qw','Qx'])    #'QG0'mean expect point 0(index)
index = hull.vertices
import matplotlib.pyplot as plt
fig = plt.figure()
plt.plot(points[:,0], points[:,1], 'o')
for simplex in hull.simplices:
     plt.plot(points[simplex, 0], points[simplex, 1], 'k-')
    #print(points[simplex, 0], points[simplex, 1])
plt.show()


# %% ConvexHull 90%

persent = 0.9
from itertools import combinations




base_v = persent*hull.area


ndim = hull.ndim

result = pd.DataFrame({'ndim':[],'comb':[]})

comb = []
ndim = []

if ndim<2:
    pass
else:
    for npoint in range(3,len(hull.vertices)):
        temp = list(combinations(range(0,len(points)), npoint))
        for i in temp:
            point_temp = np.append([points[temp[0]]],[(points[temp[1]])],axis = 0)
            hull = ConvexHull(points,qhull_options=['Qt','QJ','Pp','Qw','Qx'])

if ndim > 2:
    pass








# %% 3D
points = np.random.rand(30, 3)   # 30 random points in 2-D

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(points[:,0], points[:,1], points[:,2], 'o')

for simplex in hull.simplices:
     ax.plot(points[simplex, 0], points[simplex, 1], points[simplex, 2],'k-')
plt.show()



