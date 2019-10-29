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

from cobra.flux_analysis import flux_variability_analysis

import cobra.test
model = cobra.test.create_test_model("textbook")
flux_variability_analysis(model, model.reactions[:10])