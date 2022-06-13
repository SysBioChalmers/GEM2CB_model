#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 14:52:18 2018

@author: Stefan Vet (stefan.vet@vub.be)
"""

# Packages

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi


""" Initialization """

# Time
t_end = 50.    # Total integration time
t_start = 0.
t_step = 0.05

x_0 = np.array([1.31,0.11,0.28])    # RI, FP, BH
sx_0 = 30.
s_0 = np.array([47.4,1.13,4.19,1.78,sx_0])    # fructose, formate, acetate, butyrate, unknown compound, H2, CO2 (treated separately due to higher sampling rate)
q_0 = 10.**np.array([0,-5,-1.])    # lag_RI, lag_FP, lag_BH

y_init = np.append(x_0,np.append(s_0,q_0))


n_x = len(x_0)
n_s = len(s_0)
n_q = len(q_0)

dim = n_x + n_s         # Dimension of the system (without lag phase)
dim_tot = dim + n_q    # total different dimensions


# Maximum growth rate. RI, FP, BH
mu = np.ones(3)    

# Monod constants, only non-zero values for elements that correspond to nutrient being consumed by a species. Rows: RI, FP, BH. Columns: 0:fructose, 1:formate, 2:acetate, 3:butyrate, 4:unknown compound, 5:H2, 6:CO2
K = np.zeros((n_x,n_s))

# Weights for additive growth dependence on nutrients. row 1: RI, row 2: FP - column 1: fructose/fructose+unknownComp, column 2: acetate
w = np.ones(3)

# Consumption matrix: positive values correspond to consumption, negative values correspond to production of nutrients.
v = np.zeros((n_s,n_x))     # 1st-2nd-3rd column: RI, FP, BH 


# Maximum growth rate. RI, FP, BH
mu[0],mu[1],mu[2] =  2.125,2.397,1.823    #([0.976,1.43,2.])

# Monod constants, only non-zero values for elements that correspond to nutrient being consumed by a species. Rows: RI, FP, BH. Columns: 0:fructose, 1:formate, 2:acetate, 3:butyrate, 4:unknown compound, 5:H2, 6:CO2
K[0,0],K[0,2] = 205,  71. # RI
K[1,0],K[1,2],K[1,4] = 12,41,62    # FP
K[2,0],K[2,1] = 94,  413    # BH

# Weights for additive growth dependence on nutrients. row 1: RI, row 2: FP - column 1: fructose/fructose+unknownComp, column 2: acetate
w[0],w[1],w[2] = 0.85,0.237,1. #np.array([[1.,1.2],[0.1,1.]])

# Consumption matrix: positive values correspond to consumption, negative values correspond to production of nutrients.
v[0,:] = [0.364,   1.962,   0.389]     # Fructose
v[1,:] = [  0.012,  -1.684,   2.516]    # Formate
v[2,:] = [1.03 ,  11.443,  -0.617]        # Acetate
v[3,:] = [ -0.528,  -2.263,   0.  ]        # Butyrate
v[4,:] = [0.,2.485,0.]            # Unknown nutrient Sx




""" Definitions for the dynamics """ 


def lagPhase(q):
     """ Lag phases for each species"""
     return q/(1+q)

def Rate(s):
    """ Growth rates for separate nutrients. These terms are necessary for facultative consumption, e.g. of acetate (S2)"""
    rate = np.zeros(3)
    #RI
    rate_RI = np.array([1.,w[0]*s[2]/(K[0,2]+s[2])])*s[0]/(K[0,0]+s[0])
    rate[0] = sum(rate_RI)
    
    #FP
    rate_FP = np.array([1.,w[1]*s[2]/(K[1,2]+s[2])])*s[4]/(K[1,4]+s[4])*s[0]/(K[1,0]+s[0])
    rate[1] = sum(rate_FP)

    #BH
    rate_BH = np.array([s[0]/(K[2,0]+s[0]),w[2]*s[1]/(K[2,1]+s[1])])    #Additive growth rate: choice between fructose and/or CO2 (neglect influence of formate)
    rate[2] = sum(rate_BH)
    
    return [rate,rate_RI,rate_FP,rate_BH]


def growthRate(s,q):
    """ Growth rates for the species, depends on the nutrients s and the lag phases q. """

    rate,rate_RI,rate_FP,rate_BH = Rate(s)    # 3x2-matrix: rows: bacteria - columns: the terms that need to multiplied by the weights and then added. -> matrix product with w_matrix.
    lag = lagPhase(q)
    
    # General growth rates
    rate = lag*mu*rate
    rate_RI = lag[0]*mu[0]*rate_RI
    rate_FP = lag[1]*mu[1]*rate_FP
    rate_BH = lag[2]*mu[2]*rate_BH
        
    return [rate,rate_RI,rate_FP,rate_BH]


def dynamics(t,y):
    """ Definitions of the ODE's """

    x=y[0:n_x]
    s=y[n_x:dim]
    q=y[dim:]
           
    # Growth Rates and lag phases
    growth,growth_RI,growth_FP,growth_BH = growthRate(s,q)    # General growth rate
          
    # the model equations 
    dyn_x = growth*x    # RI, FP, BH at once        
    dyn_s = -v.dot(dyn_x)    # this equations takes off too much (growth rates composed of two parts), this will need to be added
    
    # For S0, S1 and S2: only distract seperate rate for nutrient for facultative consumption.
    dyn_s[0] = -v[0,0]*(growth_RI[0]*x[0] + growth_RI[1]*x[0])  -v[0,1]*(growth_FP[0]*x[1] + growth_FP[1]*x[1]) -v[0,2]*growth_BH[0]*x[2]
    dyn_s[2] =  -v[2,0]*growth_RI[1]*x[0]  -v[2,1]*growth_FP[1]*x[1] -v[2,2]*(growth_BH[0]*x[2]+growth_BH[1]*x[2]) 
    dyn_s[1] = -v[1,0]*(growth_RI[0]*x[0] + growth_RI[1]*x[0])  -v[1,1]*(growth_FP[0]*x[1] + growth_FP[1]*x[1]) -v[1,2]*growth_BH[1]*x[2] 
    
    dyn_q = mu*q*(q<10.**5)    # Lag phase: ignore this once q > 10.**5
    
    dyn = np.append(dyn_x,np.append(dyn_s,dyn_q))
                  
    return dyn


def simulate():
    """ Simulate the ODE's """

    # Solve ODE
    ode = spi.ode(dynamics)        #R = odeint(dynamics,x_init,t_grid)

    # BDF method suited to stiff systems of ODEs
    ode.set_integrator('vode',nsteps=500,method='bdf')
    ode.set_initial_value(y_init,t_start)

    ts = []
    ys = []

    while ode.successful() and ode.t < t_end:
        ode.integrate(ode.t + t_step)
        ts.append(ode.t)
        ys.append(ode.y)

    t_grid = np.vstack(ts)
    y_total = np.vstack(ys).T    #y = (x,s,q)
    
    return [t_grid,y_total]



def figure(t_grid,y_total,name):
    """ Make figure with subfigs for the time trace of species and nutrients """ 
    
    labels_x = ["RI","FP","BH"]
    color_x = ["red","blue","green"]
    labels_s = ["fructose","formate","acetate","butyrate","unknown"]
    color_s = np.array(['red','green','magenta', 'midnightblue','coral'])

    x_total = y_total[0:n_x,:]
    s_total = y_total[n_x:dim,:]

    plt.figure(1)
    ax = plt.subplot2grid((2,6), (0,0), colspan=5)    # Bacteria

    for i, color in enumerate(color_x,start=0):
        ax.plot(t_grid,x_total[i,:],'-',color=color,label=labels_x[i])
        
    ax.legend(loc='center right', bbox_to_anchor=(1.02, 0., .2, 1.),
              ncol=1, mode="expand", borderaxespad=0., fancybox=True)
    ax.grid()
    ax.set_ylabel(r"X ($10^8$ Counts/mL)",fontsize=12)
    ax.set_xlim(0,t_end)
    ax.set_ylim(0,150)

    ax2 = plt.subplot2grid((2,6), (1,0), colspan=5)    # Nutrients    
    
    for i,color in enumerate(color_s, start=0):
        ax2.plot(t_grid,s_total[i,:],'-',color=color,label=labels_s[i])
    ax2.legend(loc='center right', bbox_to_anchor=(1.02, 0., .3, 1.),
              ncol=1, mode="expand", borderaxespad=0., fancybox=True)
    ax2.grid()
    ax2.set_ylabel("S (mM)",fontsize=12)
    ax2.set_xlabel("Time (hours)",fontsize=12)
    ax2.set_xlim(0,t_end)
    ax2.set_ylim(0,100)

    plt.savefig(name)
    plt.show()
    plt.close()

    return


t_grid,y_total = simulate()
figure(t_grid,y_total,"timeTrace.png")
