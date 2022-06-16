#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 12:18:04 2022

@author: nwx
"""

import numpy as np
import matplotlib.pyplot as plt


# assuming boundary condition selected s.t. shift factor is eliminated for s wave but not others!


ac = 6.7e-15 # m or 6.7 femtometers
M = 63
m = 1
hbar = 6.582119569e-16 # eV-s
c = 2.99792458e8 # m/s
m_eV = 939.565420e6 # eV/c^2

# !!! energy in eV !!!


def k_wavenumber(E):
    k = 1/hbar * M/(m+M) * np.sqrt(2*m_eV*E) * 1/c
    return k

def PS_recursive(E, ac, orbital_angular_momentum):
    
    rho = k_wavenumber(E)*ac
    
    S_vector = np.zeros([orbital_angular_momentum+1,len(E)])
    P_vector = np.ones([orbital_angular_momentum+1,len(E)])
    P_vector[0,:] *= rho
    
    for l in range(1,orbital_angular_momentum+1):
        S = (rho**2*(l-S_vector[l-1]) / ((l-S_vector[l-1])**2 + P_vector[l-1]**2)) - l
        P = rho**2*P_vector[l-1] / ((l-S_vector[l-1])**2 + P_vector[l-1]**2)
        
        S_vector[l,:]=S; P_vector[l,:]=P
        
    return [S_vector, P_vector]


def PS_explicit(E,ac,orbital_angular_momentum):
    rho = k_wavenumber(E)*ac
    if orbital_angular_momentum == 0:
        P = rho
        S = np.zeros(len(E))
    if orbital_angular_momentum == 1:
        P = rho**3/(1+rho**2)
        S = -1/(1+rho**2)
    if orbital_angular_momentum == 2:
        P = rho**5/(9+3*rho**2+rho**4)
        S = -(18+3*rho**2)/(9+3*rho**2+rho**4)
    if orbital_angular_momentum == 3:
        P = rho**7/(255+45*rho**2+6*rho**4+rho**6)
        S = -(675+90*rho**2+6*rho**4)/(225+45*rho**2+6*rho**4+rho**6)
        
    return [S, P]
        
        
        


        
    



