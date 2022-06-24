#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 10:34:17 2022

@author: nwx
"""

import numpy as np



def t_to_e(t,d,rel):
    if rel:
        mn = 939.56542052e6 # eV/c2
        c = 299792458 # m/s
        E = mn*(1/np.sqrt(1-(d/t/c)**2)-1)
    else:
        mn = 1.674927498e-27 #kg
        jev = 1.6022e-19 # J/eV
        E = 0.5*mn*(d/t)**2 /jev # eV
    return E


def e_to_t(E,d,rel):
    if rel:
        mn = 939.56542052e6 # eV/c2
        c = 299792458 # m/s
        t = d/c * 1/np.sqrt(1-1/(E/mn+1)**2)
    else:
        jev = 1.6022e-19 # J/eV
        mn = 1.674927498e-27 #kg
        t = d/np.sqrt(E*jev*2/mn)
    return t






