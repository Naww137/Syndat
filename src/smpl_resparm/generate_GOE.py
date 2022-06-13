#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 13:00:08 2022

@author: noahwalton
"""

import numpy as np

def generate_GOE(N):
    A = np.random.rand(N,N)/np.sqrt(2*N)
    X = A + np.transpose(A)
    return X

  