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


def get_covT(tof, Bi, dc,dC, sys_unc, a,b, k,K, c,C, b0,B0, alpha):
    
    """
    Calculates the output covariance matrix of transmission from input uncertainties.

    This function uses the covariance sandwhich rule:
    .. math:: C_y = J^T*C_x*J
    to propagate input variance-covariance to transmission data

    Parameters
    ----------
    tof : array-like
        Array of time of flight values for each data point - corresponds to energy.
    Bi : float or int
        Time dependent background shape function.
    dc : float
        Uncertainty in the count rate for sample-in.
    dC : array-like
        Uncertainty in the count rate for sample-out.
    Davg : array-like
        Nested list of average level spacing for each spin group number. First 
        list is for negative parity (J-) second is for positive parity (J+).
    Gavg : array-like
        Nested list of average widths for each spin group number. First 
        list is for negative parity (J-) second is for positive parity (J+).
    Gavg_swave : float
        Average width used to sample agregate capture widths. **Unsure of this value.
    print_out : bool
        User option to print out quantum spin (J) mapping to console.
    save_csv : bool
        User option to save resonance ladders to csv.
    sammy_run_folder : str
        Folder in which the csv(s) containing resparm ladders will be saved.

    Notes
    -----
    Unsure of the average capture width for Gg sampling.
    
    Returns
    -------
    Jn_df : DataFrame
        Pandas DataFrame conatining a resonance parameter ladder for each 
        quantum spin group with negative parity (all J-).
    Jp_df : DataFrame
        Pandas DataFrame conatining a resonance parameter ladder for each 
        quantum spin group with positive parity (all J+).
    """

    # derivatives
    D = alpha[2]*C - alpha[3]*K*Bi - B0
    N = alpha[0]*c - alpha[1]*k*Bi - b0
    
    def dTi_dci(i):
        return alpha[0]/D[i]
    def dTi_dCi(i):
        return N[i]*alpha[2]/D[i]
    def dT_dsys(i):
        dTi_da = -(k*alpha[1]*D[i]+K*alpha[3]*N[i])*np.exp(-b*tof[i]) / (D[i]**2)
        dTi_db = (k*alpha[1]*D[i])*Bi[i]*tof[i] / (D[i]**2)
        dTi_dk = -alpha[1]*Bi[i]/D[i]**2
        dTi_dK = -N[i]*alpha[3]*Bi[i]/D[i]**2
        dTi_db0 = -1/D[i]
        dTi_dB0 = N[i]/D[i]**2
        dTi_dalpha = [ c[i]/D[i], -k*Bi[i]/D[i], -C[i]*N[i]/D[i]**2, K*Bi[i]*N[i]/D[i]**2 ]
        return np.append([dTi_da, dTi_db, dTi_dk, dTi_dK, dTi_db0, dTi_dB0], dTi_dalpha)
    
    # construct statistical covariance and jacobian
    Cov_stat = np.zeros([len(tof*2),len(tof*2)])
    #samplein = True; sampleout = False
    for i in range(len(tof)):
        for j in range(len(tof)):
            if i == j:
                Cov_stat[i,j] = dc[i] 
    for i in range(len(tof),len(tof*2)):
        for j in range(len(tof),len(tof*2)):
            if i == j:
                Cov_stat[i,j] = dC[i]
    
    Jac_stat = np.zeros([len(tof*2),len(tof*2)])
    for i in range(len(tof)):
        for j in range(len(tof)):
            if i == j:
                Jac_stat[i,j] = dTi_dci(i)
    for i in range(len(tof),len(tof*2)):
        for j in range(len(tof),len(tof*2)):
            if i == j:
                Jac_stat[i,j] = dTi_dCi(i)
    
    # construct systematic covariance and jacobian

    Cov_sys = np.zeros([len(sys_unc),len(sys_unc)])
    for i in range(len(sys_unc)):
        for j in range(len(sys_unc)):
            if i == j:
                Cov_sys[i,j] = sys_unc[i]
    Cov_sys[0,1] = sys_unc[0]*sys_unc[1]  
    Cov_sys[1,0] = sys_unc[1]*sys_unc[0]        
    
    Jac_sys = np.zeros([len(sys_unc),len(tof)])
    for i in range(len(sys_unc)):
        for j in range(len(tof)):
            Jac_sys[i,j] = dT_dsys(j)[i]
                
    # calculate covariance of output
    CovT_stat = Jac_stat.T @ Cov_stat @ Jac_stat
    CovT_sys = Jac_sys.T @ Cov_sys @ Jac_sys
    
    CovT = CovT_stat + CovT_sys
    
    return CovT




