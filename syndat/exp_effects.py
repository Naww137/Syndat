#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 10:34:17 2022

@author: nwx
"""

import numpy as np
import scipy.stats as stat


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

def gaus_noise(vector, std_vec):
    """
    Samples gaussian noise around a vector of mean values.

    Parameters
    ----------
    vector : array-like
        Vector of mean values.
    std_vec : array-like
        Vector of standard deviations (standard errors).
    
    Returns
    -------
    Noisy vector sampled as guassian around each mean/std.
    """
    # scale (std) = sqrt(mean) resembles almost identically the poisson distribution with a number of counts>20
    noise = np.random.default_rng().normal(loc=0.0, scale=std_vec, size=len(vector)) #np.mean(vector)*level
    return vector + noise

def pois_noise(vector):
    noise = []
    for counts in vector:
        noise.append(np.random.default_rng().poisson(lam=counts))
    return vector + noise

def generate_open_counts(energy, flux_mag, mean, std):
    """
    Generate open (sample out) raw count data from a wide gaussian wrt energy.

    Parameters
    ----------
    energy : array-like
        Array of energy values for each data point - corresponds to tof.
    flux_mag : float
        Magnitude scaling factor applied to flux shape, how many counts!
    mean : float
        Average value for gaussian shape.
    std : float
        Standard deviation for gaussian shape.
    
    Returns
    -------
    Open counts (sample out).
    """
    cts_o = stat.norm.pdf(energy, loc=mean, scale=std)*flux_mag # gaussian in energy, std=range of energy
    return cts_o

def cts_to_ctr(cts, d_cts, bw, trig):
    """
    Converts counts to count rate and propagates uncertainty.

    Parameters
    ----------
    cts : array-like
        Array of count data corresponting to each tof bin.
    d_cts : array-like
        Array of uncertainty on each count data point corresponting to each tof bin.
    bw : array-like
        Array of tof bin widths.
    trig : float or int
        Number of linac pulses.
        
    Notes
    _____
    Uncertainty propagation with sandwich rule (JxCxJ.T) is over 1000x slower. 
    A more simple error propagtion is used because there is no covariance 
    between the statistical variables.
    
    Returns
    -------
    ctr : array-like
        Array of count rates corresponding to each tof bin.
    d_nctr : array-like
        Array of propagated uncertainty on each count rate point.
    """
    ctr = cts/(bw*trig)
    partial = 1/(bw*trig) 
# Matrix multiplication takes for ever here! no covariance, so take advantage of that!
# =============================================================================
#     Cin = np.diag(d_cts**2)
#     J = np.diag(np.ones(len(d_cts))*partial)
#     Cout = J.T @ Cin @ J
#     alt_d_nctr = np.sqrt(np.diag(Cout))
# =============================================================================
    d_nctr = [np.sqrt((partial[i]**2)*dc**2) for i,dc in enumerate(d_cts)]
# =============================================================================
#     if sum(d_nctr-alt_d_nctr) > 1e-10:
#         print('Warning: JxCxJ.T != nsqrt((d_dx*dx**2))')
# =============================================================================
    
    return ctr, d_nctr

    

def generate_raw_count_data(energy, T_theo, C, bw, trig, k,K, Bi, b0,B0, alpha):
    """
    Generates raw count data for sample-in given a theoretical tranmission. 
    
    This function performs the inverse of the reduction process, calculating raw count data
    from a theoretical transmission. This process requires the assumption of know,
    true underlying reduction parameters.

    Parameters
    ----------
    energy : array-like
        Vector of energy values for each data point - corresponds to tof.
    T_theo : array-like
        Vector of theoretical transmission values that the detector will see, (must be experimentally corrected).
    C : array-like
        Open count data - sample out counts.
    bw : float
        Width in time a given channel is open, bin width.
    trig : int
        Number of times the LINAC is fired, corresponding to the number of times each channel is openned for counts.
    k : float
        Background normalization for sample in.
    K : float
        Background normalization for sample out.
    Bi : array-like
        Background shape function stored as a vector.
    b0 : float
        Constant background for sample in.
    B0 : float
        Constant background for sample out.
    alpha : array-like
        Vector of monitor stability factors [m1,m2,m3,m4]
    
    Returns
    -------
    nc : array-like
        Noisy raw count data for sample in - noise sampled from Poisson wrt true underlying counts.
    dnc : array-like
        Uncertainty (standard error) associated with nc from Poisson.
    """
    #calculate open count rates
    Cr, dCr = cts_to_ctr(C, np.sqrt(C), bw, trig) # cts_o/(bw*trig)
    
    # calculate sample in count rate from theoretical transmission, bkg, m,k, and open count rate
    [m1,m2,m3,m4] = alpha
    cr = (T_theo*(m3*Cr - m4*K*Bi - B0) + m2*k*Bi + b0)/m1
    
    # calculate sample in counts, noise, and uncertainty
    c = cr*bw*trig 
    dc = np.sqrt(c)
    nc = gaus_noise(c,dc) # will create some negative counts, force to zero
    nc = np.where(nc<0, 0, nc) # replace negative counts with 0
    dnc = np.sqrt(nc)

    return nc, dnc, c, dc





def get_covT(tof, c,C, dc,dC, a,b, k,K, Bi, b0,B0, alpha, sys_unc):
    """
    Calculates the output covariance matrix of transmission from input uncertainties.

    This function uses the covariance sandwhich rule:
    .. math:: C_y = J^T*C_x*J
    to propagate input variance-covariance to transmission data

    Parameters
    ----------
    tof : array-like
        Array of time of flight values for each data point - corresponds to energy.
    c : float
        Count rate for sample in.
    C : float
        Count rate for sample out.
    dc : float
        Uncertainty in the count rate for sample-in.
    dC : array-like
        Uncertainty in the count rate for sample-out.
    a : float
        Shaping parameter for exponential background function.
    b : float
        Shaping parameter for exponential background function.
    k : float
        Background normalization for sample in.
    K : float
        Background normalization for sample out.
    Bi : array-like
        Background shape function stored as a vector.
    b0 : float
        Constant background for sample in.
    B0 : float
        Constant background for sample out.
    alpha : array-like
        Vector of monitor stability factors [m1,m2,m3,m4]
    sys_unc : array-like
        Vector of systematic uncertainties: [da,db,dk_i,dk_o,dB0_i,dB0_o,m1,m2,m3,m4].
    
    Notes
    -----
    Background function must be of form Bi = a*exp(-b). Explicitly coded derivatives for Jacobian.
    
    Returns
    -------
    Output covaraiance matrix for transmission.
    """

# =============================================================================
# make inputs the proper type
# =============================================================================
    
    dc = np.array(dc)
    dC = np.array(dC)

    # derivatives
    D = alpha[2]*C - alpha[3]*K*Bi - B0
    N = alpha[0]*c - alpha[1]*k*Bi - b0
    
# =============================================================================
#     def dTi_dci(i):
#         return alpha[0]/D[i]
#     def dTi_dCi(i):
#         return N[i]*alpha[2]/D[i]
# =============================================================================
    def dT_dsys(i):
        dTi_da = -(k*alpha[1]*D[i]+K*alpha[3]*N[i])*np.exp(-b*tof[i]) / (D[i]**2)
        dTi_db = (k*alpha[1]*D[i])*Bi[i]*tof[i] / (D[i]**2)
        dTi_dk = -alpha[1]*Bi[i]/D[i]**2
        dTi_dK = N[i]*alpha[3]*Bi[i]/D[i]**2
        dTi_db0 = -1/D[i]
        dTi_dB0 = N[i]/D[i]**2
        dTi_dalpha = [ c[i]/D[i], -k*Bi[i]/D[i], -C[i]*N[i]/D[i]**2, K*Bi[i]*N[i]/D[i]**2 ]
        return np.append([dTi_da, dTi_db, dTi_dk, dTi_dK, dTi_db0, dTi_dB0], dTi_dalpha)
    
    
    # construct statistical covariance and jacobian
    dc_dC = np.append(dc,dC)
    Cov_stat = np.diag(dc_dC**2)
# =============================================================================
#     Cov_stat = np.zeros([len(tof*2),len(tof*2)])
#     #samplein = True; sampleout = False
#     for i in range(len(tof)):
#         for j in range(len(tof)):
#             if i == j:
#                 Cov_stat[i,j] = dc[i] 
#     for i in range(len(tof),len(tof*2)):
#         for j in range(len(tof),len(tof*2)):
#             if i == j:
#                 Cov_stat[i,j] = dC[i]
# =============================================================================
    dTi_dci = alpha[0]/D; dTi_dci = np.diag(dTi_dci**2)
    dTi_dCi = N*alpha[2]/D**2; dTi_dCi = np.diag(dTi_dCi**2)
    Jac_stat = np.vstack((dTi_dci,dTi_dCi))
# =============================================================================
#     Jac_stat = np.zeros([len(tof)*2,len(tof)])
#     for i in range(len(tof)):
#         for j in range(len(tof)):
#             if i == j:
#                 Jac_stat[i,j] = dTi_dci(i)
#     for i in range(len(tof),len(tof*2)):
#         for j in range(len(tof)):
#             if i == j:
#                 Jac_stat[i,j] = dTi_dCi(i)
# =============================================================================
    
    # construct systematic covariance and jacobian
    Cov_sys = np.diag(sys_unc**2)
# =============================================================================
#     Cov_sys = np.zeros([len(sys_unc),len(sys_unc)])
#     for i in range(len(sys_unc)):
#         for j in range(len(sys_unc)):
#             if i == j:
#                 Cov_sys[i,j] = sys_unc[i]
# =============================================================================
    # print("WARNING: Need to update getCov function to take a/b covariances, currently it says cov = var*var")
    Cov_sys[0,1] = 1.42405866e-01 # sys_unc[0]*sys_unc[1]  
    Cov_sys[1,0] = 1.42405866e-01 #sys_unc[1]*sys_unc[0]        
    
    Jac_sys = np.zeros([len(sys_unc),len(tof)])
    for j in range(len(tof)):
        Jac_sys[:,j] = dT_dsys(j)
                
    # calculate covariance of output
    CovT_stat = Jac_stat.T @ Cov_stat @ Jac_stat
    CovT_sys = Jac_sys.T @ Cov_sys @ Jac_sys
    
    CovT = CovT_stat + CovT_sys
    
    return CovT



def transmission(cr,Cr, Bi, k,K, b0,B0, alpha):
    [m1,m2,m3,m4] = alpha
    return (m1*cr - m2*k*Bi - b0) / (m3*Cr - m4*K*Bi - B0) 

    
def reduce_raw_count_data(tof, c,C, dc,dC, bw, trig, a,b, k,K, Bi, b0,B0, alpha, sys_unc):
    """
    Reduces raw count data to transmission data with propagated uncertainty.

    This function uses the covariance sandwhich rule:
    .. math:: C_y = J^T*C_x*J
    to propagate input variance-covariance from both statistical uncertainties
    and systematic uncertainties to transmission data.

    Parameters
    ----------
    tof : array-like
        Array of time of flight values for each data point - corresponds to energy.
    c : float
        Count rate for sample in.
    C : float
        Count rate for sample out.
    dc : float
        Uncertainty in the count rate for sample-in.
    dC : array-like
        Uncertainty in the count rate for sample-out.
    bw : float
        Width in time a given channel is open, bin width.
    trig : int
        Number of times the LINAC is fired, corresponding to the number of times each channel is openned for counts.
    a : float
        Shaping parameter for exponential background function.
    b : float
        Shaping parameter for exponential background function.
    k : float
        Background normalization for sample in.
    K : float
        Background normalization for sample out.
    Bi : array-like
        Background shape function stored as a vector.
    b0 : float
        Constant background for sample in.
    B0 : float
        Constant background for sample out.
    alpha : array-like
        Vector of monitor stability factors [m1,m2,m3,m4]
    sys_unc : array-like
        Vector of systematic uncertainties: [da,db,dk_i,dk_o,dB0_i,dB0_o,m1,m2,m3,m4].
    
    Notes
    -----
    Background function must be of form Bi = a*exp(-b)
    
    Returns
    -------
    Output covaraiance matrix for transmission.
    """
    # calculate count rate and propagate uncertainty
    Cr, dCr = cts_to_ctr(C, dC, bw, trig) 
    cr, dcr = cts_to_ctr(c, dc, bw, trig)
    
    # calculate transmission
    Tn = transmission(cr,Cr, Bi, k,K, b0,B0, alpha)
    # propagate uncertainty to transmission
    CovT = get_covT(tof, cr,Cr, dcr,dCr, a,b, k,K, Bi, b0,B0, alpha, sys_unc)
    dT = np.sqrt(np.diagonal(CovT))
    
    return Tn, dT, CovT

