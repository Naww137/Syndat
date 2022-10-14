#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 10:34:17 2022

@author: noahwalton
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

# def pois_noise(vector):
#     noise = []
#     for counts in vector:
#         noise.append(np.random.default_rng().poisson(lam=counts))
#     return vector + noise
def pois_noise(vector):
    """
    Samples poissonian noise around a vector of values.
    Parameters
    ----------
    vector : array-like
        vector of expected values.
    Returns
    -------
    noisy_vector
        Noisy vector sampled as poissonian around each expected value
    """
    noisy_vector = np.random.default_rng().poisson(lam=vector)
    return noisy_vector

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
    
    
    ctr = np.array(cts/(bw*trig))
    partial = np.array(1/(bw*trig))
    d_nctr = [np.sqrt((partial[i]**2)*dc**2) for i,dc in enumerate(d_cts)]
    
    return ctr, d_nctr

    

def generate_raw_count_data(sample_df, open_df, add_noise, trigo,trigs, k,K, Bi, b0,B0, alpha):
    """
    Generates raw count data for sample-in given a theoretical tranmission. 
    
    This function performs the inverse of the reduction process, calculating raw count data
    from a theoretical transmission. This process requires the assumption of known,
    true underlying reduction parameters and open count data.

    Parameters
    ----------
    sample_df : pandas.DataFrame
        Sample in dataframe with a column for theoretical tranmission ['theo_trans'] and energy ['E'].
    open_df : pandas.DataFrame
        Open dataframe, columns ['E'], ['bw']
    trigo : int
        Number of times the LINAC is fired for the open counts, corresponding to the number of times each channel is openned for counts.
    trigs : int
        Number of times the LINAC is fired for the sample-in counts, corresponding to the number of times each channel is openned for counts.
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
    sample_df : pandas.DataFrame
        Dataframe containing data for sample in.
    open_df : pandas.DataFrame
        Dataframe containing data for sample out.
    """
    # calculate open count rates
    Cr, dCr = cts_to_ctr(open_df.c, open_df.dc, open_df.bw, trigo) # cts_o/(bw*trig)
    # open_df['cps'] = Cr; open_df['dcps'] = dCr
    
    # calculate sample in count rate from theoretical transmission, bkg, m,k, and open count rate
    [m1,m2,m3,m4] = alpha
    cr = (sample_df.theo_trans*(m3*Cr - m4*K*Bi - B0) + m2*k*Bi + b0)/m1
    
    # calculate sample in counts from count rate
    c = cr*open_df.bw*trigs

    # correct for dead time
    # c = c*


    sample_df['theo_cts'] = c
    
    dc = np.sqrt(c)
    if add_noise:
        c = pois_noise(c)
        assert(c.all() >= 0)
        dc = np.sqrt(c)
        
    sample_df['c'] = c
    sample_df['dc'] = dc

    return sample_df





def get_covT(tof, c,C, dc,dC, a,b, k,K, Bi, b0,B0, alpha, sys_unc, ab_cov):
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
    ab_cov  : float
        Covariance between a & b background function parameters.
    
    Notes
    -----
    Background function must be of form Bi = a*exp(-b). Explicitly coded derivatives for Jacobian.
    
    Returns
    -------
    Output covaraiance matrix for transmission.
    """
    
    # cast all into numpy arrays
    c = np.array(c); dc = np.array(dc)
    C = np.array(C); dC = np.array(dC)
    tof = np.array(tof)
    Bi = np.array(Bi)
    sys_unc = np.array(sys_unc)

    # numerator and denominator
    D = alpha[2]*C - alpha[3]*K*Bi - B0
    N = alpha[0]*c - alpha[1]*k*Bi - b0
    
    # construct statistical covariance and jacobian
    dTi_dci = alpha[0]/D
    dTi_dCi = N*alpha[2]/D**2
    diag = (dc**2)*(dTi_dci**2) + (dC**2)*(dTi_dCi**2)
    CovT_stat = np.diag(diag)

    # construct systematic covariance and jacobian
    Cov_sys = np.diag(sys_unc**2)
    # print("WARNING: Need to update getCov function to take a/b covariances, currently it says cov = var*var")
    Cov_sys[0,1] = ab_cov #1.42659922e-1 # 1.42405866e-01 # sys_unc[0]*sys_unc[1]  
    Cov_sys[1,0] = ab_cov #1.42659922e-1 # 1.42405866e-01 #sys_unc[1]*sys_unc[0]        
    
    # systematic derivatives
    dTi_da = -1*(k*D+K*N)/(a*D**2)  #   -(k*alpha[1]*D+K*alpha[3]*N)*np.exp(-b*tof) / (D**2)
    dTi_db = (k*D+K*N)*Bi*tof/D**2  #   (k*alpha[1]*D)*Bi*tof / (D**2)
    dTi_dk = -alpha[1]*Bi/D       
    dTi_dK = N*alpha[3]*Bi/D**2
    dTi_db0 = -1/D
    dTi_dB0 = N/D**2
    dTi_dalpha = [ c/D, -k*Bi/D, -C*N/D**2, K*Bi*N/D**2 ]
    
    Jac_sys = np.array([dTi_da, dTi_db, dTi_dk, dTi_dK, dTi_db0, dTi_dB0, dTi_dalpha[0], dTi_dalpha[1],dTi_dalpha[2],dTi_dalpha[3]])
    
                
    # calculate covariance of output
    CovT_sys = Jac_sys.T @ Cov_sys @ Jac_sys
    
    CovT = CovT_stat + CovT_sys
    
    return CovT, CovT_stat, CovT_sys



def transmission(cr,Cr, Bi, k,K, b0,B0, alpha):
    [m1,m2,m3,m4] = alpha
    return (m1*cr - m2*k*Bi - b0) / (m3*Cr - m4*K*Bi - B0) 

    
def reduce_raw_count_data(tof, c,C, dc,dC, bw, trigo,trigs, a,b, k,K, Bi, b0,B0, alpha, sys_unc, ab_cov):
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
    ab_cov  : float
        Covariance between a & b background function parameters.
    
    Notes
    -----
    Background function must be of form Bi = a*exp(-b)
    
    Returns
    -------
    Output covaraiance matrix for transmission.
    """
    # calculate count rate and propagate uncertainty
    Cr, dCr = cts_to_ctr(C, dC, bw, trigo) 
    cr, dcr = cts_to_ctr(c, dc, bw, trigs)
    rates = Cr,dCr, cr,dcr
    
    # calculate transmission
    Tn = transmission(cr,Cr, Bi, k,K, b0,B0, alpha)
    # propagate uncertainty to transmission
    CovT, CovT_stat, CovT_sys = get_covT(tof, cr,Cr, dcr,dCr, a,b, k,K, Bi, b0,B0, alpha, sys_unc, ab_cov)
    dT = np.sqrt(np.diagonal(CovT))
    
    return Tn, dT, CovT, CovT_stat, CovT_sys, rates

