#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 15:16:30 2022

@author: noahwalton
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import syndat
import scipy.stats as stat



#%%

plt.rcParams['figure.dpi'] = 500

sammy_directory = os.getcwd()

# =============================================================================
# # jesse's work
# trans = pd.read_csv(os.path.join(sammy_directory,'walton-trans.twenty'), sep='\s+', header=None)
# plt.errorbar(trans[0], trans[1], yerr=trans[2], ecolor='k', color='k', capsize=4, elinewidth=1.5, ms=2, fmt=".", zorder=0)
# #plt.close()
# 
# #theoretical from sammy
# dat = pd.read_csv(os.path.join(sammy_directory,'r-mat-test.dat'), sep='\s+', header=None)
# =============================================================================

sammy_lst = pd.read_csv(os.path.join(sammy_directory,'synthetic_data/SAMMY_jb.LST'), sep='\s+', names=['E','exp_xs','exp_xs_unc','theo_xs','theo_xs_bayes','exp_trans','exp_trans_unc','theo_trans', 'theo_trans_bayes'])
energy = sammy_lst['E']

plt.plot(energy, sammy_lst['theo_trans'], zorder=2, color='royalblue')
plt.xscale('log')
plt.show(); plt.close()


# next steps
# take this transmission data, add statistical noise and consider backgroud correction/correlation

# =============================================================================
# noise = np.random.default_rng().normal(loc=0.0, scale=np.mean(xs_theoretical)*1e-4, size=len(xs_theoretical))
# xs_experimental = xs_theoretical + noise
# 
# 
# plt.rcParams['figure.dpi'] = 500
# 
# plt.plot(energy,xs_theoretical, lw=0.5, label='$\sigma_{exp}$')
# plt.scatter(energy,xs_experimental, s=0.1, c='r', label='$\sigma_{exp}$')
# 
# plt.legend()
# plt.xlabel('Energy'); plt.ylabel('$\sigma$')
# plt.yscale('log'); plt.xscale('log')
# plt.show(); plt.close()
# =============================================================================
#%% saved plotting functions

def plot1(energy,theo,exp,label1,label2):
    
    plt.plot(energy, theo, label=label1, zorder=2)
    plt.scatter(energy,exp, label=label2, s=1, c='k', zorder=1)
    
    plt.legend()
    #plt.yscale('log'); 
    plt.xscale('log')
    plt.show();plt.close()
    

def plot2(energy,theo,exp,exp_unc, title):
    
    fig, (ax1,ax2) = plt.subplots(2,1, sharex=True, figsize=(6,4),gridspec_kw={'height_ratios': [2, 1]}) # , figsize=(12,5)
    plt.rcParams['figure.dpi'] = 500
    
    ax1.plot(energy,theo, lw=0.5, color='b', label='$T_{theo}$', zorder=2)
    #ax1.scatter(energy,exp, s=0.1, c='r', label='$T_{exp}$')
    ax1.errorbar(energy, exp, yerr=exp_unc, color='k',ecolor='k',elinewidth=1,capsize=2, fmt='.', ms=3, label='$T_{exp}$', zorder=0)
    
    ax1.legend()
    ax1.set_ylabel('T') #('$\sigma$')
    #ax1.set_yscale('log'); 
    ax1.set_xscale('log')
    ax1.set_ylim([0,max(exp)+0.1])
    
    rel_se = (exp-theo)/theo
    ax2.scatter(energy, rel_se, s=2)
    ax2.set_ylim([-.5,.5])
    ax2.set_xlabel('ToF (s)'); ax2.set_ylabel('L1 Norm (relative)')
    
    plt.suptitle(title)
    plt.tight_layout()
    plt.show(); plt.close()


#%%

def gaus_noise(vector, std_vec):
    # scale (std) = sqrt(mean) resembles almost identically the poisson distribution with a number of counts>20
    noise = np.random.default_rng().normal(loc=0.0, scale=std_vec, size=len(vector)) #np.mean(vector)*level
    return vector + noise

def pois_noise(vector):
    
    noise = []
    for counts in vector:
        noise.append(np.random.default_rng().poisson(lam=counts))
        
    return vector + noise
     
    

#path_to_lst = '/Users/noahwalton/Library/Mobile Documents/com~apple~CloudDocs/Research Projects/Resonance Fitting/summer_2022/SAMMY_jb.LST'
#sammy_lst = pd.read_csv(path_to_lst, sep='\s+', names=['E','exp_xs','exp_xs_unc','theo_xs','theo_xs_bayes','exp_trans','exp_trans_unc','theo_trans','theo_trans_bayes'])           

energy = np.array(sammy_lst['E'])
xs_theoretical = np.array(sammy_lst['theo_xs'])



# =============================================================================
#  simple uncertainty propagation from counts to transmission (no background)
# =============================================================================

# =============================================================================
# experimentally unique values
# =============================================================================
n = .12
trig = 10# number of linac pulses
bw = 1e-3 # bin width
flux_mag = 1e5 # what is a reasonable flux magnitude??
detector_efficiency = 1
tof_dist = 100 # m

# assuming 1% uncertainty of these values
m1 = 1; dm1 = m1*0.01
m2 = 1; dm2 = m2*0.01
m3 = 1; dm3 = m3*0.01
m4 = 1; dm4 = m4*0.01
alpha = [m1,m2,m3,m4]; d_alpha = [dm1,dm2,dm3,dm4]

a = 1; da = a*0.05
b = 1; db = b*0.05

# sample in
k = 0.1; dk = k*0.05
b0 = 1; db0 = b0*0.05
# sample out
K = 0.1; dK = K*0.05
B0 = 1; dB0 = B0*0.05

tof = syndat.exp_effects.e_to_t(energy,tof_dist,True)

# background function
def bkg_func(ti,a,b):
    return a*ti**-b
Bi = bkg_func(tof,a,b)






def generate_count_data(energy, xs_theo, flux_mag, bw, trig, n, Bi, k,K, b0,B0, alpha):
    
    [m1,m2,m3,m4] = alpha
    
    # =============================================================================
    # open counts, sample out- do we want this to be noisey or clean here?
    # =============================================================================
    cts_o = stat.norm.pdf(energy, loc=50, scale=100)*flux_mag # gaussian in energy, std=range of energy
    d_cts_o = np.sqrt(cts_o)
    #ncts_o = gaus_noise(cts_o,d_cts_o) # noisey open counts
    
    # Plot open flux spectra
    # =============================================================================
    # fig, (ax1,ax2) = plt.subplots(2,1, gridspec_kw={'height_ratios': [1, 1]})
    # ax1.plot(energy, C); ax1.set_xlabel('energy (eV)'); ax1.set_ylabel('C(E)')
    # ax2.plot(tof,C); ax2.set_xlabel('tof (s)'); ax1.set_ylabel('C(t)')
    # =============================================================================
    
    ctr_o = cts_o/(bw*trig)
    #nctr_o = ncts_o/(bw*trig) # noisy, open count rate
    
    # =============================================================================
    # find the number of counts given theoretical transmission
    # =============================================================================
    
    trans_theo = np.exp(-n*xs_theo) # sammy can also just output this transmission value
    
    ctr_i = (trans_theo*(m3*ctr_o - m4*K*Bi - B0) + m2*k*Bi + b0)/m1
    
    cts_i = ctr_i*bw*trig 
    d_cts_i = np.sqrt(cts_i)
    ncts_i = gaus_noise(cts_i,d_cts_i)

    return ncts_i






# =============================================================================
# now reduce the raw count data - can use the same, or different functions for C_open and background
# =============================================================================



# get noisey open count rate (sample out)
cts_o = stat.norm.pdf(energy, loc=50, scale=100)*flux_mag # gaussian in energy, std=range of energy
d_cts_o = np.sqrt(cts_o)
ncts_o = gaus_noise(cts_o,d_cts_o) # noisey open counts


ncrs_o = ncts_o/(bw*trig) # noisy open count rate

# generate noisey count rate for sample in
ncts_i = generate_count_data(energy,xs_theoretical, flux_mag,bw,trig, n, Bi, k,K, b0,B0, alpha)
d_ncts_i = np.sqrt(ncts_i)




def transmission(cr,Cr, Bi, k,K, b0,B0, alpha):
    [m1,m2,m3,m4] = alpha
    return (m1*c - m2*k*Bi - b0) / (m3*C - m4*K*Bi - B0) 

crn = cn/(bw*trig)

T = transmission(crn,Crn, Bi, k,K, b0,B0, alpha)




# =============================================================================
# sys_unc = np.append([da,db,dk,dK,db0,dB0], d_alpha)
# 
# CovT = syndat.exp_effects.get_covT(tof, Bi, dc,dC, sys_unc, a,b, k,K, c,C, b0,B0, alpha)
# 
# dT = np.diagonal(CovT)
# 
# plot2(energy,T,T_noise,dT, 'Fully Correlated Uncertainty')
# =============================================================================




#%%
# flat flux for sample out
C = np.ones((len(tof)))*flux_mag
dC = np.sqrt(C) # poisson counting statistics
C_noise = gaus_noise(C,dC)

c = C*np.exp(-xs_theoretical*n) * detector_efficiency
dc = np.sqrt(c) #poisson counting statistics
c_noise = gaus_noise(c,dc)

# need a proper function to convert to count rate and correct for deadtime
# the function for transmission below assumes c/C are count rates
# are these calculations included in the uncertainty propagation? technically the poisson unc is on the raw count #
# Cr = C/bin_width
# cr = c/bin_width


# calculate noisey transmission points and background 
def bkg_func(ti,a,b):
    return a*ti**-b

def transmission(alpha, c,C, k,K,  Bi, b0,B0):
    return (alpha[0]*c - alpha[1]*k*Bi - b0) / (alpha[2]*C - alpha[3]*K*Bi - B0) 
 
Bi = bkg_func(tof,a,b)

T = transmission(alpha, c,C, k,K, Bi, b0,B0)
T_noise = transmission(alpha, c_noise,C_noise, k,K, Bi, b0,B0)


# =============================================================================
# uncertainty propagation from poisson uncertainty on count rates
# =============================================================================





sys_unc = np.append([da,db,dk,dK,db0,dB0], dalpha)

CovT = syndat.exp_effects.get_covT(tof, Bi, dc,dC, sys_unc, a,b, k,K, c,C, b0,B0, alpha)

dT = np.diagonal(CovT)

plot2(energy,T,T_noise,dT, 'Fully Correlated Uncertainty')


inv_dT = np.linalg.inv(np.sqrt(np.diag(dT)))
CorrT = inv_dT @ CovT @ inv_dT
plt.matshow(CorrT)
cb = plt.colorbar()
#cb.ax.tick_params(labelsize=10)
plt.title('Correlation of Transmission', fontsize=10);
plt.show(); plt.close()


#%%
# =============================================================================
# Covariance matrix from ALEX manual
#   - derivation is done for channel to channel, when in actuality each data point I am looking at is a sum over multiple channels
#   - interesting how dc and dC don't show up, is this correct? see ALEX manual
#       - this is because dc = sqrt(c) and dc only shows up as dc^2 -> c
#   - in the ALEX manual, why is there not a subscript i bor the background spectra?
# =============================================================================

# =============================================================================
# def var_sig_i(sigi, c, C, d, D, n,dn, m,dm, M,dM, g,dg, G,dG):
#     var = 1/n * ( (sigi*dn**2) + (dM/M)**2 + (dm/m)**2 \
#                  + (dg/(d*c-g))**2 + (dG/(D*C-G))**2 \
#                  + d**2*c/(d*c-g)**2 + D**2*C/(D*C-G)**2 )
#     return var
#     
# def cov_sig_ij(sigi, sigj, di,ci,Di,Ci, dj,cj,Dj,Cj, n,dn, M,dM, m,dm, g,dg, G,dG):
#     cov_ij = 1/n**2 * ( sigi*sigj*dn**2 + (dM/M)**2 + (dm/m)**2 \
#                        + dg**2/((di*ci-g)*(dj*cj-g)) \
#                        + dG**2/((Di*Ci-G)*(Dj*Cj-G)) )
#     return cov_ij
#         
# energy = np.linspace(1,100,10)
# dn=0
# m = 1; dm = 0
# M = 1; dM = 0 
# g = 1; dg = .4
# G = 1; dG = .4
# 
# d = np.ones([len(energy)])
# D = np.ones([len(energy)])
# 
# n = .12
# time = 10
# flux_mag = 1e4
# detector_efficiency = 1
# 
# 
# C = np.ones((len(energy)))*flux_mag
# seC = np.sqrt(C) # poisson counting statistics
# C_noise = gaus_noise(C,seC)
# 
# c = C*np.exp(-xs_theoretical*n) * detector_efficiency
# sec = np.sqrt(c) #poisson counting statistics
# c_noise = gaus_noise(c,sec)
# 
# 
# 
# covmat = np.zeros([len(energy),len(energy)])
# for i in range(len(energy)):
#     for j in range(len(energy)):
#         
#         if i == j:
#             covmat[i,j] = var_sig_i(xs_theoretical[i],c[i],C[i],d[i],D[i], \
#                                     n,dn, m,dm, M,dM, g,dg, G,dG)
#         else:
#             covmat[i,j] = cov_sig_ij(xs_theoretical[i], xs_theoretical[j], \
#                                      d[i],c[i],D[i],C[i], \
#                                      d[j],c[j],D[j],C[j], \
#                                      n,dn, M,dM, m,dm, g,dg, G,dG)
#                 
#                 
# T = M*(d*c-g)/m*(D*C-G)
# T_noise = M*(d*c_noise-g)/m*(D*C_noise-G)
# 
# seT = np.sqrt((c/C**2) + (c**2/C**3))
# 
# Texp = gaus_noise(T,seT)
# 
# plot2(energy, T, T_noise)
# =============================================================================

# =============================================================================
# to sample a covariance matrix, but I don't think we want to do this
# covL = np.tril(covmat)
# uncor_unitvar = np.random.default_rng().normal(loc=0.0, scale=1, size=len(energy))
# 
# test = np.dot(uncor_unitvar, uncor_unitvar.T)
# =============================================================================









#%%

# open flux spectra measurement
# do we need an appropriate function for the incident neutron flux spectra
# =============================================================================
# flux0 = np.ones((len(energy)))*flux_mag
# se_0 = np.sqrt(flux0)
# 
# # detector measurement with sample in
# fluxdet = flux0*np.exp(-xs_theoretical*n)
# trans = fluxdet # APPLY SENSITYIVITY Of the detector - cps/flux - assumed to be 1
# se_trans = np.sqrt(trans)
# 
# 
# #se_net = np.sqrt((trans/time) + (flux0/time))
# trans_exp = gaus_noise(trans, se_trans)
# 
# #plot1(energy, flux0, trans, 'no sample', 'sample')
# 
# plot1(energy,trans,trans_exp, 'trans', 'exp')
# 
# plot2(energy, trans, trans_exp)
# =============================================================================


#%%
# =============================================================================
# 
# xs_experimental = gaus_noise(xs_theoretical,1e-2)
# 
# fig, (ax1,ax2) = plt.subplots(2,1, sharex=True, figsize=(5,5)) # , figsize=(12,5)
# plt.rcParams['figure.dpi'] = 500
# 
# 
# ax1.plot(energy,xs_theoretical, lw=0.5, label='$\sigma_{theo}$')
# ax1.scatter(energy,xs_experimental, s=0.1, c='r', label='$\sigma_{exp}$')
# 
# ax1.legend()
# ax1.set_ylabel('$\sigma$')
# ax1.set_yscale('log'); ax1.set_xscale('log')
# 
# 
# 
# ax2.scatter(energy, (xs_experimental-xs_theoretical), s=2)
# ax2.set_xlabel('Energy'); ax2.set_ylabel('L1 Norm')
# 
# plt.tight_layout()
# plt.show(); plt.close()
# 
# 
# 
# # if you provide a more course energy grid, you can resolve the structure of the cov/corr better
# #energy = np.linspace(0.1,100,50)
# 
# 
# # very simple background function from sammy
# def bkg_func(E, a,b,c,d,f):
#     return a + b/np.sqrt(E) + c*np.sqrt(E) + d*np.exp(-f/np.sqrt(E))
# 
# 
# #plt.plot(energy, bkg_func(energy, 0,10,0.2,10,1))
# 
# background = bkg_func(energy, 0,10,0.2,10,1)
# plt.plot(energy, background)
# plt.title('Background Function')
# plt.show(); plt.close()
# 
# def repeat_function(vector):
#     return np.repeat([vector,vector],len(vector)*2,axis=0)
# 
# 
# # if you make this perfectly correlated by just repeating the background function and taking corr, 
# # you get a perfectly 1:1 correlation between all values - is this what we want?
# # or do we want to sample noise around the background function? 
# #bkg_cov = np.corrcoef(repeat_function(background))
# 
# random_bkg = [gaus_noise(bkg_sample, 0.01) for bkg_sample in repeat_function(background)]
# bkg_cov = np.cov(random_bkg)
# 
# plt.matshow(bkg_cov)
# cb = plt.colorbar()
# cb.ax.tick_params(labelsize=10)
# plt.title('Cov of bkg function', fontsize=10);
# 
# plt.show(); plt.close()
# 
# =============================================================================



#%%











