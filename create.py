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
    

def plot2(energy,theo,exp):
    
    fig, (ax1,ax2) = plt.subplots(2,1, sharex=True, figsize=(5,5)) # , figsize=(12,5)
    plt.rcParams['figure.dpi'] = 500
    
    ax1.plot(energy,theo, lw=0.5, label='$T_{theo}$')
    ax1.scatter(energy,exp, s=0.1, c='r', label='$T_{exp}$')
    
    ax1.legend()
    ax1.set_ylabel('T') #('$\sigma$')
    #ax1.set_yscale('log'); 
    ax1.set_xscale('log')
    
    rel_se = (exp-theo)/theo
    ax2.scatter(energy, rel_se, s=2)
    ax2.set_ylim([-.1,.1])
    ax2.set_xlabel('Energy'); ax2.set_ylabel('L1 Norm (relative)')
    
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
n = .12
time = 10
flux_mag = 1e5
detector_efficiency = 1


C = np.ones((len(energy)))*flux_mag
seC = np.sqrt(C) # poisson counting statistics

c = C*np.exp(-xs_theoretical*n) * detector_efficiency
sec = np.sqrt(c) #poisson counting statistics

T = c/C
seT = np.sqrt((c/C**2) + (c**2/C**3))

Texp = gaus_noise(T,seT)

plot2(energy, T, Texp)


#%%
# =============================================================================
# rather than propagating uncertainty to T, I'm going to build a covariance matrix
# and sample from it s.t. the propagted uncertainties will be correlated with one another
# also, this will be done for the cross section sigmat because the cov derivation is greatly simplified by the ln()
#   - derivation is done for channel to channel, when in actuality each data point I am looking at is a sum over multiple channels
#   - interesting how dc and dC don't show up, is this correct? see ALEX manual
#       - this is because dc = sqrt(c) and dc only shows up as dc^2 -> c
#   - in the ALEX manual, why is there not a subscript i bor the background spectra?
# =============================================================================

def var_sig_i(sigi, c, C, d, D, n,dn, m,dm, M,dM, g,dg, G,dG):
    var = 1/n * ( (sigi*dn**2) + (dM/M)**2 + (dm/m)**2 \
                 + (dg/(d*c-g))**2 + (dG/(D*C-G))**2 \
                 + d**2*c/(d*c-g)**2 + D**2*C/(D*C-G)**2 )
    return var
    
def cov_sig_ij(sigi, sigj, di,ci,Di,Ci, dj,cj,Dj,Cj, n,dn, M,dM, m,dm, g,dg, G,dG):
    cov_ij = 1/n**2 * ( sigi*sigj*dn**2 + (dM/M)**2 + (dm/m)**2 \
                       + dg**2/((di*ci-g)*(dj*cj-g)) \
                       + dG**2/((Di*Ci-G)*(Dj*Cj-G)) )
    return cov_ij
        
energy = np.linspace(1,100,10)
dn=0
m = 1; dm = 0
M = 1; dM = 0 
g = 1; dg = .4
G = 1; dG = .4

d = np.ones([len(energy)])
D = np.ones([len(energy)])


covmat = np.zeros([len(energy),len(energy)])
for i in range(len(energy)):
    for j in range(len(energy)):
        
        if i == j:
            covmat[i,j] = var_sig_i(xs_theoretical[i],c[i],C[i],d[i],D[i], \
                                    n,dn, m,dm, M,dM, g,dg, G,dG)
        else:
            covmat[i,j] = cov_sig_ij(xs_theoretical[i], xs_theoretical[j], \
                                     d[i],c[i],D[i],C[i], \
                                     d[j],c[j],D[j],C[j], \
                                     n,dn, M,dM, m,dm, g,dg, G,dG)
                
                

covL = np.tril(covmat)
uncor_unitvar = np.random.default_rng().normal(loc=0.0, scale=1, size=len(energy))

test = np.dot(uncor_unitvar, uncor_unitvar.T)










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











