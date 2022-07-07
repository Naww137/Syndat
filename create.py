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
    
    ax2.scatter(energy, (exp-theo)/theo, s=2)
    ax2.set_ylim([-1,1])
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



n = .12
time = 10
flux_mag = 1e2
detector_efficiency = 1

C = np.ones((len(energy)))*flux_mag
seC = np.sqrt(C) # poisson counting statistics

c = C*np.exp(-xs_theoretical*n) * detector_efficiency
sec = np.sqrt(c) #poisson counting statistics

T = c/C
seT = np.sqrt((c/C**2) + (c**2/C**3))

Texp = gaus_noise(T,seT)

plot2(energy, T, Texp)


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











