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



def gaus_noise(vector, level):
    noise = np.random.default_rng().normal(loc=0.0, scale=np.mean(vector)*level, size=len(vector))
    return vector + noise


path_to_lst = '/Users/noahwalton/Library/Mobile Documents/com~apple~CloudDocs/Research Projects/Resonance Fitting/summer_2022/SAMMY_jb.LST'



sammy_lst = pd.read_csv(path_to_lst, sep='\s+', names=['E','exp_xs','exp_xs_unc','theo_xs','theo_xs_bayes','exp_trans','exp_trans_unc','theo_trans','theo_trans_bayes'])           

energy = sammy_lst['E']
xs_theoretical = sammy_lst['theo_xs']

noise = gaus_noise(xs_theoretical,1e-2)
xs_experimental = xs_theoretical + noise




#%%

fig, (ax1,ax2) = plt.subplots(2,1, sharex=True, figsize=(5,5)) # , figsize=(12,5)
plt.rcParams['figure.dpi'] = 500


ax1.plot(energy,xs_theoretical, lw=0.5, label='$\sigma_{theo}$')
ax1.scatter(energy,xs_experimental, s=0.1, c='r', label='$\sigma_{exp}$')

ax1.legend()
ax1.set_ylabel('$\sigma$')
ax1.set_yscale('log'); ax1.set_xscale('log')



ax2.scatter(energy, (xs_experimental-xs_theoretical), s=2)
ax2.set_xlabel('Energy'); ax2.set_ylabel('L1 Norm')

plt.show(); plt.close()


#%%

# if you provide a more course energy grid, you can resolve the structure of the cov/corr better
#energy = np.linspace(0.1,100,50)


# very simple background function from sammy
def bkg_func(E, a,b,c,d,f):
    return a + b/np.sqrt(E) + c*np.sqrt(E) + d*np.exp(-f/np.sqrt(E))

# =============================================================================
# def pois_noise(vector, lam):
#     noise = random.Generator.poisson(lam=1.0, size=None)
# =============================================================================

#plt.plot(energy, bkg_func(energy, 0,10,0.2,10,1))

background = bkg_func(energy, 0,10,0.2,10,1)
plt.plot(energy, background)
plt.title('Background Function')
plt.show(); plt.close()

def repeat_function(vector):
    return np.repeat([vector,vector],len(vector)*2,axis=0)


# if you make this perfectly correlated by just repeating the background function and taking corr, 
# you get a perfectly 1:1 correlation between all values - is this what we want?
# or do we want to sample noise around the background function? 
#bkg_cov = np.corrcoef(repeat_function(background))

random_bkg = [gaus_noise(bkg_sample, 0.01) for bkg_sample in repeat_function(background)]
bkg_cov = np.cov(random_bkg)

plt.matshow(bkg_cov)
cb = plt.colorbar()
cb.ax.tick_params(labelsize=10)
plt.title('Cov of bkg function', fontsize=10);

plt.show(); plt.close()
















