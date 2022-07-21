#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 13:14:35 2022

@author: nwx
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import syndat
import scipy.stats as stat
import matplotlib as mpl
#import syndat.plot as myplt
import sys


sammy_directory =  os.path.realpath('./synthetic_data/Ta181')



#%% redux parms ---  experimentally unique values

# =============================================================================
# # constants and uncertainties
# =============================================================================
n = .12 # need to pull this from the endf evaluation
trig = 9760770# number of linac pulses
flux_mag = 1e5 # what is a reasonable flux magnitude??
detector_efficiency = 1
tof_dist = 35 # m   !!! not the proper tof_distance, but I need the erange to match up

# assuming 1% uncertainty of these values
m1 = 1; dm1 = m1*0.016
m2 = 1; dm2 = m2*0.008
m3 = 1; dm3 = m3*0.018
m4 = 1; dm4 = m4*0.005
alpha = [m1,m2,m3,m4]; d_alpha = [dm1,dm2,dm3,dm4]

a = 582.8061256946647; da = 1.14174241e+03
b = 0.0515158865500879; db = 2.18755273e-05
# a,b COV 1.42405866e-01

# sample in
k_i = 0.563; dk_i = k_i*0.0427
b0_i = 10.0; db0_i = b0_i*0.01
# sample out
K_o = 1.471; dK_o = K_o*0.0379
B0_o = 13.4; dB0_o = B0_o*0.01


# =============================================================================
# # open count rate 
# =============================================================================

# # estimate true underlying, raw, open count data with a wide gaussian flux
# cts_o_true = syndat.exp_effects.generate_open_counts(energy, flux_mag, 50, 100)

# or: import open count rate from RPI Ta-181 experiment:
C_o = pd.read_csv(os.path.join(sammy_directory,'ta181opencountrate.dat'), sep=',')
C_o['E'] = syndat.exp_effects.t_to_e(C_o['tof']*1e-6, tof_dist, True) 
# C_o['E'] = pd.read_csv(os.path.join(sammy_directory,'ta-181-12mm.twenty'), sep='\s+', names=['E','1','2'])['E']

ctr_o_true = C_o['co'][0:-5] # C_o.loc[C_o['E'] < 330, 'co']
tof = C_o['tof'][0:-5] # C_o.loc[C_o['E'] < 330, 'tof'] *1e-6  #!!! tof must be in microseconds to have the proper magnitude that Jesse normed to
energy = C_o['E'][0:-5]

# get ben width vector from tof
bw = np.flipud(np.diff(np.flipud(tof)))
bw = np.insert(bw,0,bw[0]) # assumes given leading tof edge

# background function - Jesse has normalized alread, tof mus be micro-seconds
def bkg_func(ti,a,b):
    return a*np.exp(ti*-b)
Bi = bkg_func(tof,a,b)


# =============================================================================
# 
# #plt.plot(C_o['tof'], C_o['co'])
# plt.plot(tof, ctr_o_true)
# plt.xlabel('tof'); plt.ylabel('countrate')
# plt.yscale('log'); plt.xscale('log')
# plt.show(); plt.close()
# =============================================================================


# get open counts - no need to propagate uncertainty bc we are un-reducing - taking these input values as truth
cts_o_true = ctr_o_true*bw*trig


#%% read sammy lst for experimentally corrected theoretical cross section

sam = syndat.sammy_interface.readlst(os.path.join(sammy_directory,'SAMMY.LST'))

# could test that energy structure lines up the same
# energy = sam['E'] 

T_theo = sam['theo_trans']

#%%
# =============================================================================
# # generate noisy, raw, sample in count data with statistical unc from a true underlying transmission
# =============================================================================

noisy_cts_i, noisy_cts_i_se, cts_i, cts_i_se = syndat.exp_effects.generate_raw_count_data(energy, T_theo, cts_o_true, bw, trig, k_i,K_o, Bi, b0_i,B0_o, alpha)

# take a noisy measurement of raw open count data with uncertainty
cts_o = cts_o_true # syndat.exp_effects.generate_open_counts(energy, flux_mag, 50, 100)
cts_o_se = np.sqrt(cts_o) # statistical error from poisson counting stats
noisy_cts_o = syndat.exp_effects.gaus_noise(cts_o, cts_o_se)
noisy_cts_o_se = np.sqrt(noisy_cts_o)


#%% Plotting only

# this section is for plotting only, the next section does these calculations under the hood
# but we want out count rates for plotting

ctr_o, ctr_o_se = syndat.exp_effects.cts_to_ctr(cts_o, cts_o_se, bw, trig) 
ctr_i, ctr_i_se = syndat.exp_effects.cts_to_ctr(cts_i, cts_i_se, bw, trig)
Tn = syndat.exp_effects.transmission(ctr_i,ctr_o, Bi, k_i,K_o, b0_i,B0_o, alpha)

plt.plot(tof, ctr_i, lw=1, label='cr_s')
plt.plot(tof, ctr_o, lw=1, label='cr_o')
plt.plot(tof, T_theo, label='Theoretical T', c='k', alpha=0.5)

plt.legend()
plt.ylim([1e-3,1e5])
#plt.xlim([2e2, 2.5e2])
plt.xscale('log'); plt.yscale('log')

plt.title('Count Rate comparison with overlayed theoretical transmission')
plt.show(); plt.close()


#%% Reduce the noisy, raw count data

# compile systematic uncertainties
sys_unc = np.append([da,db,dk_i,dK_o,db0_i,dB0_o], d_alpha)
# reduce raw, noisy count data with statistical uncertainties to transmission data with propagated uncertainty
Tn, dT, CovT = syndat.exp_effects.reduce_raw_count_data(tof, noisy_cts_i,noisy_cts_o, noisy_cts_i_se,noisy_cts_o_se, \
                                                        bw, trig, a,b, k_i,K_o, Bi, b0_i,B0_o, alpha, sys_unc)
# calculate the correlation matrix
CorT = np.corrcoef(CovT)

#%%

plt.matshow(CorT)
ax = plt.gca()

# ax.set_xticklabels(round(ax.get_xticks()), rotation = 45)
cb = plt.colorbar()
plt.title('Correlation')

plt.show(); plt.close()

#%%

syndat.plot.exp_theo(tof, Tn, dT, T_theo)
syndat.plot.unc_noise(tof, dT, T_theo, Tn, ctr_o)







