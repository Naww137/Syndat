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
from matplotlib.pyplot import *


sammy_directory =  os.path.realpath('../synthetic_data/Ta181')

#'/Users/nwx/Documents/GitHub/nuc_syndat/synthetic_data/Ta181'

#"/Users/noahwalton/Documents/GitHub/nuc_syndat/synthetic_data/Ta181"
# #os.path.realpath('./synthetic_data/Ta181')





#%% redux parms ---  experimentally unique values

# =============================================================================
# # constants and uncertainties
# =============================================================================
n = .12 # need to pull this from the endf evaluation
trig_o = 9760770 # number of linac pulses for open 
# trig_i = 18479997 
flux_mag = 1e5 # what is a reasonable flux magnitude??
detector_efficiency = 1
tof_dist = 35.185 # m   !!! not the proper tof_distance, but I need the erange to match up
t0 = 3.326

# assuming 1% uncertainty of these values
m1 = 1; dm1 = m1*0.016
m2 = 1; dm2 = m2*0.008
m3 = 1; dm3 = m3*0.018
m4 = 1; dm4 = m4*0.005
alpha = [m1,m2,m3,m4]; d_alpha = [dm1,dm2,dm3,dm4]

a = 582.8061256946647; da = np.sqrt(1.14174241e+03)
b = 0.0515158865500879; db = np.sqrt(2.18755273e-05)
# a,b COV 1.42405866e-01

# sample in
k_i = 0.563; dk_i = k_i*0.0427
# b0_i = 10.0; db0_i = b0_i*0.01
b0_i,db0_i = 9.9,0.1
# sample out
K_o = 1.471; dK_o = K_o*0.0379
# B0_o = 13.4; dB0_o = B0_o*0.01
B0_o,dB0_o = 13.4,0.7


#%%

# =============================================================================
# true open count rate 
# =============================================================================
# # estimate true underlying, raw, open count data with a wide gaussian flux
# cts_o_true = syndat.exp_effects.generate_open_counts(energy, flux_mag, 50, 100)

# or: import open count rate from RPI Ta-181 experiment:
odat = pd.read_csv(os.path.join(sammy_directory,'rpi-open-ta181.csv'), sep=',') #pd.read_csv(os.path.join(sammy_directory,'ta181opencountrate.dat'), sep=',')
odat = odat[odat.tof >= t0]
odat.reset_index(drop=True, inplace=True)
odat['E'] = syndat.exp_effects.t_to_e((odat.tof+t0)*1e-6, tof_dist, True) 

# cts_o_true = odat.counts[4:]
# tof = odat.tof[4:]  #!!! tof must be in microseconds to have the proper magnitude that Jesse normed to
# energy = odat.E[4:]

# get bin width vector from tof o
# bw = np.flipud(np.diff(np.flipud(tof)))
# bw = np.insert(bw,0,bw[0]) # assumes given leading tof edge
odat['bw'] = odat.bin_width*1e-6 # must we convert to microseconds
# bw = odat.bw[4:]

# background function - Jesse has normalized alread, tof mus be micro-seconds
def bkg_func(ti,a,b):
    return a*np.exp(ti*-b)
Bi = bkg_func(odat.tof,a,b)
tof=odat.tof


# syndat.sammy_interface.write_estruct_file(energy, "/Users/noahwalton/Library/Mobile Documents/com~apple~CloudDocs/Research Projects/Resonance Fitting/summer_2022/Ta181_sammy/estruct")

#%% read sammy lst for experimentally corrected theoretical cross section

sam = syndat.sammy_interface.readlst(os.path.join(sammy_directory,'SAMMY.LST'))
# sam = pd.read_csv("/Users/noahwalton/research_local/resonance_fitting/summer_2022/ogsam.csv", sep=',', names=['E','exp_xs','exp_xs_unc','theo_xs','theo_xs_bayes','exp_trans','exp_trans_unc','theo_trans', 'theo_trans_bayes'])
# could test that energy structure lines up the same
# energy = sam['E'] 

T_theo = np.flipud(sam.theo_trans)
sdat = pd.DataFrame()
sdat['theo_trans'] = sam.theo_trans #T_theo
sdat['E'] = sam.E
sdat['tof'] = syndat.exp_effects.e_to_t(sdat.E, tof_dist, True)*1e6+t0
sdat.sort_values('tof', axis=0, ascending=True, inplace=True)
sdat.reset_index(drop=True, inplace=True)

plt.figure()
plt.plot(sdat.tof,sdat.theo_trans)
plt.xscale('log'); plt.yscale('log')

# could test that energy structure lines up the same
# =============================================================================
# # energy = sam['E'] 
# sdat = pd.DataFrame()
# sdat['E'] = sam.E
# sdat['theo_trans'] = sam.theo_trans
# sdat['tof'] = syndat.exp_effects.e_to_t(sdat.E, tof_dist, True)*1e6+t0
# sdat.sort_values('tof', axis=0, inplace=True)
# 
# #T_theo = np.flipud(sam.theo_trans)
# #sdat['theo_trans'] = T_theo
# 
# plt.plot(sdat.tof, sdat.theo_trans)
# plt.xscale('log'); plt.yscale('log')
# =============================================================================

#%%


sdat, odat = syndat.exp_effects.generate_raw_count_data(sdat, odat, trig_o, k_i,K_o, Bi, b0_i, B0_o, alpha)


plt.figure()
plt.plot(tof,sdat.c, label='cts in'); 
plt.plot(tof,odat.counts, label='cts out')
plt.plot(tof, k_i*Bi+b0_i, label='bkg_i')
plt.plot(tof, K_o*Bi+B0_o, label='bkg_o')
plt.xscale('log'); plt.yscale('log')
plt.legend()


#%% Plotting only

# this section is for plotting only, the next section does these calculations under the hood
# but we want out count rates for plotting

# calculated sample in count rate clean and noisey

sdat['cps'], sdat['dcps'] = syndat.exp_effects.cts_to_ctr(sdat.c, sdat.dc, odat.bw, trig_o)
sdat['ncps'], sdat['dncps'] = syndat.exp_effects.cts_to_ctr(sdat.nc, sdat.dnc, odat.bw, trig_o)

figure()
plt.plot(tof,odat.cps,label='cps_o')
plt.plot(tof,sdat.ncps, label='cps_i')
xscale('log'); yscale('log')
legend()



#%% Reduce the noisy, raw count data

# compile systematic uncertainties #[dks,dko,db0s,db0o,dalpha1,dalpha2,dalpha3,dalpha4]
sys_unc = np.append([da,db,dk_i,dK_o,db0_i,dB0_o], d_alpha)
# %time
# #%%timeit
# # reduce raw, noisy count data with statistical uncertainties to transmission data with propagated uncertainty
Tn, dT, CovT = syndat.exp_effects.reduce_raw_count_data(tof, sdat.nc, odat.counts, sdat.dnc, np.sqrt(odat.counts), \
                                                        odat.bw, trig_o, a,b, k_i,K_o, Bi, b0_i,B0_o, alpha, sys_unc)
# # calculate the correlation matrix
CorT = np.corrcoef(CovT)

#%% reduce clean raw count data to see if I can reconstruct theoretical transmission

# compile systematic uncertainties #[dks,dko,db0s,db0o,dalpha1,dalpha2,dalpha3,dalpha4]
sys_unc = np.append([da,db,dk_i,dK_o,db0_i,dB0_o], d_alpha)
# %time
# #%%timeit
# # reduce raw, noisy count data with statistical uncertainties to transmission data with propagated uncertainty
T_clean, dT_clean, CovT_clean = syndat.exp_effects.reduce_raw_count_data(tof, sdat.c, odat.counts, sdat.dc, np.sqrt(odat.counts), \
                                                        odat.bw, trig_o, a,b, k_i,K_o, Bi, b0_i,B0_o, alpha, sys_unc)

    #%%
    
figure()
plot(tof,T_clean, label='Reconstructed')
plot(tof,sdat.theo_trans, label='theo')
legend()
xscale('log')
xlim([1e2,1e3])
ylim([-0.05,0.6])
axhline(y=0.0, ls='--')
#yscale('log')
show(); close()


#%%


figure()
plot(tof,np.sqrt(np.diag(CovT)))
xlim([1e1,2e3])
xscale('log'); yscale('log')
plt.show(); plt.close()

#%%
# figure()
matshow(CorT)
colorbar()
title('Correlation')
plt.show(); plt.close()

#%%

#syndat.plot.exp_theo(tof, Tn, dT, T_theo)
#syndat.plot.unc_noise(tof, dT, T_theo, Tn, ctr_o, ctr_i)

fig, (ax1, ax3) = plt.subplots(2, gridspec_kw={'height_ratios': [1, 1]}, sharex=True) # , figsize=(12,5)
# plt.rcParams['figure.dpi'] = 500
# ax1 = ax[0,0]; ax2=ax[1,0]; ax3=ax[0,1]; ax4=ax[1,1]

ax1.scatter(tof, dT/Tn*100, lw=0.5, color='b', s=0.5, zorder=2)
#ax1.set_ylim([0,2])
ax1.set_yscale('log')
ax1.set_ylabel('% Unc'); #('$\sigma$')('$\delta$T')

# ax2.plot(tof,ctr_o, lw= 0.5, c='orange', label=r'$ctr_{out}$')
# ax2.plot(tof,ctr_i, lw= 0.5, c='cornflowerblue', label=r'$ctr_{in}$')
# ax2.plot(tof,T_theo, lw= 0.5, c='g', label=r'$T_{theo}$')
# ax2.set_yscale('log')
# ax2.set_ylim([1e-5,3e5])
# ax2.set_ylabel(r'$ctr_o$')
# ax2.legend()
# =============================================================================
#     ax25 = plt.twinx(ax2)
#     ax25.plot(tof,T_theo, lw= 0.5, c='g', label=r'$T_{theo}$')
#     ax25.set_ylabel(r'$T_{theo}$')
#     ax25.set_yscale('log')
#     ax25.set_ylim([1e-3,1e1])
# =============================================================================

rel_se = (Tn-T_theo)/T_theo
ax3.scatter(tof, rel_se, s=0.5, c='b')
#ax3.set_ylim([0,1])
ax3.set_ylabel('Relative Noise')
ax3.set_yscale('log')

#plt.xlim([1e2,2e3])

plt.xscale('log')
plt.xlabel('ToF (s)');
plt.suptitle('Uncertainty and Noise on Transmission')
plt.tight_layout()
#lt.show(); plt.close()






