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
import matplotlib as mpl
#import syndat.plot as myplt
import sys


#%% need to pull estructure from the open count rate file

plt.rcParams['figure.dpi'] = 500

sammy_directory =  os.path.realpath('./synthetic_data/Ta181')

sam = pd.read_csv(os.path.join(sammy_directory,'SAMMY.LST'), sep='\s+', names=['E','exp_xs','exp_xs_unc','theo_xs','theo_xs_bayes','exp_trans','exp_trans_unc','theo_trans', 'theo_trans_bayes'])
#samndf = pd.read_csv(os.path.join(sammy_directory,'SAMMY_endf.LST'), sep='\s+', names=['E','exp_xs','exp_xs_unc','theo_xs','theo_xs_bayes','exp_trans','exp_trans_unc','theo_trans', 'theo_trans_bayes'])
energy = np.array(sam['E']);  tof = syndat.exp_effects.e_to_t(energy, 35, False)
xs_theoretical = np.array(sam['theo_xs'])
trans_theo = sam['theo_trans']

plt.rcParams['figure.dpi'] = 500
#plt.plot(energy,sam['theo_trans'], lw=1, label='$\sigma_{exp}$')
plt.plot(tof,sam['theo_trans'], lw=1, label='$\sigma_{exp}$')
#plt.plot(energy,samndf['theo_trans'], lw=0.5, c='r', label='$\sigma_{endf}$')
plt.legend()
plt.xlabel('tof'); plt.ylabel('$\sigma$')
plt.yscale('log'); plt.xscale('log')
plt.show(); plt.close()


#%%
# =============================================================================
# tof = 1.2320000000000002*1e-6 # 2363.1648*1e-6
# e = 763412.3553447116 #1.1619998393
# 
# 
# jev = 1.6022e-19 # J/eV
# mn = 1.674927498e-27 #kg
# #t = d/np.sqrt(E*jev*2/mn)
# 
# d = tof*np.sqrt(e*jev*2/mn)
# =============================================================================

#%%

     
#path_to_lst = '/Users/noahwalton/Library/Mobile Documents/com~apple~CloudDocs/Research Projects/Resonance Fitting/summer_2022/SAMMY_jb.LST'
#sammy_lst = pd.read_csv(path_to_lst, sep='\s+', names=['E','exp_xs','exp_xs_unc','theo_xs','theo_xs_bayes','exp_trans','exp_trans_unc','theo_trans','theo_trans_bayes'])           




# =============================================================================
# experimentally unique values
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
# a,b COV 
# =============================================================================
# [[1.14174241e+03 1.42405866e-01]
# [1.42405866e-01 2.18755273e-05]]
# =============================================================================

# sample in
k_i = 0.563; dk_i = k_i*0.0427
b0_i = 10.0; db0_i = b0_i*0.01
# sample out
K_o = 1.471; dK_o = K_o*0.0379
B0_o = 13.4; dB0_o = B0_o*0.01







#%% import Open Count rate

# # estimate true underlying, raw, open count data with a wide gaussian flux
# cts_o_true = syndat.exp_effects.generate_open_counts(energy, flux_mag, 50, 100)

# or: import open count rate from RPI Ta-181 experiment:
C_o = pd.read_csv(os.path.join(sammy_directory,'ta181opencountrate.dat'), sep=',')
# =============================================================================
# C_o_tof = np.array(C_o['tof'])*1e-6 # microseconds to s
# C_o_Epts = syndat.exp_effects.t_to_e(C_o_tof, tof_dist, True) 
# =============================================================================

#C_o['E'] = syndat.exp_effects.t_to_e(C_o['tof']*1e-6, tof_dist, True) 
C_o['E'] = pd.read_csv(os.path.join(sammy_directory,'ta-181-12mm.twenty'), sep='\s+', names=['E','1','2'])['E']

ctr_o_true = C_o['co'][0:-5] # C_o.loc[C_o['E'] < 330, 'co']
tof = C_o['tof'][0:-5]*1e-6 # C_o.loc[C_o['E'] < 330, 'tof'] *1e-6
energy = C_o['E'][0:-5]

#plt.plot(C_o['tof'], C_o['co'])
plt.plot(tof, ctr_o_true)
plt.xlabel('tof'); plt.ylabel('countrate')
plt.yscale('log'); plt.xscale('log')
plt.show(); plt.close()
# plt.plot(np.diff(C_o['tof']))


# syndat.sammy_interface.write_estruct_file(C_o_Epts, "./synthetic_data/Ta181/estruct")

#%%
# =============================================================================
# bleh = pd.read_csv(os.path.join(sammy_directory,'ta-181-12mm.twenty'), sep='\s+', names=['E','1','2'])
# 
# plt.plot(bleh['E'],C_o['co'])
# 
# plt.xscale('log'); plt.yscale('log')
# =============================================================================


#%% read sammy lst

sam = syndat.sammy_interface.readlst(os.path.join(sammy_directory,'SAMMY.LST'))
# energy = sam['E']

#tof = syndat.exp_effects.e_to_t(energy, tof_dist, True)

bw = np.flipud(np.diff(np.flipud(tof)))
bw = np.insert(bw,0,bw[0]) # assumes given leading tof edge

# background function
def bkg_func(ti,a,b):
    return a*np.exp(ti*-b)
Bi = bkg_func(tof,a,b)*1e-2


#%%
# get open counts - no need to propagate uncertainty bc we are un-reducing
# but, does uncertainty propagate backwards and some out the same?

cts_o_true = ctr_o_true*bw*trig


# plot counts vs count rate for open spectra
# =============================================================================
# plt.plot(tof,cts_o_true, label='cts')
# plt.plot(tof,cts_o_true/bw/trig, label='ctr')
# plt.xlabel('tof (s)')
# plt.yscale('log'); plt.xscale('log')
# 
# plt.legend()
# plt.show(); plt.close()
# 
# print()
# print("WARNING: Counts are fractional, must be integers")
# print()
# =============================================================================

#%% compare resolution broading on transmission vs cross section - sammy output of transmission vs cross section
# these might actually be the same thing, I bet that sammy still resolution broadens transmission, but then converts back to xs
# there is a flag in the sammy input to resolution broaden the xs rather than trans

T_theo = sam['theo_trans']

# =============================================================================
# T_theo_wrong = np.exp(-n*sam['theo_xs']); T_theo_wrong[np.where(T_theo_wrong <= 1e-10)[0]] = 0 # set very small values to zero
# 
# plt.plot(energy, T_theo, label='samcalc')
# plt.scatter(energy, T_theo_wrong, label='wrong', s=2, color='r')
# 
# plt.xscale('log'); plt.yscale('log')
# plt.legend()
# plt.show(); plt.close()
# 
# print(); print('The difference in resolution broadening transmission vs xs is:')
# print(sum(T_theo-T_theo_wrong))
# =============================================================================



#%%

# =============================================================================
# # generate noisy, raw, sample in count data with statistical unc from a true underlying transmission
# =============================================================================

noisy_cts_i, noisy_cts_i_se, cts_i, cts_i_se = syndat.exp_effects.generate_raw_count_data(energy, T_theo, cts_o_true, bw, trig, k_i,K_o, Bi, b0_i,B0_o, alpha)

if np.array([noisy_cts_i<0]).any():
    print('Warning: negative count values')
    sys.exit()
elif np.isnan(noisy_cts_i_se).any():
    print('Warning: nan in Noisy counts SE')
    
    

#%%

# take a noisey measurement of raw open count data with uncertainty
cts_o = cts_o_true # syndat.exp_effects.generate_open_counts(energy, flux_mag, 50, 100)
cts_o_se = np.sqrt(cts_o) # statistical error from poisson counting stats
noisy_cts_o = syndat.exp_effects.gaus_noise(cts_o, cts_o_se)
noisy_cts_o_se = np.sqrt(noisy_cts_o)


# calculate count rate and propagate uncertainty
ctr_o, ctr_o_se = syndat.exp_effects.cts_to_ctr(cts_o, cts_o_se, bw, trig) 
ctr_i, ctr_i_se = syndat.exp_effects.cts_to_ctr(cts_i, cts_i_se, bw, trig)

# calculate transmission
Tn = syndat.exp_effects.transmission(ctr_i,ctr_o, Bi, k_i,K_o, b0_i,B0_o, alpha)
    

    

#%% plot noisy counts and clean counts with errorbars

# =============================================================================
# plt.plot(tof,Tn)
# 
# plt.xscale('log'); plt.yscale('log')
# =============================================================================

# =============================================================================
# plt.errorbar(tof, noisy_cts_i, yerr=noisy_cts_i_se, label='in', color='b',ecolor='b',elinewidth=1,capsize=2, fmt='.', ms=3)
# plt.errorbar(tof, cts_i, yerr=cts_i_se, label='in', color='r',ecolor='k',elinewidth=1,capsize=2, fmt='.', ms=3)
# #plt.plot(tof,cts_i, c='b', ms=1)
# #plt.scatter(tof,noisy_cts_i, c='r', s=1)
# 
# plt.xscale('log'); plt.yscale('log')
# plt.xlim([1.5e-4,1.9e-4]);
# plt.ylim([1.3e2,1.4e2])
# =============================================================================



#%% plot generated counts and count rate for sample in/sample out


# =============================================================================
# plt.scatter(tof, cts_i, s=1, label='c_s')
# plt.scatter(tof, cts_o, s=1, label='c_o')
# #plt.plot(tof, T_theo, label='trans', c='k')
# #plt.plot(tof, Bi, label='Bi', c='k')
# 
# plt.legend()
# #plt.xlim([1e-4,1e-3]);
# plt.ylim([1e-2,1e7])
# plt.xscale('log'); plt.yscale('log')
# 
# #plt.title('Count comparison with overlayed theoretical transmission')
# 
# plt.scatter(tof, ctr_i, s=1, label='cr_s')
# plt.scatter(tof, ctr_o, s=1, label='cr_o')
# plt.plot(tof, T_theo, label='Theoretical T', c='k', alpha=0.5)
# 
# plt.legend()
# # =============================================================================
# # plt.xlim([1.5e-4,2.1e-4]);
# # plt.ylim([4e2,1e3])
# # =============================================================================
# plt.xscale('log'); plt.yscale('log')
# 
# plt.title('Count Rate comparison with overlayed theoretical transmission')
# plt.show(); plt.close()
# 
# 
# #plt.plot(tof,cts_o_true*np.exp(-n*xs_theoretical))
# =============================================================================

    
    
#%% plots open counts


# =============================================================================
# fig, (ax1,ax2) = plt.subplots(2,1,sharex=True, constrained_layout=True)
# ax1.plot(tof, Bi); plt.xlabel('tof'); ax1.set_ylabel('Bi(t)')
# ax2.plot(tof,cts_o_true); plt.xlabel('tof'); ax2.set_ylabel('open counts'); ax2.set_yscale('log');ax2.set_xscale('log')
# 
# plt.show(); plt.close()
# =============================================================================


#plot2(tof,T_theo,Tn,dT, 'Fully Correlated Uncertainty')



#%%


# =============================================================================
# plt.errorbar(tof,noisy_cts_o,yerr=noisy_cts_o_se, label='open', color='k',ecolor='k',elinewidth=1,capsize=2, fmt='.', ms=3)
# plt.errorbar(tof,noisy_cts_i,yerr=noisy_cts_i_se, label='in', color='b',ecolor='b',elinewidth=1,capsize=2, fmt='.', ms=3)
# 
# plt.xscale('log');plt.yscale('log')
# plt.legend()
# plt.xlabel('tof'); plt.ylabel('counts')
# =============================================================================


#%% reduce raw count data to Transmission and propagate uncertainty

# =============================================================================
# now reduce the raw count data - can use the same, or different reduction parameters, bkg/open counts
# =============================================================================

print(np.array([noisy_cts_i<0]).any())
#print(np.isnan(np.array(noisy_cts_o).any()))

print(np.isnan(np.array(noisy_cts_i).any()))

#%%
#systematic uncertainties
sys_unc = np.append([da,db,dk_i,dK_o,db0_i,dB0_o], d_alpha)

# reduce raw, noisy count data with statistical uncertainties to transmission data with propagated uncertainty
Tn, dT, CovT = syndat.exp_effects.reduce_raw_count_data(tof, noisy_cts_i,noisy_cts_o, noisy_cts_i_se,noisy_cts_o_se, \
                                                        bw, trig, a,b, k_i,K_o, Bi, b0_i,B0_o, alpha, sys_unc)


    
#%%


plt.scatter(tof, Tn, label='Experimental', s=1, c='k')
plt.plot(tof, T_theo, label='Theoretical', c='g', lw=0.25)

plt.legend()
plt.ylim([1e-3,1e1])
plt.xscale('log');plt.yscale('log')

#%%

fig, ax = plt.subplots(2,2, gridspec_kw={'height_ratios': [1, 1], 'width_ratios':[2,1]}) # , figsize=(12,5)
plt.rcParams['figure.dpi'] = 500
ax1 = ax[0,0]; ax2=ax[1,0]; ax3=ax[0,1]; ax4=ax[1,1]

ax1.scatter(tof, dT, lw=0.5, color='b', s=2, zorder=2)
ax1.set_ylim([0,10])
ax1.set_ylabel('$\delta$T'); ax1.set_xscale('log')#('$\sigma$')

rel_se = (Tn-T_theo)/T_theo
ax2.scatter(tof, (Tn-T_theo), s=2)
ax2.set_xlabel('ToF (s)'); ax2.set_ylabel('Noise')

# =============================================================================
# ax3.matshow(CovT)
# ax4.matshow(np.corrcoef(CovT))
# =============================================================================
plt.xscale('log')

plt.suptitle('Uncertainty and Noise on Transmission')
plt.tight_layout()
plt.show(); plt.close()













#%% Ignore below this: came from previous work, incorrect method


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











