


import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import syndat
import scipy.stats as stats
import matplotlib as mpl




#%% need to pull estructure from the open count rate file
# =============================================================================
# 
# plt.rcParams['figure.dpi'] = 500
# 
# sammy_directory =  os.path.realpath('./synthetic_data/Ta181')
# 
# sam = pd.read_csv(os.path.join(sammy_directory,'SAMMY.LST'), sep='\s+', names=['E','exp_xs','exp_xs_unc','theo_xs','theo_xs_bayes','exp_trans','exp_trans_unc','theo_trans', 'theo_trans_bayes'])
# #samndf = pd.read_csv(os.path.join(sammy_directory,'SAMMY_endf.LST'), sep='\s+', names=['E','exp_xs','exp_xs_unc','theo_xs','theo_xs_bayes','exp_trans','exp_trans_unc','theo_trans', 'theo_trans_bayes'])
# energy = np.array(sam['E']);  tof = syndat.exp_effects.e_to_t(energy, 35, True)
# xs_theoretical = np.array(sam['theo_xs'])
# trans_theo = sam['theo_trans']
# 
# plt.rcParams['figure.dpi'] = 500
# plt.plot(energy,sam['theo_trans'], lw=1, label='$\sigma_{exp}$')
# #plt.plot(tof,sam['theo_trans'], lw=1, label='$\sigma_{exp}$')
# #plt.plot(energy,samndf['theo_trans'], lw=0.5, c='r', label='$\sigma_{endf}$')
# plt.legend()
# plt.xlabel('E'); # plt.ylabel('$\sigma$')
# plt.yscale('log'); plt.xscale('log')
# #plt.show(); plt.close()
# 
# 
# 
# #% import Open Count rate
# 
# # # estimate true underlying, raw, open count data with a wide gaussian flux
# # cts_o_true = syndat.exp_effects.generate_open_counts(energy, flux_mag, 50, 100)
# 
# # or: import open count rate from RPI Ta-181 experiment:
# C_o = pd.read_csv(os.path.join(sammy_directory,'ta181opencountrate.dat'), sep=',')
# # =============================================================================
# # C_o_tof = np.array(C_o['tof'])*1e-6 # microseconds to s
# # C_o_Epts = syndat.exp_effects.t_to_e(C_o_tof, tof_dist, True) 
# # =============================================================================
# 
# #C_o['E'] = syndat.exp_effects.t_to_e(C_o['tof']*1e-6, tof_dist, True) 
# C_o['E'] = pd.read_csv(os.path.join(sammy_directory,'ta-181-12mm.twenty'), sep='\s+', names=['E','1','2'])['E']
# 
# ctr_o_true = C_o['co'][0:-5] # C_o.loc[C_o['E'] < 330, 'co']
# tof_ = C_o['tof'][0:-5]*1e-6 # C_o.loc[C_o['E'] < 330, 'tof'] *1e-6
# energy_ = syndat.exp_effects.t_to_e(tof_, 35, False) # C_o['E'][0:-5] 
# 
# #plt.plot(tof_, ctr_o_true)
# plt.plot(energy_, ctr_o_true, label='c_o')
# plt.xlabel('E (eV)'); #plt.ylabel('countrate')
# plt.yscale('log'); plt.xscale('log')
# 
# 
# #plt.xlim([1,500])
# plt.ylim([1e-4, 1e5])
# plt.legend()
# plt.show(); plt.close()
# # plt.plot(np.diff(C_o['tof']))
# =============================================================================

#%%


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


# =============================================================================
# experimentally unique values
# =============================================================================
n = .12 # need to pull this from the endf evaluation
trig = 9760770 # number of linac pulses
# flux_mag = 1e5 # what is a reasonable flux magnitude??
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
tof = C_o['tof'][0:-5]#*1e-6 # C_o.loc[C_o['E'] < 330, 'tof'] *1e-6
energy = C_o['E'][0:-5]

#plt.plot(C_o['tof'], C_o['co'])
plt.plot(tof, ctr_o_true)
plt.xlabel('tof'); plt.ylabel('countrate')
plt.yscale('log'); plt.xscale('log')
plt.show(); plt.close()
# plt.plot(np.diff(C_o['tof']))



#%% read sammy lst

sam = syndat.sammy_interface.readlst(os.path.join(sammy_directory,'SAMMY.LST'))
# energy = sam['E']

#tof = syndat.exp_effects.e_to_t(energy, tof_dist, True)

bw = np.flipud(np.diff(np.flipud(tof)))
bw = np.insert(bw,0,bw[0]) # assumes given leading tof edge

# background function
def bkg_func(ti,a,b):
    return a*np.exp(ti*-b)
Bi = bkg_func(tof,a,b) #*bw*trig*1e-2


plt.plot(tof,Bi*K_o+[B0_o]*len(tof))
#plt.plot(tof,Bi*k_i+b0_i)
plt.xlabel('tof')
plt.title('Background Function')
plt.xscale('log'); plt.yscale('log')
plt.show(); plt.close()

#%%
# get open counts - no need to propagate uncertainty bc we are un-reducing
# but, does uncertainty propagate backwards and some out the same?

cts_o_true = ctr_o_true*bw*trig

T_theo = sam['theo_trans']


#%%

# =============================================================================
#  generate raw count data
# =============================================================================
Cr, dCr = syndat.exp_effects.cts_to_ctr(cts_o_true, np.sqrt(cts_o_true), bw, trig) # cts_o/(bw*trig)

# calculate sample in count rate from theoretical transmission, bkg, m,k, and open count rate
[m1,m2,m3,m4] = alpha
cr = (T_theo*(m3*Cr - m4*K_o*Bi - B0_o) + m2*k_i*Bi + b0_i)/m1

#print(sum(T_theo - ( (m1*cr-m2*k_i*Bi-b0_i)/(m3*Cr-m4*K_o*Bi-B0_o))))

# calculate sample in counts, noise, and uncertainty
c = cr*bw*trig 
dc = np.sqrt(c)
# =============================================================================
# nc = gaus_noise(c,dc) # will create some negative counts, force to zero
# nc = np.where(nc<0, 0, nc) # replace negative counts with 0
# dnc = np.sqrt(nc)
# =============================================================================
#print(np.array([cr<0]).any())



#%% plot it

plt.scatter(tof, cr, s=1, label='c_s')
plt.scatter(tof, Cr, s=1, label='c_o')
#plt.plot(tof, T_theo, label='trans', c='k')
plt.plot(tof, Bi*K_o+B0_o, c='k', label='Bi*K_o-B0_o')
plt.plot(tof, Bi*k_i+b0_i, c='r', alpha=0.5, label='Bi*K_s-B0_s')

plt.legend()
#plt.xlim([1e-4,1e-3]);
#plt.ylim([1e1,1e5])
plt.xscale('log'); plt.yscale('log')

#plt.title('Count comparison with overlayed theoretical transmission')

# =============================================================================
# plt.scatter(tof, ctr_i, s=1, label='cr_s')
# plt.scatter(tof, ctr_o, s=1, label='cr_o')
# plt.plot(tof, T_theo, label='Theoretical T', c='k', alpha=0.5)
# =============================================================================

plt.legend()
# =============================================================================
# plt.xlim([1.5e-4,2.1e-4]);
# plt.ylim([4e2,1e3])
# =============================================================================
plt.xscale('log'); plt.yscale('log')

plt.title('Count comparison with overlayed Background')
plt.xlabel('tof (s)');

plt.show(); plt.close()


#%%

# =============================================================================
# dc = noisy_cts_i_se
# dC = noisy_cts_o_se 
# 
# 
# # =============================================================================
# # dc = dc.flatten()
# # dC = dC.flatten()
# # =============================================================================
# 
# 
# 
# import time
# start = time.time()
# 
# dcC = np.append(dc,dC)
# Cov_stat1 = np.diag(dcC)
# 
# # =============================================================================
# # Cov_stat = np.zeros([len(tof)*2,len(tof)*2])
# # #samplein = True; sampleout = False
# # for i in range(len(tof)):
# #     for j in range(len(tof)):
# #         if i == j:
# #             Cov_stat[i,j] = dc[i] 
# # for i in range(len(tof),len(tof)*2):
# #     for j in range(len(tof),len(tof*2)):
# #         if i == j:
# #             Cov_stat[i,j] = dC[i]
# # 
# # =============================================================================
# 
# end = time.time()
# print(end - start)
# 
# =============================================================================

# =============================================================================
# Jac_stat = np.zeros([len(tof*2),len(tof*2)])
# for i in range(len(tof)):
#     for j in range(len(tof)):
#         if i == j:
#             Jac_stat[i,j] = dTi_dci(i)
# for i in range(len(tof),len(tof*2)):
#     for j in range(len(tof),len(tof*2)):
#         if i == j:
#             Jac_stat[i,j] = dTi_dCi(i)
# 
# # construct systematic covariance and jacobian
# 
# Cov_sys = np.zeros([len(sys_unc),len(sys_unc)])
# for i in range(len(sys_unc)):
#     for j in range(len(sys_unc)):
#         if i == j:
#             Cov_sys[i,j] = sys_unc[i]
# Cov_sys[0,1] = sys_unc[0]*sys_unc[1]  
# Cov_sys[1,0] = sys_unc[1]*sys_unc[0]        
# 
# Jac_sys = np.zeros([len(sys_unc),len(tof)])
# for i in range(len(sys_unc)):
#     for j in range(len(tof)):
#         Jac_sys[i,j] = dT_dsys(j)[i]
#             
# # calculate covariance of output
# CovT_stat = Jac_stat.T @ Cov_stat @ Jac_stat
# CovT_sys = Jac_sys.T @ Cov_sys @ Jac_sys
# 
# CovT = CovT_stat + CovT_sys
# 
# =============================================================================

 
 
