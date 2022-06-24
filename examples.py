#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 11:48:10 2022

@author: nwx
"""

import syndat
import os
import numpy as np
import matplotlib.pyplot as plt
# import nuctools as nuc
import pandas as pd



#%% create energy grid linear in ToF

jev = 1.6022e-19 # J/eV
mn = 1.674927498e-27 #kg

# min and max RRR energies for Cu-63 - in actuality we will want to go above and below these values a bit
RRR_speeds = np.array([np.sqrt(600*jev*2/mn), np.sqrt(300000*jev*2/mn)])

# 400 m flightpath at GELINA
distance = np.array([400, 400]) 

RRR_ToF = distance/RRR_speeds
dt = 1e-7 # linear time between neutron pts... this should also be better informed
time_points = np.arange(min(RRR_ToF),max(RRR_ToF),dt)

e_pts = syndat.exp_effects.t_to_e(time_points,400, True)
Estruct = np.flipud(e_pts)

plt.plot(Estruct,np.flipud(time_points))
plt.title('ToF vs Energy')
plt.xlabel('Energy (eV)'); plt.ylabel('ToF (s)')
plt.show(); plt.close()



#%% 



sammy_directory = os.path.realpath('./synthetic_data') #!!! give an appropriate directory


# =============================================================================
#   calculated average parameters for each spin group in endf defined in samndf.par file
# =============================================================================
average_parameters = syndat.sammy_interface.read_SAMNDF_PAR(os.path.join(sammy_directory,'SAMNDF.PAR'))




# =============================================================================
# # Other nuclear data
# =============================================================================
I = 1.5
i = 0.5
l_wave_max = 1
print_out = True
save_csv = False

RRR_Erange = [400, 300000]

# from mug
Davg_swave = 722
Davg_pwave = 404
Gavg_swave = 0.5
Gavg_pwave = 0.26

# synthetic
# =============================================================================
# Davg_neg = [1600]*5 #[700, 700, 700, 700, 700]
# Gavg_neg = [0.5, 0.5, 0.5, 0.5, 0.5]
# 
# Davg_pos = [700]*5 #[404, 404, 404, 404]
# Gavg_pos = [0.26, 0.26, 0.26, 0.26]
# =============================================================================
        

Davg = [list(average_parameters.dE[1:2]), list(average_parameters.dE[3:6])]
Ggavg = [list(average_parameters.Gg[1:2]), list(average_parameters.Gg[3:6])]
Gnavg = [list(average_parameters.Gn[1:2]), list(average_parameters.Gn[3:6])]


#%%

# =============================================================================
# # sample a single resonance ladder (all spin groups) and create necessary sammy files
# =============================================================================

Jn_ladders, Jp_ladders = syndat.spin_groups.sample_all_Jpi(I, i, l_wave_max,  
                    RRR_Erange, 
                    Davg, Ggavg, Gnavg, 
                    print_out,
                    save_csv, 
                    sammy_directory)

syndat.sammy_interface.create_sammyinp(os.path.join(sammy_directory,'sammy.inp'))
syndat.sammy_interface.create_sammypar(Jn_ladders, Jp_ladders,os.path.join(sammy_directory,'sammy.par'))
syndat.sammy_interface.write_estruct_file(Estruct, os.path.join(sammy_directory,'estruct'))


#%%

# Or! create sammy files for a number of resonance ladder realizations with 
# the Module for Multiple DataSet Acquisition (MMDA)

# =============================================================================
# case_directory= "/Users/nwx/work/synthetic_data/"  #!!! give an appropriate case directory
# number_of_realizations = 5
# 
# syndat.MMDA.wrapped_sammy_file_creator(number_of_realizations, case_directory, Estruct, \
#                                I, i, l_wave_max,
#                                RRR_Erange,  
#                                Davg, Gavg,  
#                                Gavg_swave,  
#                                print_out,  
#                                save_csv)
# =============================================================================



#%%


# somehow, somehwere, somebody can run the each of the created sammy's 
# in the realization directories <case_directory>/realization_#/

# I have set the directory structure to read a sammy LST I have run and included in the github repo



#%%

# Then we read out the LST file as the theoretical, experimentally corrected cross section
# from their we add experimental noise to each data point

sammy_lst = pd.read_csv(os.path.join(sammy_directory,'SAMMY_high.LST'), sep='\s+', names=['E','exp','exp_unc','xs','xs_bayes','1','2','3'])

energy = sammy_lst['E']
xs_theoretical = sammy_lst['xs']

noise = np.random.default_rng().normal(loc=0.0, scale=np.mean(xs_theoretical)*1e-4, size=len(xs_theoretical))
xs_experimental = xs_theoretical + noise


plt.rcParams['figure.dpi'] = 500
plt.plot(energy,xs_theoretical, lw=0.5, label='$\sigma_{exp}$')
plt.scatter(energy,xs_experimental, s=0.1, c='r', label='$\sigma_{exp}$')

plt.legend()
plt.xlabel('Energy'); plt.ylabel('$\sigma$')
plt.yscale('log'); plt.xscale('log')
plt.show(); plt.close()





