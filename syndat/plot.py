#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 13:16:02 2022

@author: nwx
"""

import matplotlib.pyplot as plt
import numpy as np



def plot1(energy,theo,exp,label1,label2):
    
    plt.plot(energy, theo, label=label1, zorder=2)
    plt.scatter(energy,exp, label=label2, s=1, c='k', zorder=1)
    
    plt.legend()
    #plt.yscale('log'); 
    plt.xscale('log')
    plt.show();plt.close()
    

def plot2(x,theo,exp,exp_unc, title):
    
    fig, (ax1,ax2,ax3) = plt.subplots(3,1, sharex=True, constrained_layout=True, gridspec_kw={'height_ratios': [2, 1, 1]}) # , figsize=(12,5)
    plt.rcParams['figure.dpi'] = 500
    
    ax1.plot(x,theo, lw=0.5, color='b', label='$T_{theo}$', zorder=2)
    #ax1.scatter(energy,exp, s=0.1, c='r', label='$T_{exp}$')
    ax1.errorbar(x, exp, yerr=exp_unc, color='k',ecolor='k',elinewidth=1,capsize=2, fmt='.', ms=3, label='$T_{exp}$', zorder=0)
    
    ax1.legend()
    ax1.set_ylabel('T') #('$\sigma$')
    #ax1.set_yscale('log'); 
    #ax1.set_xscale('log')
    ax1.set_ylim([0,max(exp)+0.1])
    
    rel_se = np.sqrt((exp-theo)**2) #/theo
    ax2.scatter(x, rel_se, s=2)
    #ax2.set_ylim([-.5,.5])
    ax2.set_ylabel('L2 Norm'); #ax2.set_ylabel('L1 Norm (relative)')
    
    ax3.scatter(x, exp_unc, lw=0.5, color='b', s=2, zorder=2)
    ax3.set_ylabel('$\delta$T') #('$\sigma$')
    ax3.set_xlabel('ToF (s)');
    
    plt.suptitle(title)
    plt.tight_layout()
    plt.show(); plt.close()
    
    
def exp_theo(tof, Tn, dT, T_theo):
    plt.errorbar(tof,Tn, yerr=dT,color='r',ecolor='k',elinewidth=1,capsize=2, fmt='.', ms=3)
    #plt.scatter(tof, Tn, label='Experimental', s=1, c='k')
    plt.plot(tof, T_theo, label='Theoretical', c='g', lw=0.25)
    
    plt.legend()
    #plt.ylim([1e-5,1e1])
    #plt.xlim([1e2,1e3])
    plt.xscale('log');plt.yscale('log')
    

def unc_noise(tof, dT, T_theo, Tn):
    #fig, ax = plt.subplots(2,2, gridspec_kw={'height_ratios': [1, 1], 'width_ratios':[2,1]}) # , figsize=(12,5)
    fig, (ax1, ax2, ax3) = plt.subplots(3, gridspec_kw={'height_ratios': [1, 1, 1]}, sharex=True) # , figsize=(12,5)
    plt.rcParams['figure.dpi'] = 500
    # ax1 = ax[0,0]; ax2=ax[1,0]; ax3=ax[0,1]; ax4=ax[1,1]
    
    ax1.scatter(tof, dT, lw=0.5, color='b', s=1, zorder=2)
    #ax1.set_ylim([0,2])
    ax1.set_yscale('log')
    ax1.set_ylabel('$\delta$T'); #('$\sigma$')
    
    ax2.plot(tof,T_theo, lw= 0.5, c='g')
    ax2.set_ylabel(r'$T_{theo}$')
    ax2.set_yscale('log')
    #ax2.set_ylim([1e-10,1e1])
    
    rel_se = (Tn-T_theo)/T_theo
    ax3.scatter(tof, rel_se, s=1)
    #ax3.set_ylim()
    ax3.set_ylabel('Noise')
    ax3.set_yscale('log')
    
    #plt.xlim([1e2,2e3])
    plt.xscale('log')
    plt.xlabel('ToF (s)');
    plt.suptitle('Uncertainty and Noise on Transmission')
    plt.tight_layout()
    plt.show(); plt.close()