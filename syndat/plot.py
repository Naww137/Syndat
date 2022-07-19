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