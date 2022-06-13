#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 11:08:47 2022

@author: noahwalton
"""


import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.stats as stats


def sample_chisquare(N_samples, DOF):
    samples = np.random.chisquare(DOF, size=N_samples)
    return samples
    
def chisquare_PDF(x, DOF, avg_reduced_width_square):
    x = x/avg_reduced_width_square
    y = stats.chi2.pdf(x, DOF)
    y_norm = y/avg_reduced_width_square
    return y_norm

def sample_resonance_widths(DOF, N_levels, avg_reduced_width_square):
    
    reduced_widths_square = avg_reduced_width_square*sample_chisquare(N_levels, DOF)
    partial_widths = 0  # add function with penetrability =2*P(E)*red_wid_sqr
    
    return reduced_widths_square, partial_widths


def compare_pdf_to_samples(reduced_widths_square_vector, avg_reduced_width_square, dof):
    
    fig = plt.figure(num=1,frameon=True); ax = fig.gca()
    
    x = np.linspace(0,max(reduced_widths_square_vector),10000)
    plt.plot(x, chisquare_PDF(x,dof,avg_reduced_width_square), color='r', label='$\chi^2$ PDF', zorder=10)
    
    if avg_reduced_width_square != 1:
        print(); print('WARNING: ')
        print('pdf has not been transformed for <D> other than 1 - will not match samples'); print()
        
    plt.hist(reduced_widths_square_vector, bins=75, density=True, ec='k', linewidth=0.75,color='cornflowerblue', zorder=2, label='samples')
    
    ax.set_facecolor('whitesmoke'); ax.grid(color='w', linestyle='-', linewidth=2, zorder=1, alpha=1)
    ax.set_xlabel('Reduced Widths Squared ($\gamma^2$)'); ax.set_ylabel('Normalized Frequency'); plt.title('Reduced Widths Squared ($\gamma^2$)')
    plt.legend()
    plt.show(); plt.close()
    
    return
