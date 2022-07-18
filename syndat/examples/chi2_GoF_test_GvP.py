#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 11:38:21 2022

@author: nwx
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats



def sample_poisson(m):
    return np.random.default_rng().poisson(lam=m)

def sample_gaussian(m):
    return np.random.default_rng().normal(loc=m,scale=np.sqrt(m))


def test_dist(test_dist, m, s, disp_bool):
    
    possible_counts = np.arange(min(test_dist), max(test_dist), 1)
    expected_frequencies = stats.poisson.pmf(possible_counts, m, loc=0)*s
    obs_hist = np.histogram(test_dist, bins=possible_counts)
    obs_freq = obs_hist[0]
    
    # =============================================================================
    # lump count bins where expected frequency is < 5
    # =============================================================================
    sumvec = np.cumsum(expected_frequencies)
    
    index_bin = []; index_bins = []
    previous_value = 0
    for index, value in enumerate(sumvec):
        if value < previous_value + 5:
            index_bin.append(index)
        if value > previous_value + 5:
            previous_value = value
            index_bins.append(index_bin)
            index_bin = [index]
            
    lumped_exp_freq = []; lumped_obs_freq = []
    lumped_count_bins = []; bin_edges = [min(possible_counts)]
    mid_bin = []
    for index in index_bins:
        lumped_exp_freq.append(sum(expected_frequencies[index]))
        lumped_obs_freq.append(sum(obs_freq[index]))
        lumped_count_bins.append(possible_counts[index])
        bin_edges.append(max(possible_counts[index])+1)
        mid_bin.append(np.median(possible_counts[index]))
     
        
    # =============================================================================
    # calculate chisquare statistic and compare
    # =============================================================================
    if disp_bool:
        plt.scatter(mid_bin,lumped_obs_freq, c='b', s=1, label='Observed Frequency')
        plt.plot(mid_bin,lumped_exp_freq,  c='orange',label='Expected Frequency')
        plt.legend()
        plt.show();plt.close()
    
    lumped_obs_freq = np.array(lumped_obs_freq)
    lumped_exp_freq = np.array(lumped_exp_freq)
    chi2_stat = np.sum((lumped_obs_freq-lumped_exp_freq)**2/lumped_exp_freq)
    crit_val = stats.chi2.ppf(1-0.05, df=len(lumped_exp_freq)-1)
    
    if chi2_stat < crit_val:
        reject = 0
    elif chi2_stat > crit_val:
        reject = 1
    
    if disp_bool:
        print()
        print(f'Chi-square statistic\n{chi2_stat}')
        print(f'Chi-square critical value for dof=N-1 & sig=0.05\n{crit_val}')
        print()
        if chi2_stat < crit_val:
            print('Fail to reject null hypothesis')
        elif chi2_stat > crit_val:
            print('Reject null hypothesis')
            
    return reject


def test_rounded_gaus(mean_num_counts, sample_size, disp_bool):
    
    m = mean_num_counts
    s = sample_size
    
    gaus = []; pois = []
    for i in range(s):
        gaus.append(round(sample_gaussian(m)))
        #gaus.append(sample_gaussian(m))
        pois.append(sample_poisson(m))
        
    reject_gaus = test_dist(gaus, m, s, disp_bool)
    reject_pois = test_dist(pois, m, s, disp_bool)
    
    return [reject_gaus, reject_pois]



mean_num_counts=1
sample_size = int(1e4)

mean_counts = [10, 1e2, 1e3, 1e4]

rejects = []
for each in mean_counts:
    rejects.append(test_rounded_gaus(each, sample_size, True))







