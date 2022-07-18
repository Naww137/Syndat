
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

import scipy.stats as stats
import matplotlib as mpl





def sample_poisson(m):
    return np.random.default_rng().poisson(lam=m)

def sample_gaussian(m):
    return np.random.default_rng().normal(loc=m,scale=np.sqrt(m))

m=1e2

# =============================================================================
# test = sample_gaussian(m)
# print(test)
# print(round(test))
# =============================================================================

s = int(1e4)

gaus = []; pois = []
for i in range(s):
    gaus.append(round(sample_gaussian(m)))
    #gaus.append(sample_gaussian(m))
    pois.append(sample_poisson(m))
    
test_dist = gaus
 

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
# plot the observed frequency vs expected from poisson
# =============================================================================

plt.scatter(mid_bin,lumped_obs_freq, c='b', s=1, label='Observed Frequency')
plt.plot(mid_bin,lumped_exp_freq,  c='orange',label='Expected Frequency')
plt.legend()

plt.show();plt.close()



lumped_obs_freq = np.array(lumped_obs_freq)
lumped_exp_freq = np.array(lumped_exp_freq)

chi2_stat = np.sum((lumped_obs_freq-lumped_exp_freq)**2/lumped_exp_freq)
crit_val = stats.chi2.ppf(1-0.05, df=len(lumped_exp_freq)-1)

print()
print(f'Chi-square statistic\n{chi2_stat}')
print(f'Chi-square critical value for dof=N-1 & sig=0.05\n{crit_val}')
print()
if chi2_stat < crit_val:
    print('Fail to reject null hypothesis')
elif chi2_stat > crit_val:
    print('Reject null hypothesis')


#%%
# =============================================================================
# perform chisquare goodness of fit test
# =============================================================================

result = stats.chisquare(lumped_obs_freq, f_exp=lumped_exp_freq, ddof=0)
stat = result[0]; p = result[1]

print(stats.chi2.ppf(1-0.05, df=len(lumped_exp_freq)-1))

print()
print(stat); print(p)
if stat > p:
    print('Reject null hypothesis: The observed data does not represent the expected')
elif stat < p:
    print('Fail to reject null hypothesis: The observed data may represent the expected')
print()





#%%
# =============================================================================
# 
# mybins = np.append(possible_values,possible_values[-1]+1)
# 
# pois_func = stats.poisson.pmf(possible_values, m, loc=0)
# 
# ghist = plt.hist(gaus, bins=mybins, density=True, ec='k', linewidth=0.75,color='cornflowerblue', label='gauss')
# phist = plt.hist(pois, bins=mybins, density=True, ec='k', linewidth=0.75,color='pink', label='poisson')
# plt.legend()
# 
# gn = ghist[0]; gb = ghist[1]
# pn = phist[0]; pb = phist[1]
# 
# #%%
# 
# result = stats.chisquare(gn, f_exp=pois_func, ddof=0, axis=0)
# stat = result[0]; p = result[1]
# 
# print(stats.chi2.ppf(1-0.05, df=len(possible_values)-1))
# 
# print()
# print(stat); print(p)
# if stat > p:
#     print('Reject null hypothesis: The observed data does not represent the expected')
# elif stat < p:
#     print('Fail to reject null hypothesis: The observed data may represent the expected')
# print()
# 
# =============================================================================



# =============================================================================
# shapiro_test = stats.shapiro(pois)
# 
# p = shapiro_test.pvalue
# sts = shapiro_test.statistic
# print(); print(f'{p}')
# if p < 0.05:
#     print('Reject null hypothesis: sample does not look normal')
# elif p >0.05:
#     print('Fail to reject null hypothesis: sample looks normal')
# print()
# 
# result = stats.anderson(pois)
# print(f'{result.statistic}')
# print(f'{result.critical_values}')
# if result.statistic > result.critical_values[2]:
#     print("Reject Null Hypothesis - Anderson")
# print()
# =============================================================================
