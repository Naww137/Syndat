


import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import syndat
import scipy.stats as stats
import matplotlib as mpl





sammy_directory = os.path.realpath('./synthetic_data/Ta181/')

samexp = pd.read_csv(os.path.join(sammy_directory,'SAMMY.LST'), sep='\s+', names=['E','exp_xs','exp_xs_unc','theo_xs','theo_xs_bayes','exp_trans','exp_trans_unc','theo_trans', 'theo_trans_bayes'])
samndf = pd.read_csv(os.path.join(sammy_directory,'SAMMY_endf.LST'), sep='\s+', names=['E','exp_xs','exp_xs_unc','theo_xs','theo_xs_bayes','exp_trans','exp_trans_unc','theo_trans', 'theo_trans_bayes'])

energy = samexp['E']

plt.rcParams['figure.dpi'] = 500
plt.plot(energy,samexp['theo_trans'], lw=3, label='$\sigma_{exp}$')
plt.plot(energy,samndf['theo_trans'], lw=0.5, c='r', label='$\sigma_{endf}$')

plt.legend()
plt.xlabel('Energy'); plt.ylabel('$\sigma$')
plt.yscale('log'); plt.xscale('log')
plt.show(); plt.close()

