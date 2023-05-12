
# %%
import numpy as np
import syndat
import pandas as pd
import os
from matplotlib.pyplot import *
import h5py


# %%
ac = 0.81271  # scattering radius in 1e-12 cm 
M = 180.948030  # amu of target nucleus
m = 1           # amu of incident neutron
I = 3.5         # intrinsic spin, positive parity
i = 0.5         # intrinsic spin, positive parity
l_max = 0       # highest order l-wave to consider


Ta_pair = syndat.particle_pair( ac, M, m, I, i, l_max,
                                input_options={})

# create an energy domain, min/max
E_min_max = [100, 120]

### or just give the min/max and the experiment object will do the above
energy_grid = E_min_max

input_options = {'Add Noise': True,
                'Calculate Covariance': True,
                'Compression Points':[],
                'Grouping Factors':None}

experiment_parameters = {'bw': {'val':0.3,    'unc'   :   0}}

# initialize experimental setup
exp = syndat.experiment(energy_grid, 
                        input_options=input_options, 
                        experiment_parameters=experiment_parameters)


spin_groups = [ (3.0,1,0)] # , (4.0,1,[0]) ]

average_parameters = pd.DataFrame({ 'dE'    :   {'3.0':8.79, '4.0':4.99},
                                    'Gg'    :   {'3.0':46.4, '4.0':35.5},
                                    'gn2'    :   {'3.0':64.0, '4.0':64.0}  })
                                    
resonance_ladder = pd.DataFrame({'E':[120], 'Gg':[0], 'gnx2':[1], 'J':[3], 'chs':[1], 'lwave':[[0]], 'J_ID':[None]})

xs_tot, xs_scat, xs_cap = syndat.scattering_theory.SLBW(exp.energy_domain, Ta_pair, resonance_ladder)

# convert to transmisison and put in an appropriate dataframe
n = 0.067166 # atoms per barn or atoms/(1e-12*cm^2)
trans = np.exp(-n*xs_tot)
trans_thick = np.exp(-0.5*xs_tot)
theoretical_df = pd.DataFrame({'E':exp.energy_domain, 'theo_trans':trans})

# exp.run(theoretical_df)


# case_file = './MMDA_data'  # if NOT using hdf5
case_file = './MMDA_data.hdf5'  # if using hdf5

dataset_range = (0, 10)

spin_groups = [ (3.0,1,0) ]
# an alternative option would be to give Ta_pair.J, as long as you give an average parameter dataframe with corresponding indices
# spin_groups = Ta_pair.J
Ta_pair = syndat.particle_pair( ac, M, m, I, i, l_max,
                                input_options={},
                                spin_groups=spin_groups,
                                average_parameters=average_parameters )    

vary_Erange = {'fullrange':(3,1000), 'maxres':5 , 'prob':0.01}

samples_not_generated = syndat.MMDA.generate(Ta_pair, exp, 
                                        'syndat_SLBW', 
                                        dataset_range, 
                                        case_file,
                                        fixed_resonance_ladder=None, 
                                        open_data=None,
                                        vary_Erange=vary_Erange,
                                        use_hdf5=True,
                                        overwrite = True
                                                                    )

# vary_Erange.keys()

# %%
# print(pd.read_hdf(case_file, 'sample_2/exp_pw'))
with h5py.File(case_file, 'a') as f:
    print(f['sample_3/exp_cov'][()])
    f.close()

# %%
test = pd.read_hdf(case_file, 'sample_2/exp_pw')
figure()
plot(test.E, test.theo_trans)

# %%


# %%



