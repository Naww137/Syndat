
# %%
import numpy as np
import syndat
import pandas as pd
import os
from matplotlib.pyplot import *



# %% Define reaction particle-pair

ac = 0.81271    # scattering radius in 1e-12 cm 
M = 180.948030  # amu of target nucleus
m = 1           # amu of incident neutron
I = 3.5         # intrinsic spin, positive parity
i = 0.5         # intrinsic spin, positive parity
l_max = 1       # highest order l-wave to consider

Ta_pair = syndat.particle_pair(ac, M, m, I, i, l_max)

Ta_pair.map_quantum_numbers(False)

# %%  define energy grid

# linear in energy
# energy_grid = np.linspace(1,100,1000)

# linear in tof
E_min_max = [10, 1000]
tof_min_max = syndat.exp_effects.e_to_t(np.array(E_min_max),35.185, True)
bin_width = 0.1e-6
tof_grid = np.arange(min(tof_min_max), max(tof_min_max), bin_width)
energy_grid = syndat.exp_effects.t_to_e(tof_grid,35.185,True)

# %% Define spin groups and average parameters

# spin_groups variable here is equivalent to Ta_pair.J
# spin_groups = [ (3.0,1,[0]), (4.0,1,[0]), (-4.0,2,[1.0, 1.0])]
spin_groups = [ (3.0,1,[0])]

average_parameters = pd.DataFrame({ 'dE'    :   {'3.0':20.0, '4.0':15.0, '-4.0':15.0},
                                    'Gg'    :   {'3.0':80.0, '4.0':55.0, '-4.0':55.0},
                                    'gn2'    :   {'3.0':50.0, '4.0':10.0, '-4.0':10.0}  })


# %% Sample a resonance ladder
resonance_ladder = Ta_pair.sample_resonance_ladder(energy_grid, spin_groups, average_parameters)

#%% ### Calculate a theoretical cross section from the resonance ladder

# SLWB with syndat
xs_tot, xs_scat, xs_cap = syndat.scattering_theory.SLBW(energy_grid, Ta_pair, resonance_ladder)

# convert to transmisison and put in an appropriate dataframe
n = 0.067166 # atoms per barn or atoms/(1e-12*cm^2)
trans = np.exp(-n*xs_tot)
theoretical_df = pd.DataFrame({'E':energy_grid, 'theo_trans':trans})

# %%
options = { 'Perform Experiment':True,
            'Add Noise': True}

new_parameters = {'trigs': {'val':18476117,    'unc'   :   0}}

exp = syndat.experiment(theoretical_df, 
                        options=options, 
                        experiment_parameters=new_parameters)
                        

# ### Output in SAMMY format

# %%
syndat.write_samdat(exp.trans,"./test_sammy.dat")
syndat.sammy_interface.write_sampar(resonance_ladder, Ta_pair, False, "./test_sammy.par")


