#%%
import os
import shutil
import syndat
import numpy as np
import pandas as pd



def sample_a_ladder(fixed_res_df):

    # constants
    midpoint = np.diff(fixed_res_df.E)/2 + min(fixed_res_df.E)
    std = (midpoint-fixed_res_df.E[0])/3
    j3_D_avg = 8.38
    j3_Gn_avg = 41.12  #milli ev
    j3_Gg_avg = 60 # milli-ev
    j2n_Gg_avg = 50
    j2n_Gn_avg = 20.8

    j3_bool = np.random.randint(0,2)
    j2n_bool = np.random.randint(0,2)

    if j3_bool == 1:
        j3_location = np.random.default_rng().normal(midpoint, std)
        # j3_S, j3_P = syndat.scattering_theory.PS_recursive(j3_location, pair.ac, pair.M, pair.m, 0)
        j3_Gg = syndat.sample_widths.sample_RRR_widths(j3_location, j3_Gg_avg/10000, 10000)
        j3_Gn = syndat.sample_widths.sample_RRR_widths(j3_location, j3_Gn_avg/7, 7)
        j3_df = pd.DataFrame({'E':j3_location, 'Gg':[j3_Gg], 'Gn':[j3_Gn], 'jspin':[1.0]}) # pd.concat([fixed_res_df,pd.DataFrame({'E':j3_location, 'Gg':[j3_Gg], 'Gn':[j3_Gn], 'jspin':1.0})])
        new_res_df = pd.concat([fixed_res_df,j3_df])
        if j2n_bool == 1:
            j2n_location = syndat.sample_levels.sample_RRR_levels([min(fixed_res_df.E), max(fixed_res_df.E)], 8.38)[0]
            j2n_Gg = syndat.sample_widths.sample_RRR_widths(j2n_location, j2n_Gg_avg/10000, 10000)
            j2n_Gn = syndat.sample_widths.sample_RRR_widths(j2n_location, j2n_Gn_avg/3, 3)
            j2n_df = pd.DataFrame({'E':j2n_location, 'Gg':j2n_Gg, 'Gn':j2n_Gn, 'jspin':np.array([3.0]*len(j2n_location))})
            new_res_df = pd.concat([new_res_df, j2n_df])
    else:
        new_res_df = fixed_res_df

    new_res_array = syndat.sammy_interface.create_samtools_array_from_DF(new_res_df, False)

    return new_res_df, new_res_array

def write_qsub_shell_script(isample, sample_directory):
    with open(os.path.join(sample_directory,f'qsub_{isample}.sh'), 'w') as f:
        f.write("""#!/bin/bash

#PBS -V
#PBS -l nodes=1:ppn=1
#PBS -q fill

cd ${PBS_O_WORKDIR}

/home/nwalton1/my_sammy/SAMMY/sammy/build/install/bin/sammy < piped_sammy_commands.sh""")


#%%

avg, fixed_res_df = syndat.read_sammy_par('/Users/noahwalton/research_local/resonance_fitting/synthetic_data/Ta181/cole_NN/SAMNDF_FixedRes.PAR')
energy = np.linspace(15,45,200)

#%%

samples = 2
case_basename = 'syndat'
run_sammy = False

if os.path.isdir(os.path.join(os.getcwd(), f'{case_basename}')):
    pass
else:
    os.mkdir(os.path.join(os.getcwd(), f'{case_basename}'))

for i in range(samples):

    sample_df, sample_array = sample_a_ladder(fixed_res_df)

    sample_name = case_basename + '_smpl_' + str(i)
    sample_directory = os.path.join(os.getcwd(), f'{case_basename}/{sample_name}')
    
    if os.path.isdir(sample_directory):
        pass
    else:
        os.mkdir(sample_directory)
    
    sammy_inp_filename = 'sammy_syndat.inp'
    sammy_par_filename = 'sammy_syndat.par'
    estruct_filename = 'estruct'
    piped_commands_filename = 'piped_sammy_commands.sh'
    
    # perform submethod actions and return sub-method warnings
    syndat.sammy_interface.write_estruct_file(energy, os.path.join(sample_directory,estruct_filename))
    syndat.sammy_interface.create_sammyinp(filename=os.path.join(sample_directory,sammy_inp_filename), template='/Users/noahwalton/Documents/GitHub/nuc_syndat/templates/sammy_template_Ta181.inp')
    syndat.sammy_interface.samtools_fmtpar(sample_array, os.path.join(sample_directory,sammy_par_filename))
    write_qsub_shell_script(i, sample_directory)
    
    # create piped_sammy_commands.sh
    with open(os.path.join(sample_directory, piped_commands_filename) , 'w') as pipefile:
        line1 = sammy_inp_filename
        line2 = sammy_par_filename
        line3 = estruct_filename
        pipefile.write(line1 + '\n' + line2 + '\n' + line3 + '\n\n')

if run_sammy:
    print();print('going to run sammy to create synthetic data'); print()
    irunsammy = syndat.MMDA.run_sammy_and_wait(os.getcwd(), case_basename, samples)
    syndat.MMDA.copy_syndat(os.getcwd(),case_basename,1,samples)
else:
    irunsammy = 0

# %%
