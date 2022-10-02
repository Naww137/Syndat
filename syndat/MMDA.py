#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 11:26:07 2022

@author: noahwalton
"""

import os
import syndat
import shutil


def wrapped_sammy_file_creator(number_of_realizations, case_directory, Estruct, \
                               I, i, l_wave_max,  
                               RRR_Erange, 
                               Davg, Gavg, 
                               Gavg_swave, 
                               print_out,
                                   save_csv):
    
    estruct_created = 0; inputs_created = 0; par_created = 0
    
    for irealize in range(1,number_of_realizations+1):
        
        realization_dir = os.path.join(case_directory, f'realization_{irealize}/')
        
        if os.path.isdir(realization_dir):
# in here I could look for existing sammy files and have an option to overwrite or keep
            _ = 0
        else:
            os.mkdir(realization_dir)
            
    #     sample resparms
        Jn_ladders, Jp_ladders = syndat.spin_groups.sample_all_Jpi(I, i, l_wave_max,  
                            RRR_Erange, 
                            Davg, Gavg, 
                            Gavg_swave, 
                            print_out,
                            save_csv, 
                            realization_dir)
    
    #   create necessary sammy files
        syndat.sammy_interface.create_sammyinp(os.path.join(realization_dir,'sammy.inp')); inputs_created+=1
        syndat.sammy_interface.create_sammypar(Jn_ladders, Jp_ladders,os.path.join(realization_dir,'sammy.par')); par_created+=1
    #   could maybe sample a paremter for energy structure, i.e. detector deadtime
        syndat.sammy_interface.write_estruct_file(Estruct, os.path.join(realization_dir,'estruct')); estruct_created+=1
    
    report_string = f'Report for wrapped sammy file creator:\n\
{estruct_created} Energy structure files created\n\
{inputs_created} sammy.inp files created\n\
{par_created} sammy.par files created'
                    
    print();print(report_string); print()
                    
    return report_string




def run_sammy_and_wait(case_directory, case_basename, number_of_cases):
        
    # delete qsub_icase.sh.* files - these files indicate that qsub job has completed
    for isample in range(1,number_of_cases+1):
        wildcard_path = os.path.join(case_directory, case_basename, f'{case_basename}_smpl_{isample}/qsub_{isample}.sh.*')
        os.system(f'rm {wildcard_path}')
        
    # run sammy with bayes for all files created
    irunsammy = 0
    for isample in range(1,number_of_cases+1):
        directory = os.path.join(case_directory, case_basename,f'{case_basename}_smpl_{isample}')
        os.system("ssh -t necluster.ne.utk.edu 'cd "+directory+f" ; qsub qsub_{isample}.sh'")
        irunsammy += 1
        
    # wait on all cases to complete running - looking for qsub_icase.sh.o file
    running_sammy = True
    print(); print('Waiting for sammy to run'); print()
    while running_sammy:
        case_run_bool = []
        for isample in range(1,number_of_cases+1):
            directory = os.path.join(case_directory, case_basename, f'{case_basename}_smpl_{isample}')
            
            idone_file = 0
            for file in os.listdir(directory):
                if file.startswith(f'qsub_{isample}.sh.o'):
                    idone_file += 1
                else:
                    _ = 0
                    
            if idone_file > 0:
                case_run_bool.append(False)
            else:
                case_run_bool.append(True)
                
        if any(case_run_bool):
            continue
        else:
            running_sammy = False
        isamples_still_running = case_run_bool.count(True)
        print(f'Waiting on {isamples_still_running} to complete') #!!! this could be done better - only prints this when all are complete for some reason
        
    return irunsammy


def copy_syndat(case_directory,case_basename,first_case,last_case):
    if os.path.isdir(os.path.join(case_directory, case_basename,'synthetic_data')):
        _ = 0
    else:
        os.mkdir(os.path.join(case_directory, case_basename,'synthetic_data'))
    run_cases = range(1,last_case+1); icopy = 0
    for i in run_cases:
        shutil.copy(os.path.join(case_directory,case_basename,case_basename+f'_smpl_{i}',f'syndat_{i}'), os.path.join(case_directory,case_basename,'synthetic_data'))
        #os.system("scp nwalton1@necluster.ne.utk.edu:/home/nwalton1/my_sammy/slbw_testing_noexp/slbw_1L_noexp_case1/syndat_{i} /Users/noahwalton/research_local/resonance_fitting/synthetic_data")
        icopy += 1
        # ssh -t necluster.ne.utk.edu 'cd /home/nwalton1/my_sammy/slbw_testing/slbw_fitting_case1/ ; /home/nwalton1/my_sammy/SAMMY/sammy/build/install/bin/sammy < slbw_fitting_case1.sh'
    print(); print(f'copied {icopy} synthetic data files'); print()



def write_qsub_shell_script(isample, sample_directory):
    with open(os.path.join(sample_directory,f'qsub_{isample}.sh'), 'w') as f:
        f.write("""#!/bin/bash

#PBS -V
#PBS -l nodes=1:ppn=1
#PBS -q fill

cd ${PBS_O_WORKDIR}

/home/nwalton1/my_sammy/SAMMY/sammy/build/install/bin/sammy < piped_sammy_commands.sh""")



def create_sammy_runfiles(case_basename, samples, energy, ladder_sample_function, inp_template_file, run_sammy,
                            run_directory=os.getcwd()):
    """
    Creates directories and SAMMY runfiles for a number of sample data cases.

    This function will create a directory and all files necessary to run SAMMY for each sample. 
    Currently this funciton is setup to generate theoretical cross sections from resonance parameters
    that can then be put through the syndat methodology to generate experimental noise.

    The directory strucure is:
    - cwd
        - case_basename
            - case_basename_smpl_#
                - sammy runfiles for sample #

    A function for gennerating a resonance ladder must be supplied allowing the user to generate very problem specific resonance ladders.

    Parameters
    ----------
    case_basename : string
        Name of the main directory for which this set of synthetic data will live. This folder is created within the directory that this script is run from.
    samples : int
        Number of sample cases
    energy : array-like
        Energy grid for the calculation
    ladder_sample_function : function
        Function that when called samples a resonance ladder and outputs (dataframe, samtools array).
    inp_template_file : string
        Full path to the template sammy.inp file
    run_sammy : bool
        Boolean option to run sammy or not.
    """
    if os.path.isdir(os.path.join(run_directory, f'{case_basename}')):
        pass
    else:
        os.mkdir(os.path.join(run_directory, f'{case_basename}'))

    for i in range(samples):

        sample_df, sample_array = ladder_sample_function()

        sample_name = case_basename + '_smpl_' + str(i)
        sample_directory = os.path.join(run_directory, f'{case_basename}/{sample_name}')
        
        if os.path.isdir(sample_directory):
            pass
        else:
            os.mkdir(sample_directory)
        
        sammy_inp_filename = 'sammy_syndat.inp'
        sammy_par_filename = 'sammy_syndat.par'
        estruct_filename = 'estruct'
        piped_commands_filename = 'piped_sammy_commands.sh'
        
        # write necessary sammy runfiles
        syndat.sammy_interface.write_estruct_file(energy, os.path.join(sample_directory,estruct_filename))
        syndat.sammy_interface.create_sammyinp(filename=os.path.join(sample_directory,sammy_inp_filename), template=inp_template_file)
        syndat.sammy_interface.samtools_fmtpar(sample_array, os.path.join(sample_directory,sammy_par_filename))
        
        # write qsub shell script and piped sammy input shell script
        write_qsub_shell_script(i, sample_directory)
        with open(os.path.join(sample_directory, piped_commands_filename) , 'w') as pipefile:
            line1 = sammy_inp_filename
            line2 = sammy_par_filename
            line3 = estruct_filename
            pipefile.write(line1 + '\n' + line2 + '\n' + line3 + '\n\n')

    if run_sammy:
        print();print('going to run sammy to create synthetic data'); print()
        irunsammy = syndat.MMDA.run_sammy_and_wait(run_directory, case_basename, samples)
        syndat.MMDA.copy_syndat(run_directory,case_basename,1,samples)
    else:
        irunsammy = 0