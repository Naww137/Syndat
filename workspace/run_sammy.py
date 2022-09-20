
#%%
import os
import shutil



def run_sammy_and_wait(case_directory, case_basename, number_of_cases):
        
    # delete qsub_icase.sh.* files - these files indicate that qsub job has completed
    for isample in range(number_of_cases):
        wildcard_path = os.path.join(case_directory, case_basename, f'{case_basename}_smpl_{isample}/qsub_{isample}.sh.*')
        if os.path.isfile(wildcard_path):
            os.system(f'rm {wildcard_path}')
        else:
            pass
        
    # run sammy with bayes for all files created
    irunsammy = 0
    for isample in range(number_of_cases):
        directory = os.path.join(case_directory, case_basename,f'{case_basename}_smpl_{isample}')
        os.system("ssh -t necluster.ne.utk.edu 'cd "+directory+f" ; qsub qsub_{isample}.sh'")
        irunsammy += 1
        
    # wait on all cases to complete running - looking for qsub_icase.sh.o file
    running_sammy = True
    print(); print('Waiting for sammy to run'); print()
    while running_sammy:
        case_run_bool = []
        for isample in range(number_of_cases):
            directory = os.path.join(case_directory, case_basename, f'{case_basename}_smpl_{isample}')
            
            idone_file = 0
            for file in os.listdir(directory):
                if file.startswith(f'qsub_{isample}.sh.o'):
                    idone_file += 1
                else:
                    pass
                    
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


def run_leftover_sammy(case_directory, case_basename, number_of_cases):
        
    irunsammy = 0
    for isample in range(number_of_cases):
        lst_path = os.path.join(case_directory, case_basename, f'{case_basename}_smpl_{isample}/SAMMY.LST')
        if os.path.isfile(lst_path):
            pass
        else:
            print(f'Re-running samle {isample}')
            directory = os.path.join(case_directory, case_basename,f'{case_basename}_smpl_{isample}')
            os.system("ssh -t necluster.ne.utk.edu 'cd "+directory+f" ; qsub qsub_{isample}.sh'")
            irunsammy += 1
            
            running_sammy = True
            while running_sammy:
                if os.path.isfile(lst_path):
                    shutil.copy(lst_path, os.path.join(case_directory,'synthetic_data',f'SAMMY_smpl_{isample}.LST'))
                    shutil.copy(os.path.join(case_directory,case_basename,f'{case_basename}_smpl_{isample}', 'sammy_syndat.par'), os.path.join(case_directory,'synthetic_data',f'SAMMY_smpl_{isample}.par'))
                    running_sammy=False
                else:
                    pass
        
    return irunsammy


def copy_syndat(case_directory,case_basename,samples):
    if os.path.isdir(os.path.join(case_directory,'synthetic_data')):
        pass
    else:
        os.mkdir(os.path.join(case_directory,'synthetic_data'))
    run_cases = range(samples); icopy = 0
    for i in run_cases:
        shutil.copy(os.path.join(case_directory,case_basename,case_basename+f'_smpl_{i}', 'SAMMY.LST'), os.path.join(case_directory,'synthetic_data',f'SAMMY_smpl_{i}.LST'))
        shutil.copy(os.path.join(case_directory,case_basename,case_basename+f'_smpl_{i}', 'sammy_syndat.par'), os.path.join(case_directory,'synthetic_data',f'SAMMY_smpl_{i}.par'))
        #os.system("scp nwalton1@necluster.ne.utk.edu:/home/nwalton1/my_sammy/slbw_testing_noexp/slbw_1L_noexp_case1/syndat_{i} /Users/noahwalton/research_local/resonance_fitting/synthetic_data")
        icopy += 1
        # ssh -t necluster.ne.utk.edu 'cd /home/nwalton1/my_sammy/slbw_testing/slbw_fitting_case1/ ; /home/nwalton1/my_sammy/SAMMY/sammy/build/install/bin/sammy < slbw_fitting_case1.sh'
    print(); print(f'copied {icopy} synthetic data files'); print()



#%%

run_sammy = False
run_leftover = False
copy_files = True
case_basename = 'syndat'
samples = 10000

if run_sammy:
    print();print('going to run sammy to create synthetic data'); print()
    irunsammy = run_sammy_and_wait(os.getcwd(), case_basename, samples)


if run_leftover:
    run_leftover_sammy(os.getcwd(), case_basename, samples)

if copy_files:
    copy_syndat(os.getcwd(),case_basename,samples)

# %%


