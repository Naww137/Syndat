#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 11:26:07 2022

@author: nwx
"""

import os
import syndat


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