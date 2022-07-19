#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 12:30:52 2022

@author: nwx
"""

import numpy as np
import os
from pathlib import Path
import syndat
import pandas as pd

# =============================================================================
# 
# =============================================================================
def readlst(filepath):
    df = pd.read_csv(filepath, sep='\s+', names=['E','exp_xs','exp_xs_unc','theo_xs','theo_xs_bayes','exp_trans','exp_trans_unc','theo_trans', 'theo_trans_bayes'])
    return df


# =============================================================================
# 
# =============================================================================
def samtools_fmtpar(a, template, filename):
    
    with open(template, 'r') as f:
        template_lines = f.readlines()
    f.close()
    
    with open(filename,'w+') as f:
        for line in template_lines:
            
            if line.startswith('%%%ResParms%%%'):
                for row in a:
                    # for some reason, the first format string has had problems with the width definition,
                    # if using a different sized energy range (digits before decimal) this may become an issue
                    f.write(f'{row[0]:0<11.4f} {row[1]:0<10f} {row[2]:0<10f} {row[3]:0<10f} {row[4]:0<10f} ')
                    f.write(f'{row[5]:1g} {row[6]:1g} {row[7]:1g} {row[8]:1g} {row[9]:1g} {row[10]:1g}\n')
# =============================================================================
#                     f.write('{:0<11.5f} {:0<10.5f} {:0<10.5f} '.format(row[0],row[1],row[2]))
#                     f.write('{:0<10.5f} {:0<10.5f} {:1d} '.format(row[3],row[4],int(row[5])))
#                     f.write('{:1d} {:1d} {:1d} '.format(int(row[6]),int(row[7]),int(row[8])))
#                     f.write('{:1d} {:1d}\n'.format(int(row[9]),int(row[10])))
# =============================================================================
            else:        
                f.write(line)
        f.close()
        
    return
        
      
# =============================================================================
# 
# =============================================================================
def write_estruct_file(Energies, filename):
    print("WARNING: if 'twenty' is not specified in sammy.inp, the data file format will change.\nSee 'sammy_interface.write_estruct_file'")
    with open(filename,'w') as f:
        for ept in Energies:
            f.write(f'{ept:0<19f} {1.0:<19} {1.0:0<7}\n')
        f.close()
    return


# =============================================================================
#         
# =============================================================================
def create_sammypar(Jn_ladders, Jp_ladders, \
                    filename='sammy.par', \
                    par_template=os.path.join(Path(os.path.dirname(__file__)).parents[0],'templates/sammy_template.par') ):
    
    J = Jn_ladders + Jp_ladders
    samtools_array = np.empty((0,11))
    for ij, j_df in enumerate(J):
        
        j_array = np.array(j_df)
    
        levels = j_array.shape[0]
        zero_neutron_widths = 5-j_array.shape[1] # accepts 3 neutron widths
        
        # zeros for additional neutron widths plus zeros for all binary "vary parameter options
        j_inp_array = np.concatenate( [j_array, np.zeros((levels, zero_neutron_widths+5)), np.full((levels,1),ij+1)] , axis=1)
        
        samtools_array  = np.concatenate([samtools_array, j_inp_array], axis=0)
    
    samtools_fmtpar(samtools_array, par_template, filename)
    
    return


# =============================================================================
# Could update this to do more to the input file, i.e. input energy range
# =============================================================================
def create_sammyinp(filename='sammy.inp', \
                    template=os.path.join(Path(os.path.dirname(__file__)).parents[0],'templates/sammy_template.inp') ):
    
    with open(template, 'r') as f:
        template_lines = f.readlines()
    f.close()
    
    with open(filename,'w+') as f:
        for line in template_lines:
            f.write(line)
        f.close()
        
    return
    
# =============================================================================
# 
# =============================================================================
def read_SAMNDF_PAR(filename):

    energies = []; spin_group = []; nwidth = []; gwidth = []
    with open(filename,'r') as f:   
        readlines = f.readlines()
        in_res_dat = False
        for line in readlines:   
            if line.startswith(' '):
                in_res_dat = False  
            if in_res_dat:
                if line.startswith('-'):
                    continue #ignore negative resonances
                else:
                    splitline = line.split()
                    energies.append(float(splitline[0]))
                    gwidth.append(float(splitline[1]))
                    nwidth.append(float(splitline[2]))
                    spin_group.append(float(splitline[-1]))
            if line.startswith('RESONANCE PARAMETERS'):
                in_res_dat = True
                
                
    Gg = np.array(gwidth); Gn = np.array(nwidth); E = np.array(energies); jspin = np.array(spin_group)
    df = pd.DataFrame([E, Gg, Gn, jspin], index=['E','Gg','Gn','jspin']); df = df.transpose()
    
    #avg_widths = df.groupby('jspin', as_index=False)['Gg','Gn'].mean() 
    gb = df.groupby('jspin')    
    list_of_dfs=[gb.get_group(x) for x in gb.groups]
    
    avg_df = pd.DataFrame(index=df['jspin'].unique(),columns=['dE','Gg','Gn'])

    for ij, jdf in enumerate(list_of_dfs):
        avg_df['dE'][ij+1]=jdf['E'].diff().mean()
        avg_df['Gg'][ij+1]=jdf['Gg'].mean()
        avg_df['Gn'][ij+1]=jdf['Gn'].mean()

    return avg_df
    
     

           


