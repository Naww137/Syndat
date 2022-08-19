#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 12:19:52 2022

@author: noahwalton
"""

import numpy as np     
import pandas as pd
import syndat
    


class generation:
    
    def __init__(self, perform_methods, default_exp, add_noise, opendat_filename, sammydat_filename):
        
        if default_exp:
            pardict = {
                'n'         :   {'val'  :   0.12,                'unc'  :   0},
                'trig'      :   {'val'  :   9760770,             'unc'  :   0},
                'tof_dist'  :   {'val'  :   35.185,              'unc'  :   0},
                't0'        :   {'val'  :   3.326,               'unc'  :   0},
                'm'         :   {'val'  :   [1,1,1,1],           'unc'  :   [0.016,0.008,0.018,0.005]},
                'm1'        :   {'val'  :   1,                   'unc'  :   0.016},
                'm2'        :   {'val'  :   1,                   'unc'  :   0.008},
                'm3'        :   {'val'  :   1,                   'unc'  :   0.018},
                'm4'        :   {'val'  :   1,                   'unc'  :   0.005},
                'a'         :   {'val'  :   582.8061256946647,   'unc'  :   np.sqrt(1.14174241e+03)},
                'b'         :   {'val'  :   0.0515158865500879,  'unc'  :   np.sqrt(2.18755273e-05)},
                'ks'        :   {'val'  :   0.563,               'unc'  :   0.563*0.0427},
                'ko'        :   {'val'  :   1.471,               'unc'  :   1.471*0.0379},
                'b0s'       :   {'val'  :   9.9,                 'unc'  :   0.1},
                'b0o'       :   {'val'  :   13.4,                'unc'  :   0.7}    }
        else:
            print("Please define a reduction parameter dictionary specific to your experiment")
            pardict = {}
            
        if perform_methods:
            # workflow
            self.redpar = pd.DataFrame.from_dict(pardict, orient='index')
            # import open data from jesse's experiment
            self.get_odat(opendat_filename)
            # vectorize the background function from jesse's experiment
            self.get_bkg()
            # get sample in data 
            self.get_sdat(sammydat_filename)
            self.generate_raw_data(add_noise)
            
        
    def get_odat(self,filename):
        
        odat = pd.read_csv(filename, sep=',') 
        odat = odat[odat.tof >= self.redpar.val.t0]
        odat.sort_values('tof', axis=0, ascending=True, inplace=True)
        odat.reset_index(drop=True, inplace=True)
        odat['E'] = syndat.exp_effects.t_to_e((odat.tof+self.redpar.val.t0)*1e-6, self.redpar.val.tof_dist, True) 
        odat['bw'] = odat.bin_width*1e-6 
        odat.rename(columns={"counts": "c", "dcounts": "dc"}, inplace=True)
        self.odat = odat
        
        
    def get_bkg(self):
        def f(ti,a,b):
            return a*np.exp(ti*-b)
        self.Bi = f(self.odat.tof,self.redpar.val.a,self.redpar.val.b)

    def get_sdat(self, filename):
        sam = syndat.sammy_interface.readlst(filename)
        T_theo = np.flipud(sam.theo_trans)
        sdat = pd.DataFrame()
        sdat['theo_trans'] = sam.theo_trans #
        sdat['E'] = sam.E
        sdat['tof'] = syndat.exp_effects.e_to_t(sdat.E, self.redpar.val.tof_dist, True)*1e6+self.redpar.val.t0
        sdat.sort_values('tof', axis=0, ascending=True, inplace=True)
        sdat.reset_index(drop=True, inplace=True)
        self.sdat = sdat
    
    
    def sample_turp(self):
        print("Update this function")
        
        # sample/wiggle each true underlying value based on the associated uncertainty 
        
        # if wiggling these values, I will need to re calculate the covariance on a/b background functions
        # if statement to possibly sample open count data based on dcounts
        # if statement to smooth open count data
        
        
    def generate_raw_data(self, add_noise):
        self.sdat, self.odat = syndat.exp_effects.generate_raw_count_data(self.sdat, self.odat, add_noise,
                                                                          self.redpar.val.trig, self.redpar.val.ks,self.redpar.val.ko, 
                                                                          self.Bi, self.redpar.val.b0s, self.redpar.val.b0o, 
                                                                          self.redpar.val.m)
        




    

