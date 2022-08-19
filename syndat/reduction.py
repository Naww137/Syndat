#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 15:10:07 2022

@author: noahwalton
"""


import numpy as np     
import pandas as pd
import syndat
    



class reduction:
    
    def __init__(self, perform_methods, generation):
        
        # build reduction parameters from generating parameters
        self.redpar = generation.redpar
        self.Bi = generation.Bi
        self.odat = generation.odat.filter(['E','tof','bw','c','dc', 'cps','dcps'])
        
        # rename noisey counts from generation as counts for reduction
        sdat = generation.sdat.filter(['E','tof','bw','c','dc'])
        #sdat.rename(columns={"nc": "c", "dnc": "dc"}, inplace=True)
        self.sdat = sdat
        
        # create transmission object
        self.trans = pd.DataFrame()
        self.trans['tof'] = self.sdat.tof
        self.trans['E'] = self.sdat.E
        self.trans['theo'] = generation.sdat.theo_trans
        
        if perform_methods:
            self.get_cps()
            self.reduce()
    
    def get_cps(self):
        self.sdat['cps'], self.sdat['dcps'] = syndat.exp_effects.cts_to_ctr(self.sdat.c, self.sdat.dc, self.odat.bw, self.redpar.val.trig)
        
    def reduce(self):
        sys_unc = self.redpar.unc[['a','b','ks','ko','b0s','b0o','m1','m2','m3','m4']].astype(float)
        self.trans['expT'], self.trans['dT'], self.CovT = syndat.exp_effects.reduce_raw_count_data(self.sdat.tof, 
                                                                self.sdat.c, self.odat.c, self.sdat.dc, self.odat.dc,
                                                                self.odat.bw, self.redpar.val.trig, self.redpar.val.a,self.redpar.val.b, 
                                                                self.redpar.val.ks, self.redpar.val.ko, self.Bi, self.redpar.val.b0s,
                                                                self.redpar.val.b0o, self.redpar.val.m, sys_unc)
        
        