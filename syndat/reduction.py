#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 15:10:07 2022

@author: noahwalton
"""






class reduction:
    
    def __init__(self, generation):
        
        # build reduction parameters from generating parameters
        self.redpar = generation.redpar
    
        sdat = generation.sdat.filter(['E','tof','nc','dnc'])
        sdat.rename(columns={"nc": "c", "dnc": "dc"}, inplace=True)
        self.sdat = sdat
        
        odat = generation.odat.filter(['E','tof','bw','c','dc'])
    
                                           
        