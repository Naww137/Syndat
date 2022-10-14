

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 12:19:52 2022

@author: noahwalton
"""

import numpy as np     
import pandas as pd
import syndat
    


class experiment:
    
    def __init__(self, theoretical_data,

                    open_data=None, energy_domain = None,

                    experiment_parameters = {} , 
                    options = { 'Perform Experiment':True, 
                                'Add Noise':True, 
                                'Sample TURP':False} , 
                                                                ):
        """
        Instantiates the experiment object

        To synthesize an experiment, you must provide open count data and a theoretical cross section. 
        Options and alternative experimental parameters are optional inputs.

        Parameters
        ----------
        theoretical_data : DataFrame or str
            Theoretical cross section expected to be seen in the laboratory setting. If DataFrame needs columns 'E' and 'theo_trans'. If string, must be filepath to sammy.lst.
        open_data : DataFrame or str, optional
            Open count data. If passing a DataFrame, needs columns 'tof', 'bw', 'c', 'dc'. If passing a string, must be filepath to csv with format from Brown, et al.
            If empty (Default), the open count spectra will be approximated with an exponential function as detailed in Walton, et al.
        energy_domain : array-like, optional
            Energy domain of interest, can be just min/max or entire energy grid. If giving min/max, must be within the domain of the theoretical data.
            If giving energy grid, it must align with the energy grid given with the theoretical data.
            If empty (Default) the entire energy grid from the theoretical data will be used.
        experiment_parameters : dict, optional
            Experimental parameters alternative to default. Default parameters are described in Walton, et al., based on work by Brown, et al.,
            any parameters given here will replace the default parameters before the experiment is synthesized, by default {}
        options : dict, optional
            Keyword options, mostly for debugging, by default { 'Perform Experiment':True, 'Add Noise': True}
        """
        

        ### Gather options
        perform_experiment = options['Perform Experiment']
        add_noise = options['Add Noise']


        ### Default experiment parameter dictionary
        default_exp = {
                        'n'         :   {'val'  :   0.12,                'unc'  :   0},
                        'trigo'     :   {'val'  :   9758727,             'unc'  :   0},
                        'trigs'     :   {'val'  :   18476117,            'unc'  :   0},
                        'tof_dist'  :   {'val'  :   35.185,              'unc'  :   0},
                        't0'        :   {'val'  :   3.326,               'unc'  :   0},
                        'm1'        :   {'val'  :   1,                   'unc'  :   0.016},
                        'm2'        :   {'val'  :   1,                   'unc'  :   0.008},
                        'm3'        :   {'val'  :   1,                   'unc'  :   0.018},
                        'm4'        :   {'val'  :   1,                   'unc'  :   0.005},
                        'a'         :   {'val'  :   582.7768594580712,   'unc'  :   np.sqrt(1.14395753e+03)},
                        'b'         :   {'val'  :   0.05149689096209191, 'unc'  :   np.sqrt(2.19135003e-05)},
                        'ab_cov'    :   {'val'  :   1.42659922e-1,       'unc'  :   None},
                        'ks'        :   {'val'  :   0.563,               'unc'  :   0.02402339737495515},
                        'ko'        :   {'val'  :   1.471,               'unc'  :   0.05576763648617445},
                        'b0s'       :   {'val'  :   9.9,                 'unc'  :   0.1},
                        'b0o'       :   {'val'  :   13.4,                'unc'  :   0.7}    }

        ### redefine experiment parameter dictionary if any new values are given
        pardict = default_exp
        for old_parameter in default_exp:
            if old_parameter in experiment_parameters:
                pardict.update({old_parameter:experiment_parameters[old_parameter]})

        ### set reduction parameter attributes from given options
        self.redpar = pd.DataFrame.from_dict(pardict, orient='index')
        ### sample true underlying resonance parameters
        self.sample_turp(default_exp)

        
    
        ### read in theoretical cross section/transmission data - defines self.sdat
        self.read_theoretical(theoretical_data)



        ### Determine energy grid
        if energy_domain is None:
            self.energy_domain = self.sdat.E
        else:
            if len(energy_domain) != len(self.sdat.E):
                self.sdat = self.sdat[(self.sdat.E>=min(energy_domain))&(self.sdat.E<=max(energy_domain))].reset_index(drop=True)
            if len(energy_domain) == 2:
                self.energy_domain = self.sdat.E
            if np.array_equal(np.array(self.energy_domain), np.array(self.sdat.E)):
                pass
            else:
                raise ValueError("An energy grid was given but it does not line up with that of the theoretical data")



        ### Decide on an open spectra
        if open_data is None:
            self.odat = self.approximate_open_spectra(self.energy_domain)
        else:
            self.read_odat(open_data)
            
        ### sample a realization of the theoretical, true, underlying open count spectra
        self.sample_true_open_spectrum()

        ### Automatically perform experiment
        if perform_experiment:
            # vectorize the background function from jesse's experiment
            self.get_bkg()
            # generate raw count data for sample in given theoretical transmission and assumed true reduction parameters/open count data
            self.generate_sdat(add_noise)
            # reduce the experimental data
            self.reduce()
            

    # ----------------------------------------------------------
    #    Begin Methods
    # ----------------------------------------------------------

    def read_theoretical(self, theoretical_data):
        """
        Reads in a theoretical cross section.

        Parameters
        ----------
        theoretical_data : DataFrame or str
            If DataFrame, must contain clumns 'E' and 'theo_trans'. If str, must be the full path to a sammy.lst file.

        Raises
        ------
        ValueError
            _description_
        ValueError
            _description_
        ValueError
            _description_
        """
        
        # check types
        if isinstance(theoretical_data, pd.DataFrame):
            theo_df = theoretical_data
            if 'E' not in theo_df.columns:
                raise ValueError("Column name 'E' not in theoretical DataFrame passed to experiment class.")
            if 'theo_trans' not in theo_df.columns:
                raise ValueError("Column name 'theo_trans' not in theoretical DataFrame passed to experiment class.")
        elif isinstance(theoretical_data, str):
            theo_df = syndat.sammy_interface.readlst(theoretical_data)
        else:
            raise ValueError("Theoretical data passed to experiment class is neither a DataFrame or path name (string).")
            
        sdat = pd.DataFrame()
        sdat['theo_trans'] = theo_df.theo_trans #
        sdat['E'] = theo_df.E
        sdat['tof'] = syndat.exp_effects.e_to_t(sdat.E, self.redpar.val.tof_dist, True)*1e6+self.redpar.val.t0
        sdat.sort_values('tof', axis=0, ascending=True, inplace=True)
        sdat.reset_index(drop=True, inplace=True)

        self.sdat = sdat


# --------------------------------------------------------------------------------------------------------------------------
    
  
    def read_odat(self,open_data):
        """
        Reads in an open count dataset.

        Parameters
        ----------
        open_data : DataFrame or str
            If DataFrame, must contain clumns 'tof','bw', 'c', and 'dc'. If str, must be the full path to a csv file.

        Raises
        ------
        ValueError
            _description_
        ValueError
            _description_
        ValueError
            _description_
        ValueError
            _description_
        ValueError
            _description_
        """
        
        # check types for what to do 
        if isinstance(open_data, pd.DataFrame):
            if 'tof' not in open_data.columns:
                raise ValueError("Column name 'tof' not in open count DataFrame passed to experiment class.")
            if 'bw' not in open_data.columns:
                raise ValueError("Column name 'bw' not in open count DataFrame passed to experiment class.")
            if 'c' not in open_data.columns:
                raise ValueError("Column name 'c' not in open count DataFrame passed to experiment class.")
            if 'dc' not in open_data.columns:
                raise ValueError("Column name 'dc' not in open count DataFrame passed to experiment class.")

            # calculate energy from tof and experiment parameters
            open_data['E'] = syndat.exp_effects.t_to_e((open_data.tof-self.redpar.val.t0)*1e-6, self.redpar.val.tof_dist, True) 
            odat = open_data
            
        elif isinstance(open_data, str):
            # -------------------------------------------
            # the below code takes open data from Jesse's csv can gets it into the correct format, rather, I want to make thid function take the proper format
            # -------------------------------------------
            odat = pd.read_csv(open_data, sep=',') 
            odat = odat[odat.tof >= self.redpar.val.t0]
            odat.sort_values('tof', axis=0, ascending=True, inplace=True)
            odat.reset_index(drop=True, inplace=True)
            odat['E'] = syndat.exp_effects.t_to_e((odat.tof-self.redpar.val.t0)*1e-6, self.redpar.val.tof_dist, True) 
            odat['bw'] = odat.bin_width*1e-6 
            odat.rename(columns={"counts": "c", "dcounts": "dc"}, inplace=True)
            # -------------------------------------------
        else:
            raise TypeError("Open data passed to experiment class is not of type DataFrame or string (pathname).")


        # filter to energy limits
        if len(odat.E) != len(self.energy_domain):
            # must round to match .LST precision
            odat = odat[(round(odat.E,10)>=min(self.energy_domain))&(round(odat.E,10)<=max(self.energy_domain))].reset_index(drop=True)
        if np.allclose(np.array(odat.E), np.array(self.energy_domain)):
            pass
        else:
            raise ValueError("The open data's energy grid does not align with the defined experiment.energy_domain")

        # Define class attribute
        self.odat = odat
        

# --------------------------------------------------------------------------------------------------------------------------

        
    def get_bkg(self):
        def f(ti,a,b):
            return a*np.exp(ti*-b)
        self.Bi = f(self.odat.tof,self.redpar.val.a,self.redpar.val.b)


# --------------------------------------------------------------------------------------------------------------------------

    def sample_true_open_spectrum(self):
        
        theo_odat = self.odat.copy()
        realization_of_true_cts = syndat.exp_effects.pois_noise(theo_odat.c)
        theo_odat['c'] = realization_of_true_cts
        theo_odat['dc'] = np.sqrt(realization_of_true_cts)

        self.theo_odat = theo_odat

    
    def sample_turp(self, default_exp):
        
        theo_redpar = default_exp
        for par in theo_redpar:
            if par == 'ab_cov':
                continue
            theo_redpar[par]['val'] = np.random.default_rng().normal(theo_redpar[par]['val'], theo_redpar[par]['unc'])
        
        self.theo_redpar = pd.DataFrame.from_dict(theo_redpar, orient='index')

# --------------------------------------------------------------------------------------------------------------------------


    def generate_sdat(self, add_noise):
        """
        Generates a set of noisy, sample in count data from a theoretical cross section via the novel un-reduction method (Walton, et al.).

        Parameters
        ----------
        add_noise : bool
            Whether or not to add noise to the generated sample in data.

        Raises
        ------
        ValueError
            _description_
        """

        if len(self.theo_odat) != len(self.sdat):
            raise ValueError("Experiment open data and sample data are not of the same length, check energy domain")

        monitor_array = [self.theo_redpar.val.m1, self.theo_redpar.val.m2, self.theo_redpar.val.m3, self.theo_redpar.val.m4]

        self.sdat = syndat.exp_effects.generate_raw_count_data(self.sdat, self.theo_odat, add_noise,
                                                                self.theo_redpar.val.trigo, self.theo_redpar.val.trigs, 
                                                                self.theo_redpar.val.ks,self.theo_redpar.val.ko, 
                                                                self.Bi, self.theo_redpar.val.b0s, self.theo_redpar.val.b0o, 
                                                                monitor_array)
        

# --------------------------------------------------------------------------------------------------------------------------


    def reduce(self):
        """
        Reduces the raw count data (sample in/out) to Transmission data and propagates uncertainty.

        """

        # create transmission object
        self.trans = pd.DataFrame()
        self.trans['tof'] = self.sdat.tof
        self.trans['E'] = self.sdat.E
        self.trans['theo_trans'] = self.sdat.theo_trans

        # filter sdat object to be experimental 
        sdat = self.sdat.filter(['E','tof','bw','c','dc', 'theo_cts'])
        self.sdat = sdat

        # get count rates for sample in data
        self.sdat['cps'], self.sdat['dcps'] = syndat.exp_effects.cts_to_ctr(self.sdat.c, self.sdat.dc, self.odat.bw, self.redpar.val.trigs)

        # define systematic uncertainties
        sys_unc = self.redpar.unc[['a','b','ks','ko','b0s','b0o','m1','m2','m3','m4']].astype(float)
        monitor_array = [self.redpar.val.m1, self.redpar.val.m2, self.redpar.val.m3, self.redpar.val.m4]

        self.trans['exp_trans'], self.trans['exp_trans_unc'], self.CovT, self.CovT_stat, self.CovT_sys, rates = syndat.exp_effects.reduce_raw_count_data(self.sdat.tof, 
                                                                                                        self.sdat.c, self.odat.c, self.sdat.dc, self.odat.dc,
                                                                                                        self.odat.bw, self.redpar.val.trigo, self.redpar.val.trigs, self.redpar.val.a,self.redpar.val.b, 
                                                                                                        self.redpar.val.ks, self.redpar.val.ko, self.Bi, self.redpar.val.b0s,
                                                                                                        self.redpar.val.b0o, monitor_array, sys_unc, self.redpar.val.ab_cov)
        

# --------------------------------------------------------------------------------------------------------------------------


    def approximate_open_spectra(self, energy_grid):

        def open_count_rate(tof):
            return (2212.70180199 * np.exp(-3365.55134779 * tof*1e-6) + 23.88486286) 

        tof = syndat.exp_effects.e_to_t(energy_grid,35.185,True)*1e6 # microseconds

        # calculate a tof count rate spectra, convert to counts, add noise 
        cps_open_approx = open_count_rate(tof)
        bin_width = abs(np.append(np.diff(tof), np.diff(tof)[-1])*1e-6)
        cts_open_approx = cps_open_approx*bin_width*self.redpar.val.trigo
        cts_open_measured = syndat.exp_effects.pois_noise(cts_open_approx)
        cps_open_noisy = cts_open_measured/bin_width/self.redpar.val.trigo

        open_dataframe = pd.DataFrame({'tof'    :   tof,
                                        'bw'    :   bin_width,
                                        'c'     :   cts_open_measured,
                                        'dc'    :   np.sqrt(cts_open_measured)})

        open_dataframe['E'] = syndat.exp_effects.t_to_e((open_dataframe.tof-self.redpar.val.t0)*1e-6, self.redpar.val.tof_dist, True) 

        return open_dataframe


    

