#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 12:18:04 2022

@author: noahwalton
"""

# import os
import numpy as np
# import matplotlib.pyplot as plt
import pandas as pd
import json

from syndat import sample_widths
from syndat import sample_levels
# from syndat import scattering_theory

class particle_pair:
    """
    _summary_

    _extended_summary_

    Methods
    -------
    quant_vec_sum: 
        Calculates the quantum vector sum of two angular momenta.
    map_quantum_numbers:
        Maps the possible quantum numbers for pair.
    sample_all_Jpi:
        Samples a full resonance parameter ladder for each possible spin group.
    """

    def __init__(self, settings:str):
        """
        Initialization of particle pair object for a given reaction.

        The particle_pair class houses information about the incident and target particle for a reaction of interest. 
        The methods for this class include functions to calculate the open channels 

        Parameters
        ----------
        settings : str
            JSON file name for the settings.
        """

        with open(settings, 'r') as file:
            pair_data = json.load(file)

        # assuming boundary condition selected s.t. shift factor is eliminated for s wave but not others!

        # Particle pair data:
        print(pair_data.keys())
        self.I = pair_data['target']['i']
        assert self.I % 0.5 == 0
        
        self.i = pair_data['projectile']['i']
        assert self.i % 0.5 == 0
        
        assert pair_data['l_max'] % 1 == 0
        self.l_max = int(pair_data['l_max'])
        
        self.M = pair_data['target']['m']     # amu
        self.m = pair_data['projectile']['m'] # amu
        if 'ac' in pair_data['target'].keys():
            self.ac = (1.23*self.M**(1/3))+0.8 # fermi or femtometers
        elif pair_data['target']['ac'] < 1e-7:
            print("WARNING: scattering radius seems to be given in m rather than sqrt(barns) a.k.a. cm^-12")
        else:
            self.ac = pair_data['target']['ac'] # 6.7e-15 # m or 6.7 femtometers

        # Mean parameters:
        self.spin_groups = pair_data['spin_groups']
        all_groups = pair_data['mean_parameters'].keys()
        self.average_parameters = pd.DataFrame({param: {sg: pair_data['mean_parameters'][sg][param] for sg in all_groups} \
                                                for param in ('dE', 'Gg', 'gn2')})

        # Fudge:
        if pair_data['resonances']['fudge'] == 'auto':
            try:
                import fudge
            except ModuleNotFoundError:
                self.use_fudge = True
            else:
                self.use_fudge = False
        elif pair_data['resonances']['fudge'] in (True, False):
            self.use_fudge = pair_data['resonances']['fudge']
        else:
            raise ValueError(f'Unknown value for fudge usage. Possible values:\nTrue, False, \'auto\'')

        # Generate resonances:
        ensembles = ['wigner', 'picket fence', 'poisson', 'goe', 'auto']
        if pair_data['resonances']['ensemble'].lower() in ensembles:
            self.ensemble = pair_data['resonances']['ensemble'].lower()
            if self.ensemble == 'auto':
                self.ensemble = ('goe' if self.use_fudge else 'wigner')
            elif self.ensemble != 'wigner' and not self.use_fudge:
                raise ValueError(f'"{self.ensemble}" ensemble cannot be used without fudge installed.')
        else:
            raise ValueError(f'Unknown ensemble for resonance energy generation. Possible ensembles:\n{ensembles}')

        # Constants:
        self.hbar = 6.582119569e-16 # eV-s
        self.c    = 2.99792458e8 # m/s
        self.m_eV = 939.565420e6 # eV/c^2


    def quant_vec_sum(self, a,b):
        """
        Calculates a quantum vector sum.

        This function performs a quantum vector sum, a.k.a. it maps the quantum 
        triangular relationship between two integers or half integers.

        Parameters
        ----------
        a : float or int
            a variable.
        b : float or int
            a variable.

        Returns
        -------
        numpy.ndarray
            Array of all possible quantum values.
        """
        a=abs(a); b=abs(b)
        vec = np.arange(abs(a-b), a+b+1, 1)
        return vec


    def map_quantum_numbers(self, print_out:bool):
        """
        Maps the possible quantum numbers for pair.

        This function maps out the possible quantum spin numbers (Jpi) for a given
        particle pair up to some maximum considered incident waveform (l-wave).

        Parameters
        ----------
        particle_pair : syndat object
            Particle_pair object containing information about the reaction being studied.
        print_out : bool
            User option to print out quantum spin (J) mapping to console.

        Returns
        -------
        Jn : array-like
            List containing possible J, # of contibuting channels, and contibuting 
            waveforms for negative parity. Formatted as (J,#chs,[l-wave, l-wave])
        Jp : array-like
            List containing possible J, # of contibuting channels, and contibuting 
            waveforms for positive parity. Formatted as (J,#chs,[l-wave, l-wave])
        Notes
        -----
        
        Examples
        --------
        >>> from sample_resparm import sample_spin_groups
        >>> sample_spin_groups.map_quantum_numbers(3/2,1/2,2, False)
        ([(1.0, 1, [0.0]), (2.0, 1, [0.0])],
        [(0.0, 1, [1.0]),
        (1.0, 2, [1.0, 1.0]),
        (2.0, 2, [1.0, 1.0]),
        (3.0, 1, [1.0])])
        """
        
        # define object atributes
        I = self.I
        i = self.i
        l_wave_max = self.l_max

        # now perform calculations
        # Jn = []; Jp = []; 
        Jall = []
        S = self.quant_vec_sum(I,i)
        L = range(l_wave_max+1)

        i_parity = (-1 if i<0 else 1)
        I_parity = (-1 if I<0 else 1)
        S_parity = i_parity*I_parity

        possible_Jpi = {}
        J_negative = []; J_positive = []
        J_all = []
        for i_l, l in enumerate(L):
            this_l = {}
            
            l_parity = (-1)**l
            J_parity = S_parity*l_parity
            
            for i_s, s in enumerate(S):
                js = self.quant_vec_sum(s,l)
                this_l[f's={s}'] = js
                for j in js:
                    if J_parity == 1:
                        # J_positive.append([l,s,j])
                        J_all.append([l,s,j])
                    if J_parity == -1:
                        # J_negative.append([l,s,j])
                        J_all.append([l,s,-j])
                
            possible_Jpi[f'l={l}'] = this_l
                
        if len(J_all) > 0:
            J_total = np.array(J_all)
            J_unique = np.unique(J_total[:,2])

            for j in J_unique:
                entrance_channels = np.count_nonzero(J_total[:,2] == j)
                
                ls = []; ss = []
                for i, jtot in enumerate(J_total[:,2]):
                    if jtot == j:
                        ls.append(J_total[i,0])
                        ss.append(J_total[i,1])
                        
                Jall.append((j,entrance_channels,ls))
            
            
        if print_out:
            print()
            print('The following arrays describe all possible spin groups for a each parity.\n\
    The data is given as a tuple where the first value is the integer \n\
    or half integer total quantum spin J and the second value is the \n\
    number of entrance channels for that spin group. \n\
    * See the dictionary "possible_Jpi" for a nested packing structure.')
        
            print()
            print('Spin group data for all parity\n(Jpi, #Chs, l-waves)')
            for each in Jall:
                print(each)

        # define new attributes for particle_pair object
        # self.Jn = Jn
        # self.Jp = Jp
        self.J = Jall # Jn + Jp

        return


    def sample_resonance_ladder(self, Erange):
        """
        ...
        """

        resonance_ladder = pd.DataFrame()
        if self.use_fudge:
            import sys
            if '/Users/colefritsch/opt/anaconda3/lib/python3.8/site-packages/fodge' not in sys.path:
                sys.path.append('/Users/colefritsch/opt/anaconda3/lib/python3.8/site-packages/fodge')
            from brownies.BNL.restools.resonance_generator import getFakeResonanceSet
            import xData.XYs1d as XYs1dModule

            J_ID = 0
            for ij, j in enumerate(self.spin_groups):
                J_ID += 1

                EB = (np.min(Erange), np.max(Erange))
                dE = self.average_parameters.dE[f'{j[0]}']
                N = round((EB[1]-EB[0])/dE) + 10

                # Mean capture width:
                aveWidths = {}
                Gg = self.average_parameters.Gg[f'{j[0]}']
                # raise TypeError([[EB[0], Gg], [EB[1], Gg]])
                aveWidths['captureWidth'] = XYs1dModule.XYs1d(
                                                data=[[EB[0], Gg], [EB[1], Gg]],
                                                axes=XYs1dModule.XYs1d.defaultAxes(
                                                    labelsUnits={
                                                        XYs1dModule.yAxisIndex: ('width', 'eV'),
                                                        XYs1dModule.xAxisIndex: ('excitation_energy', 'eV')}))
                # Mean neutron width:
                gn2 = self.average_parameters.gn2[f'{j[0]}']
                Gn = gn2 # NOTE: THIS MAY NOT BE CORRECT!!!!
                aveWidths['neutronWidth'] = XYs1dModule.XYs1d(
                                                data=[[EB[0], Gn], [EB[1], Gn]],
                                                axes=XYs1dModule.XYs1d.defaultAxes(
                                                    labelsUnits={
                                                        XYs1dModule.yAxisIndex: ('width', 'eV'),
                                                        XYs1dModule.xAxisIndex: ('excitation_energy', 'eV')}))
                # Level density:
                level_density = XYs1dModule.XYs1d([[EB[0], 1/dE], [EB[1], 1/dE]],
                                                axes=XYs1dModule.XYs1d.defaultAxes(
                                                    labelsUnits={
                                                        XYs1dModule.yAxisIndex: ('width', 'eV'),
                                                        XYs1dModule.xAxisIndex: ('excitation_energy', 'eV')}))
                resonances = getFakeResonanceSet(EB[0], dE, N, self.ensemble, L=j[1], J=j[0], levelDensity=level_density, aveWidthFuncs=aveWidths, DOFs={'neutronWidth': 1, 'captureWidth': 0}, domainMin=EB[0], domainMax=EB[1], widthKeys=('neutronWidth', 'captureWidth'))
                # raise TypeError(resonances.data)
                levels = resonances[:,0]
                red_nwidth = resonances[:,3] # NOTE: Check if its is the correct width
                Gwidth = resonances[:,4] # NOTE: Check if its is the correct width

                E_Gn_gnx2 = pd.DataFrame([levels, Gwidth, red_nwidth, [j[0]]*len(levels), [j[1]]*len(levels), [j[2]]*len(levels), [J_ID]*len(levels)], index=['E','Gg', 'gnx2', 'J', 'chs', 'lwave', 'J_ID'])   
                # assert len(np.unique(j[2]))==1, "Code cannot consider different l-waves contributing to a spin group"
                resonance_ladder = pd.concat([resonance_ladder, E_Gn_gnx2.T])

        else:
            J_ID = 0
            for ij, j in enumerate(self.spin_groups):
                J_ID += 1

                # sample resonance levels for each spin group with negative parity:
                [levels, level_spacing] = sample_levels.sample_RRR_levels(Erange, self.average_parameters.dE[f'{j[0]}'], method='GOE')
                
                # a single radiative capture width is sampled w/large DOF because of many 'partial' radiative transitions to ground state
                # must divide average by the 2*DOF in order to maintain proper magnitude
                red_gwidth = sample_widths.sample_RRR_widths(levels, self.average_parameters.Gg[f'{j[0]}']/2000, 1000)
                Gwidth = 2*red_gwidth # Gbar = 2*gbar b/c P~1 for gamma channels

                # sample observable width as sum of multiple single-channel width with the same average (chi2, DOF=channels)
                red_nwidth = sample_widths.sample_RRR_widths(levels, self.average_parameters.gn2[f'{j[0]}']/j[1], j[1])
                E_Gn_gnx2 = pd.DataFrame([levels, Gwidth, red_nwidth, [j[0]]*len(levels), [j[1]]*len(levels), [j[2]]*len(levels), [J_ID]*len(levels)], index=['E','Gg', 'gnx2', 'J', 'chs', 'lwave', 'J_ID'])   
                # assert len(np.unique(j[2]))==1, "Code cannot consider different l-waves contributing to a spin group"
                resonance_ladder = pd.concat([resonance_ladder, E_Gn_gnx2.T])
    
        # resonance_ladder.sort_values(by=['E'])
        resonance_ladder.reset_index(inplace=True, drop=True)
        return resonance_ladder




### legacy code 

    # def sample_all_Jpi(self,  
    #                     Erange, 
    #                     Davg, Ggavg, gnavg,
    #                     save_csv = False, 
    #                     sammy_run_folder = os.getcwd()):
    #     """
    #     Samples a full resonance parameter ladder for each possible spin group.

    #     This function samples resonance parameters (Energy and widths) for each 
    #     possible spin group (Jpi) of a given particle pair. The results can be 
    #     printed to the console and/or saved to a csv. 

    #     Parameters
    #     ----------
    #     self : syndat object
    #         Particle pair object.
    #     Erange : array-like
    #         Array of resolve resonance range energy, only requires min/max.
    #     Davg : array-like
    #         Nested list of average level spacing for each spin group number. First 
    #         list is for negative parity (J-) second is for positive parity (J+).
    #     Ggavg : array-like
    #         Nested list of average widths for each spin group number. First 
    #         list is for negative parity (J-) second is for positive parity (J+).
    #     gnavg : float
    #         Nested list of average reduced amplitudes (gn_squared) for each spin group number. First 
    #         list is for negative parity (J-) second is for positive parity (J+).
    #     print_out : bool
    #         User option to print out quantum spin (J) mapping to console.
    #     save_csv : bool
    #         User option to save resonance ladders to csv.
    #     sammy_run_folder : str
    #         Folder in which the csv(s) containing resparm ladders will be saved.

    #     Notes
    #     -----
    #     Unsure of the average capture width for Gg sampling.
        
    #     Returns
    #     -------
    #     Jn_df : DataFrame
    #         Pandas DataFrame conatining a resonance parameter ladder for each 
    #         quantum spin group with negative parity (all J-). The column E gives the energy of the level,
    #         the column Gn gives the width of the agregate capture channel, and the following columns give
    #         reduced width amplitudes for particle channels (gn^2), with the headers indicating the waveform (l-wave).
    #     Jp_df : DataFrame
    #         Pandas DataFrame conatining a resonance parameter ladder for each 
    #         quantum spin group with positive parity (all J+). The column E gives the energy of the level,
    #         the column Gn gives the width of the agregate capture channel, and the following columns give
    #         reduced width amplitudes for particle channels (gn^2), with the headers indicating the waveform (l-wave).
    #     """
        
    #     # ensure enough average parameter values were given
    #     Jn_avg_length = [len(Davg[0]), len(Ggavg[0]), len(gnavg[0])]
    #     Jp_avg_length = [len(Davg[1]), len(Ggavg[1]), len(gnavg[1])]
    #     if any(each != len(self.Jn) for each in Jn_avg_length):
    #         raise ValueError("Not enough avarage parameters given for negative parity spin groups")
    #     if any(each != len(self.Jp) for each in Jp_avg_length):
    #         raise ValueError("Not enough avarage parameters given for positive parity spin groups")
            
    # # =============================================================================
    # #     negative parity Js
    # # =============================================================================
    #     Jn_ = []
    #     if len(Davg[0]) > 0:
    #         for ij, j in enumerate(self.Jn):
                
    #             # sample resonance levels for each spin group with negative parity
    #             [levels, level_spacing] = sample_levels.sample_RRR_levels(Erange, Davg[0][ij])
                
    #             # a single radiative capture width is sampled w/large DOF because of many 'partial' radiative transitions to ground state
    #             # must divide average by the DOF in order to maintain proper magnitude
    #             red_gwidth = sample_widths.sample_RRR_widths(levels, Ggavg[0][ij]/100, 100)
    #             Gwidth = 2*red_gwidth # Gbar = 2*gbar b/c P~1 for gamma channels
                
    #             # reduced width amplitudes are sampled as single channel (PT or chi with 1 DOF) for each contributing channel then summed
    #             # while the sum will follow chi square with DOF=#channels, if you just sample the sum over all channels, you ignore
    #             # differences in the average widths and differences in the penetrability function assosciated with each width
    #             gnx=[]; gn_lwave = []
    #             for ichannel, lwave in enumerate(j[2]):      
    #                 red_nwidth = sample_widths.sample_RRR_widths(levels, gnavg[0][ij], 1)
    #                 gnx.append(red_nwidth); gn_lwave.append(lwave)
    #             gn = pd.DataFrame(gnx, index=gn_lwave)
                
    #             E_Gg = pd.DataFrame([levels, Gwidth], index=['E','Gg'])
    #             E_Gg_gnx = pd.concat([E_Gg,gn], axis=0)
    #             E_Gg_Gnx_vert = E_Gg_gnx.transpose()
                
    #             Jn_.append(E_Gg_Gnx_vert)
                
    #             if save_csv:
    #                 E_Gg_Gnx_vert.to_csv(os.path.join(sammy_run_folder, f'Jn_{j[0]}.csv'))
    #     else:
    #         print("No average level spacing given for negative parity spin groups")
                
    # # =============================================================================
    # #       positive parity Js
    # # =============================================================================
    #     Jp_ = []
    #     if len(Davg[1]) > 0:
    #         for ij, j in enumerate(self.Jp):
                
    #             # sample resonance levels for each spin group with negative parity
    #             [levels, level_spacing] = sample_levels.sample_RRR_levels(Erange, Davg[1][ij])
                
    #             # a single radiative capture width is sampled w/large DOF because of many 'partial' radiative transitions to ground state
    #             red_gwidth = sample_widths.sample_RRR_widths(levels, Ggavg[1][ij], 100)
    #             Gwidth = 2*red_gwidth # Gbar = 2*gbar b/c P~1 for gamma channels
                
    #             # reduced width amplitudes are sampled as single channel (PT or chi with 1 DOF) for each contributing channel then summed
    #             # while the sum will follow chi square with DOF=#channels, if you just sample the sum over all channels, you ignore
    #             # differences in the average widths and differences in the penetrability function assosciated with each width
    #             gnx=[]; gn_lwave = []
    #             for ichannel, lwave in enumerate(j[2]):      
    #                 red_nwidth = sample_widths.sample_RRR_widths(levels, gnavg[1][ij], 1)
    #                 gnx.append(red_nwidth); gn_lwave.append(lwave)
    #             gn = pd.DataFrame(gnx, index=gn_lwave)
                
    #             E_Gg = pd.DataFrame([levels, Gwidth], index=['E','Gg'])
    #             E_Gg_gnx = pd.concat([E_Gg,gn], axis=0)
    #             E_Gg_Gnx_vert = E_Gg_gnx.transpose()
                
    #             Jp_.append(E_Gg_Gnx_vert)
                
    #             if save_csv:
    #                 E_Gg_Gnx_vert.to_csv(os.path.join(sammy_run_folder, f'Jp_{j[0]}.csv'))
    #     else:
    #         print("No average level spacing given for positive parity spin groups")
                
        
    #     # =============================================================================
    #     #       redefine object attributes
    #     # ============================================================================= 
    #     self.Jn_resonances = Jn_
    #     self.Jp_resonances = Jp_

    











# if __name__ == '__main__':
#     import sys
#     sys.path.insert(0, '/Users/colefritsch/ENCORE/syndat/nuc_syndat/fudge')
#     from brownies.BNL.restools.resonance_generator import getFakeResonanceSet
#     import xData.XYs1d as XYs1dModule

#     l = 1; j = 0.5
#     N = 100
#     EB = (10, 10000)

#     dE = 1.0
#     Gg = 1.0
#     Gn = 1.0

#     # Mean capture width:
#     aveWidths = {}
#     aveWidths['captureWidth'] = XYs1dModule.XYs1d(
#                                     data=[[EB[0], Gg], [EB[1], Gg]],
#                                     axes=XYs1dModule.XYs1d.defaultAxes(
#                                         labelsUnits={
#                                             XYs1dModule.yAxisIndex: ('width', 'eV'),
#                                             XYs1dModule.xAxisIndex: ('excitation_energy', 'eV')}))
#     # Mean neutron width:
#     aveWidths['neutronWidth'] = XYs1dModule.XYs1d(
#                                     data=[[EB[0], Gn], [EB[1], Gn]],
#                                     axes=XYs1dModule.XYs1d.defaultAxes(
#                                         labelsUnits={
#                                             XYs1dModule.yAxisIndex: ('width', 'eV'),
#                                             XYs1dModule.xAxisIndex: ('excitation_energy', 'eV')}))
#     # Level density:
#     level_density = XYs1dModule.XYs1d([[EB[0], 1/dE], [EB[1], 1/dE]],
#                                     axes=XYs1dModule.XYs1d.defaultAxes(
#                                         labelsUnits={
#                                             XYs1dModule.yAxisIndex: ('width', 'eV'),
#                                             XYs1dModule.xAxisIndex: ('excitation_energy', 'eV')}))
#     resonances = getFakeResonanceSet(EB[0], dE, N, 'goe', L=l, J=j, levelDensity=level_density, aveWidthFuncs=aveWidths, domainMin=EB[0], domainMax=EB[1], widthKeys=('neutronWidth', 'captureWidth'))





