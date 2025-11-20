import numpy as np
np.set_printoptions(legacy='1.25')
from scipy.integrate import trapezoid
import pandas as pd
pd.set_option('display.max_rows', None)
#import matplotlib.pyplot as plt
import os, glob
import re
import csv

class CrossSectionData(object):
    __doc__="""Class to handle neutron-capture cross section data tables from 
    ENDF-B/VIII.1."""

    _ROOT = os.path.abspath(os.path.dirname(__file__))

    def __init__(self):
        from . import get_data
        self.capture_data_path = get_data('data_capture')
        capture_data_list = [x for x in glob.glob("{0}/*.csv".format(self.capture_data_path))]
        self.capture_data_dict = {}
        for c in capture_data_list:
            c_file = c.split('data_capture/')[1]
            target = c_file.split('n-capture-')[1].split('.csv')[0]
            self.capture_data_dict.update({target: c_file})
        self.capture_data_dict = dict(sorted(self.capture_data_dict.items()))

    def find_targets(self):
        return [target for (target, value) in self.capture_data_dict.items()]

    def get_MT102(self, target):
        """Retrieve point-wise cross section data for defined target as 
        DataFrame object."""
        self.target = target
        df = None
        for k,v in self.capture_data_dict.items():
            if k==self.target:
                df = pd.read_csv("{0}/{1}".format(self.capture_data_path, v))
        if df is None:
            print("No capture-gamma cross section data for target nucleus: {0}".format(target))
            return
        else:
            return df

    def sigma_ENDF(self, target):
        """Convert ENDF energy and cross section DataFrame to numpy arrays 
        for interpolation and integration."""
        self.target = target

        df = CrossSectionData.get_MT102(self,self.target)

        endf_data = df.to_numpy()
        En = endf_data[:,0]
        sigma = endf_data[:,1]

        return (En, sigma)

class NeutronFlux(CrossSectionData):
    __doc__="""Class to handle experimental neutron-flux spectra from the 
    Budapest Research Reactor (BRR) and Garching (FMR II)."""
   
    _ROOT = os.path.abspath(os.path.dirname(__file__))
   
    def __init__(self):
        from . import get_data
        self.flux_data_path = get_data('data_spectra')
        flux_data_list = [x for x in glob.glob("{0}/*.csv".format(self.flux_data_path))]
        self.flux_data_dict = {}
        for i,f in enumerate(flux_data_list):
            f_file = f.split('data_spectra/')[1]
            self.flux_data_dict.update({i: f_file})
        self.flux_data_dict = dict(sorted(self.flux_data_dict.items()))

    def find_flux(self):
        return [(i,f) for (i, f) in self.flux_data_dict.items()]

    def get_flux_df(self,flux):
        """Retrive experimental neutron-flux spectrum as DataFrame object."""
        self.flux = flux
        df = None
        for k,v in self.flux_data_dict.items():
            if k==self.flux:
                df = pd.read_csv("{0}/{1}".format(self.flux_data_path, v))
        if df is None:
            print("Spectrum not defined for argument:".format(self.flux))
            return
        else:
            return df

    def get_flux(self, flux):
        """Convert experimental flux DataFrame into numpy arrays for interpolation 
        and integration."""
        self.flux = flux

        df = NeutronFlux.get_flux_df(self, self.flux)

        flux_spectrum = df.to_numpy()
        En = flux_spectrum[:,0]
        sigma = flux_spectrum[:,1]

        return (En, sigma)

class ResonanceData(NeutronFlux):
    __doc__="""Class to handle the Breit-Wigner and Reich-Moore resonances of 
    the ENDF-B/VIII.1 library."""

    _ROOT = os.path.abspath(os.path.dirname(__file__))
    
    def __init__(self):
        from . import get_data
        self.res_data_path = get_data('data_resonances')

        # Load Breit-Wigner resonances
        breit_wigner_list = [x for x in glob.glob("{0}/BreitWigner/*.csv".format(self.res_data_path))]
        breit_wigner_dict = {}
        for bw in breit_wigner_list:
            bw_file = bw.split('BreitWigner/')[1]
            target = bw_file.split('n-res-')[1].split('.csv')[0]
            breit_wigner_dict.update({target: [bw_file, 'BreitWigner']})
        #Sort dictionary by key
        bw_sorted_dict = dict(sorted(breit_wigner_dict.items()))

        # Load Reich-Moore resonances
        reich_moore_list = [x for x in glob.glob("{0}/ReichMoore/*.csv".format(self.res_data_path))]
        reich_moore_dict = {}
        for rm in reich_moore_list:
            rm_file = rm.split('ReichMoore/')[1]
            target = rm_file.split('n-res-')[1].split('.csv')[0]
            reich_moore_dict.update({target: [rm_file, 'ReichMoore']})
        #Sort dictionary by key
        rm_sorted_dict = dict(sorted(reich_moore_dict.items()))

        # Concatenate and sort the Breit-Wigner and Reich-Moore dictionaries
        merged_res_dict = {**bw_sorted_dict, **rm_sorted_dict}
        self.res_sorted_dict = dict(sorted(merged_res_dict.items()))

    def find_resonances(self,**kwargs):
        """Find list of all 'target+n' systems with resonance parameters in 
        the ENDF-B/VIII.1 library: includes Breit-Wigner and Reich-Moore.

        Arguments:
        Returns:
        Raises:

        Example:
            Find all resonances:
            >find_resonances()
            Find Breit-Wigner resonances only:
            >find_resonances(res='BW')
            Find Reich-Moore resonances only:
            >find_resonances(res='RM')
"""
        self.kwargs = kwargs

        if self.kwargs == {} or self.kwargs is None:
            # return all targets
            return [target for (target, value) in self.res_sorted_dict.items()]
        else:
            for key in self.kwargs.keys():
                if key == 'res':
                    for res in kwargs.values():
                        if res == 'BW':
                            return [target for (target, value) in self.res_sorted_dict.items() if value[1]=='BreitWigner']
                        elif res == 'RM':
                            return [target for (target, value) in self.res_sorted_dict.items() if value[1]=='ReichMoore']
                        else:
                            print("Unknown keyword argument for resonance parametrizations.")
                            print("Use one of the following methods:")
                            print("`find_resonances()`")
                            print("`find_resonances(res='BW')`")
                            print("`find_resonances(res='RM')`")
                    else:
                        print("Unkown key: use one of the following methods:")
                        print("`find_resonances()`")
                        print("`find_resonances(res='BW')`")
                        print("`find_resonances(res='RM')`")

    def get_res_paras(self, target):
        """Extract resonance parameters for a defined target nucleus and return 
        DataFrame object."""
        self.target = target
        df = None
        PARAS_EXIST = False
        for nucleus, csv_object in self.res_sorted_dict.items():
            if target == nucleus and csv_object[1] == 'BreitWigner':
                PARAS_EXIST = True
                df = pd.read_csv("{0}/BreitWigner/{1}".format(self.res_data_path, csv_object[0]))
            elif target == nucleus and csv_object[1] == 'ReichMoore':
                PARAS_EXIST = True
                df = pd.read_csv("{0}/ReichMoore/{1}".format(self.res_data_path, csv_object[0]))
        if PARAS_EXIST == False:
            print("No resonance parameters available for defined target or target does not exist.")
                
        df_sorted = df.sort_values(by='energy')
        return df_sorted
    

    
