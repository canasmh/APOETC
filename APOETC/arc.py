#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Mon Oct  7 14:25:10 2019

@authors: Manuel H. Canas
"""

import numpy as np
from os.path import dirname, abspath
from astropy.io import ascii
from scipy.interpolate import InterpolatedUnivariateSpline

"""[Summary]

:param [ParamName]: [ParamDescription], defaults to [DefaultParamVal]
:type [ParamName]: [ParamType](, optional)
...
:raises [ErrorType]: [ErrorDescription]
...
:return: [ReturnDescription]
:rtype
"""
class Instrument:
    """This object represents the instrument used.

    :param inst_name: This is the name of the instrument used.
    :type inst_name: str

    """

    def __init__(self, inst_name):
        """The constructor for the Instrument class.

        """
        # Path to the directory containing instrument data:
        path_to_dir = dirname(abspath(__file__))+'../data/APO/'+inst_name

        #set attributes
        #--------------

        #Name of instrument
        self.inst_name = inst_name #i.e., 'Artcic'

        #Read out noise of instrument
        self.readout_noise = np.array(ascii.read(path_to_dir + '/readout_noise.dat',
                                                 format='no_header'
                                                 )
                                     )[0][0] #electrons/pix

        #Gain of instrument
        self.gain = np.array(ascii.read(path_to_dir+'/gain.dat',
                                        format='no_header'
                                        )
                             )[0][0] #gain
        #Quantum effiency of instrument
        self.sensitivity = ascii.read(path_to_dir+'/qe.dat')['col2']/100

        #Wavelength of quantum efficiency
        self.wavelength = ascii.read(path_to_dir+'/qe.dat')['col1']*10
        self.plate_scale = np.array(ascii.read(path_to_dir+'/plate_scale.dat')
                                    )[0][0]
        
    def filter(self, bandpass, Johnson = True, SDSS = False):
        """Method that returns the transmission of specified filter.

        :param bandpass: The bandpass of the filter used (i.e., 'U','B','V','R', or 'I').
        :type bandpass: str
        :param Johnson: If true, then the bandpass is referring to the Johnson-Cousin filters. Defaults to True
        :type Johnson: bool, optional
        :param SDSS: If true, then the bandpass is referring to the Johnson-Cousin filters. Defaults to False
        :type SDSS: bool, optional
        ...
        :return: The transmission of the filter interpolated over the bandpass. Also sets a filter_range attribute which
        is a tuple containing (lambda_min,lambda_max) of the filter in Angstroms.
        :rtype: Interpolated object
        """
    
        #Define the path to the filter .dat files based on the type of filter 
        #used.
    
        if Johnson:
            path_to_dir = dirname(abspath(__file__)) + '../data/APO/Filter/Johnson/'
        elif SDSS:
            path_to_dir = dirname(abspath(__file__)) + '../data/APO/Filter/SDSS/'
    
        #Get the transmission and wavelength from the .dat file
        filt_data = ascii.read(path_to_dir+bandpass+'.dat')
            
        #Get the transmission and wavelength from the .dat file
    
        filt_wavelength = filt_data[0][:]
        filt_transmission = filt_data[1][:]

        #Interpolate the filter
        filt = InterpolatedUnivariateSpline(filt_wavelength,\
                                            filt_transmission,k=3)

        #Set the range of the filter as an object attribute.
        setattr(Instrument,filter_range,(filt_wavelength[0],filt_wavelength[-1]))

        #Return the interpolated filter.
        return filt
        #
    def interpolate_efficiency(self):
        
        """Method that interpolates the quantum efficiency.

        :return: The efficiency of the instrument interpolated over the appropriate wavelenghts (in Angstroms).
        """
    
        #Interpolate the efficiency
        efficiency = InterpolatedUnivariateSpline(self.wavelength,
                                                  self.sensitivity,
                                                  k=3)

        return efficiency


class Telescope:
    """Object that represents the telescope used.

    :param obs_name: The name of the observatory used, default to 'ARC 3.5m'.
    :type obs_name: str,optional
    :param aperature: The diameter of the telescope used (in meters), default to 3.5.
    :type aperature: float, optional
    """

    def __init__(self,obs_name = 'ARC 3.5m',aperature = 3.5):
        """The constructor of the Telescope class.

        """

        #Name of observatory
        self.obs_name = obs_name

        #Aperature area
        self.area = np.pi*((aperature*100)/2)**2 #pi r^2 #cm

        #Throughput of telescope.
        self.throughput = 0.90 #Random but reasonable number
