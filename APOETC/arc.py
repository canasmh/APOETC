"""
Created on Mon Oct  7 14:25:10 2019

@authors: Manuel H. Canas
"""

import numpy as np
from os.path import dirname, abspath
from astropy.io import ascii
from scipy.interpolate import InterpolatedUnivariateSpline
from synphot.models import BlackBody1D
from synphot import units
from astropy import units as u
from synphot import SourceSpectrum

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
    """This CLass represents the instrument used.

    :param inst_name: This is the name of the instrument used
    :type inst_name: str

    """
    
    # Attributes:
    #     inst_name : str
    #         Name of the instrument used.
    #
    #     sensitivity : arr-like
    #         Quantum efficiency of the instrument (per wavelength).
    #
    #     wavelength : arr-like
    #         Wavelength range of instrument in Angstroms (coupled with efficiency).
    #
    #
    #     readout_noise : float
    #         Readout noise of the instrument (electrons/pix).
    #
    #     gain : float
    #         Gain of the instrument.
    #
    #     plate_scale : float
    #         Plate scale of the detector (arcsec/pix).

    def __init__(self, inst_name):
        """The constructor for the Instrument class.
    
        Parameter:
            inst_name : str
                Name of the instrument used.
        """
        # Path to the directory containing instrument data:
        path_to_dir = dirname(abspath(__file__))+'/data/APO/'+inst_name

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
        """Method that returns filter transmission.
    
        Parameters:
            bandpass : str
                The bandpass used (i.e., 'U','B','V','R', or 'I').
                ARC uses Johnson-Cousin and SDSS filters
    
            Johnson : bool (optional)
                If true, then filter bandpass matches the Johnson-Cousin 
                filter. Default is True.
    
            SDSS : bool (optional)
                If true, the filter bandpass matches the SDSS filters. The 
                filters used are SDSS prime filters, but ommit the " ' "symbol
                from the file name. Default is False.
    
        Returns:
            filt_transmission : obj
                The interpolated transmission of the filter.
    
        """
    
        #Define the path to the filter .dat files based on the type of filter 
        #used.
    
        if Johnson:
            path_to_dir = dirname(abspath(__file__)) + '/data/APO/Filter/Johnson/'
        elif SDSS:
            path_to_dir = dirname(abspath(__file__)) + '/data/APO/Filter/SDSS/'
    
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
        
        """
        Method that interpolates the quantum efficiency
        
        Returns:
            x, efficiency : tuple
            X is a 1D array containing 1000 points evenly spread out by
            the filter range. efficiency contains the interpolated
            values for that x range
        """
    
        #Interpolate the efficiency
        efficiency = InterpolatedUnivariateSpline(self.wavelength,
                                                  self.sensitivity,
                                                  k=3)

        return efficiency

class Sky:
    def __init__(self, lunar_phase=0, seeing=1, airmass=1, transmission = 0.90):

        self.lunar_phase = lunar_phase
        self.seeing = seeing
        self.airmass = airmass

        self.transmission()
        self.emission()

    def transmission(self):
        if self.airmass <= 1.25:
            trans_file = 'trans_1.txt'
        elif self.airmass < 1.75 and airmass > 1.25:
            trans_file = 'trans_1_5.txt'
        elif self.airmass >= 1.75 and airmass < 2.25:
            trans_file = 'trans_2.txt'
        elif self.airmass >= 2.25:
            trans_file = 'trans_2_5.txt'

        transmission = np.loadtxt('../data/sky/' + trans_file)
        self.sky_transmission = ius(transmission[:, 0] \
                                                             , transmission[:, 1] \
                                                             )

    def emission(self):
        if lunar_phase < 0.25:
            emission_file = 'moon_00.txt'
        elif lunar_phase >= 0.25 and lunar_phase < 0.75:
            emission_file = 'moon_50.txt'
        elif lunar_phase >= 0.75:
            emission_file = 'moon_100.txt'

        emission = np.loadtxt('../data/sky/' + emission_file)
        self.sky_emission = interpolate.InterpolatedUnivariateSpline(
            emission[:, 0], emission[:, 1])

class Target:
    """Object representing the target star.

    Attributes:
        magnitude : float
            The magnitude of the star you wish to observe.

        magnitude_system : str
            The magnitude used in the above attribute.

        filter_range : tuple
            The band pass of the filter (xmin,xmax).

        SED : obj
            If specified, this will contain the interpolated spectral energy distribution
            of the target star.

        temp : float
            The temperature of the star. This is used only if you wish to
            use Plank's law to obtain the SED.

    """

    def __init__(
            self,
            magnitude,
            magsystem,
            filter_range,
            SED=None,
            temp=5778
    ):
        """Creates an instance of the Target class.

        Parameters:
            magnitude : float
                The magnitude of the star you wish to observe.

            magsystem : obj
                The magnitude system used in the above paramer.

            filter_range : tuple
                The range of the filter used to observe this target.

            SED : obj (optional)
                The spectral energy distribution of the target star.
                Default is None.

            temp : float (optional)
                The temperature of the target star. This is used to create
                a black body spectrum of the star.
        """

        #Use the specified magnitude system.
        if magsystem.lower() == 'vegamag':
            magnitude_system = units.VEGAMAG
        elif magsystem.lower() == 'stmag':
            magnitude_system = u.STmag
        elif magsystem.lower() == 'abnu':
            magnitude_system = u.ABmag


        self.magnitude = magnitude
        self.magnitude_system = magnitude_system
        self.SED = SED
        self.temp = temp
        self.filter_range = filter_range

    def convert_to_flux(self):
        """Convert magnitude of target star to flux.

        Returns:
            The wavelength flux of the target in cgs units.
        """

        #Get the spectrum of Vega
        vega = SourceSpectrum.from_vega()

        #convert to flux using units.convert_flux
        flux = units.convert_flux(
            self.filter_range,
            self.magnitude*magntidue_system,
            units.FLAM,
            vegaspec=vega
        )

        #return the flux of the star
        return flux

    def blackbody_lambda(self):
        """Calculates the spectrum of a blackbody from temperature temp.

        Returns:
            The wavelength flux of the target as determined by a blackbody
        """

        #Get a black body spectrum of temperature temp
        sp = SourceSpectrum(BlackBody1D, temperature=self.temp * u.K)
        sp_new = sp / np.mean(sp(self.range * u.AA,
                                 flux_unit=units.FLAM
                                 ) / self.convert_to_flux()
                              )
        x = sp_new(range(1000, 30000) * u.AA, flux_unit=units.FLAM)
        bb_lam = ius(range(10, 30000), x)

        return bb_lam

class Telescope:
    """Object that represents the telescope used.

    Attributes:
        obs_name : str
            Name of the telescope used.

        area : float
            The light collecting area of the telescope.

        throughput : float
            The throughput of the telescope.
    """

    def __init__(self,obs_name = 'ARC 3.5m',aperature = 3.5):
        """Creates an instance of the ARC 3.5m telescope

        Parameters:
            obs_name : str (optional)
                The name of the telescope used.
                Default is ARC 3.5m telescope.

            aperature : float
                The aperature of the telescope (in meters).
                Default is 3.5m

        Returns:
            None
        """

        #Name of observatory
        self.obs_name = obs_name

        #Aperature area
        self.area = np.pi*((aperature*100)/2)**2 #pi r^2 #cm

        #Throughput of telescope.
        self.throughput = 0.90 #Random but reasonable number
