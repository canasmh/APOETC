#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import ephem
from os.path import dirname, abspath
import tools
from scipy.interpolate import InterpolatedUnivariateSpline


class Moon:

    """Function that calculates the position and emission of the moon.

    Parameters
    ----------
    lunar_phase : float
        Number from 0 to 1 where 0 is a new moon and 1 is the full moon.

    utc : str
        The Coordinated Universal Time at which you wish to do your observation. Follow the format 'YYY/MM/DD HH:MM'

    lat : str, optional
        The geocentric latitude at which you are doing your observation. Defaults to '32.7803'

    lon : str, optional
        The geocentric longitude at which you are doing your observation. Defaults to '-105.8203'

    Attributes
    ----------
    lunar_phase : float
        Number from 0 to 1 where 0 is a new moon and 1 is the full moon.

    utc : str
        The Coordinated Universal Time at which you wish to do your observation. Follow the format 'YYY/MM/DD HH:MM'

    lat : str, optional
        The geocentric latitude at which you are doing your observation. Defaults to '32.7803'

    lon : str, optional
        The geocentric longitude at which you are doing your observation. Defaults to '-105.8203'

    """

    def __init__(self,
                 lunar_phase,
                 utc,
                 lat='32.7803',
                 lon='-105.8203'
                 ):

        """Constructor method for the Moon class.

        """

        self.lunar_phase = lunar_phase
        self.UTC = utc
        self.lat = lat
        self.lon = lon

    def emission(self):
        """Method to determine the interpolated emission of the moon

        This method takes the lunar phase and uses a data file to determine the emission over a range of wavelengths.

        Returns
        -------
        emission : Interpolated Object
            The emission of the moon.

        """

        path_to_dir = dirname(abspath(__file__))+'/data/Sky/'

        if lunar_phase < 0.25:
            emission_file = 'moon_00.txt'
        elif 0.25 <= lunar_phase < 0.75:
            emission_file = 'moon_50.txt'
        elif lunar_phase >= 0.75:
            emission_file = 'moon_100.txt'

        emission = np.loadtxt(path_to_dir + emission_file)
        emission = interpolate.InterpolatedUnivariateSpline(emission[:, 0], emission[:, 1])
        setattr(Moon, 'range', [emission[:, 0][0],emission[:, 0][-1]])

        return emission

    def coord(self):
        """Method that returns the horizontal coordinates of the moon.

        This method takes the latitude, longitude, and Coordinated Universal Time to calculate the position of the moon.

        Returns
        --------
            alt, ax : tuple
                A tuple containing the altitude and azimuthal angles of the moon in degrees.

        """

        apo = ephem.Observer()
        apo.lat = self.lat
        apo.lon = self.lon
        apo.date = self.UTC
        moon = ephem.Moon()
        moon.compute(apo)

        return np.degrees(moon.alt), np.degrees(moon.az)


class Target:

    """Object representing the target star.

    This object is only functional for unresolved or point sources. You might get away with using it for planets, but
    I definitely wouldn't recommend using this package if you're observing something like a galaxy.

    Parameters
    -----------
    coord : str
        The equatorial coordinate of the star you wish to observe. Right ascension in hours and declination in degrees.
        'HH:MM:SS.S DD:MM:SS.S'.

    local_sidereal_time : str
        The local sidereal time of the observation in the format 'HH:MM:SS.S'.

    mdt : bool
        Mountain Daylight Time. True if you are observing between March 10 - November 03.

    magnitude : float
        The magnitude of the target star.

    magnitude_system : str
        The magnitude system used to define ``magnitude``.

    range : tuple
        The range of interest. If using a filter, then this should be the filter range. If using a spectrograph, then
        refer to the dispersion of the spectrograph to determine your range.

    temp = float, optional
        The temperature of the target star. This is used to calculate a black body spectrum of the same temperature.

    Attributes
    ----------
    alt : float
        The altitude of the star in degrees.

    az : float
        The azimuthal angle of the star in degrees.

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

    airmass : float
        The airmass of the star. This is calculated using the relationship X = sec(z), where
        the zenith angle is 1 - alt.

    """

    def __init__(
            self,
            coord,
            local_sidereal_time,
            magnitude,
            magnitude_system,
            filter_range,
            temp=5778
            ):

        """Constructor method of the Target class.

        """

        # Split right ascension and declination.

        ra = coord.split(' ')[0]
        dec = coord.split(' ')[1]

        # Use the tools.equatorial_to_horizontal() function to convert from Equatorial coordinates to Horizontal
        # coordinates
        alt, az = tools.equatorial_to_horizontal(dec,ra,local_sidereal_time)

        # Use the specified magnitude system.
        if magsystem.lower() == 'vegamag':
            magnitude_system = units.VEGAMAG
        elif magsystem.lower() == 'stmag':
            magnitude_system = u.STmag
        elif magsystem.lower() == 'abnu':
            magnitude_system = u.ABmag

        # Set the attributes
        self.alt = alt
        self.az = az
        self.magnitude = magnitude
        self.magnitude_system = magnitude_system
        self.temp = temp
        self.filter_range = filter_range
        self.airmass = 1/np.cos(1-np.radians(self.alt))

    def convert_to_flux(self):
        """Convert magnitude of target star to flux.

        The spectrum of this flux will be identical to the spectrum of Vega.

        Returns
        -------
        flux :
            The wavelength flux of the target in cgs units.
        """

        # Get the spectrum of Vega
        vega = SourceSpectrum.from_vega()

        # Convert to flux using units.convert_flux
        flux = units.convert_flux(
            self.filter_range,
            self.magnitude*magntidue_system,
            units.FLAM,
            vegaspec=vega)

        # Return the flux of the star
        return flux

    def blackbody_lambda(self):

        """Calculates the spectrum of a blackbody from temperature temp.

        Returns:
            The wavelength flux of the target as determined by a blackbody
        """

        # Get a black body spectrum of temperature temp
        sp = SourceSpectrum(BlackBody1D, temperature=self.temp * u.K)
        sp_new = sp / np.mean(sp(self.range * u.AA,
                                 flux_unit=units.FLAM
                                 ) / self.convert_to_flux()
                              )
        x = sp_new(range(1000, 30000) * u.AA, flux_unit=units.FLAM)
        bb_lam = ius(range(10, 30000), x)

        return bb_lam

    def transmission(self):
        """Determines the amount of transmission allowed by the sky based on the targets altitude.

         This method takes the airmass of the star and searches through the package data files to find the file that
         most closely resembles that airmass. It returns the transmission of the sky interpolated over the newly defined
         attribute ``Target.transmission_range`` which is a tuple containing the wavelength range in angstroms contained
         in the data file.

         Returns
         -------
            interpolated_transmission : Interpolated Object
                The transmission of the filter interpolated.

         """

        if self.airmass <= 1.25:
            trans_file = 'trans_1.txt'
        elif self.airmass < 1.75 and airmass > 1.25:
            trans_file = 'trans_1_5.txt'
        elif self.airmass >= 1.75 and airmass < 2.25:
            trans_file = 'trans_2.txt'
        elif self.airmass >= 2.25:
            trans_file = 'trans_2_5.txt'

        path_to_file = dirname(abspath(__file__)) +'/data/Sky/' + trans_file

        transmission = np.loadtxt(path_to_file)

        wavelength = transmission[:, 0]*10
        trans = transmission[:, 1]

        interpolated_transmission = InterpolatedUnivariateSpline(wavelength, trans)
        setattr(Target, 'transmission_range', [wavelength[0], wavelength[-1]])

        return interpolated_transmission




#
# class Observation:
#     def __init__(self,  target, sky, instrument, telescope=None):
#
#        # telescope_transm = telescope.transmission
#         self. telescope_area = (175**2)*np.pi
#         self.source = target.F_lambda
#         self.skySED = sky.sky_emission
#         self.skyTransmission = sky.sky_transmission
#
#
#         self.counts(self.source, instrument, 'Source')
#         self.counts(self.skySED, instrument, 'Sky')
#         self.Npix(sky, instrument)
#
#     def Npix(self, sky, instrument):
#         self.Npix = np.pi*((sky.seeing/2)**2)/(instrument.scale**2)
#
#     def counts(self, source, instrument, SourceOrSky):
#         att = dir(instrument)
#         self.detector_qe = instrument.efficiency
#         for row in att:
#             if row.find('filter') > 0:
#
#                 filter_profile = getattr(instrument, row)
#                 integrate_range = getattr(instrument,row.replace('filter', 'range'))
#                 interpolationrange = range(integrate_range[0], integrate_range[1])
#                 h= 6.626*10**(-27) #ergs*s
#                 c=2.9979*10**(18) #A/s
#                 s_integrade = s_integradeInterpolate([source, self.detector_qe, self.skyTransmission, filter_profile], interpolationrange)
#
#                 s_prime = self.telescope_area*(1/(h*c))*s_integrade.integral(integrate_range[0], integrate_range[1])
#                 count_name = row.replace('_filter', '') +'_'+SourceOrSky+'countrate'
#                 setattr(Observation, count_name, s_prime)
#
#
#     def SNfromTime(self, exptime):
#
#        if row.find('filter') > 0:
#         for filter in filterlist:
#             Sprimefilter = sprimefilter
#             BprimeAfilter = bprimeafilter
#             self.rdnoise = telescope.readnoise
#
#         if row.find('filter') > 0.5:
#             for filter in filterlist:
#                 Sprimefilter = sprimefilter
#                 BprimeAfilter = bprimeafilter
#                 self.rdnoise = telescope.readnoise
#                 self.t = exptime
#
#                 filter+"SN" = (Sprimefilter*self.T*self.t)/np.sqrt(Sprimefilter*self.T*self.t + BprimeAfilter*self.T*self.t + self.Npix*self.rdnoise**2)
#
#         if row.find('filter') < 0.5:
#             Sprime = sprime
#             Bprime = bprime
#             self.rdnoise = []
#             self.t = exptime
#
#             for i in inst_range:
#                 self.rednoise.append(inst_rdnoise)
#
#             # PLOT SHIT HERE
#
#
#     def TimefromSN(self, SN):
#         if row.find('filter') > 0.5:
#             for filter in filterlist:
#                 Sprimefilter = sprimefilter
#                 Bprimefilter = bprimefilter
#                 self.rdnoise = RDnoise
#                 SN = signaltonoise
#
#                 t = (1./(2.*Sprimefilter**2))*(SN**2*(Sprimefilter+Bprimefilter)+np.sqrt(SN**4*(Sprimefilter+Bprimefilter)**2+4.*self.Npix*(Sprimefilter*SN*self.rdnoise)**2))
#
#         if row.find('filter') < 0.5:
#             Sprime = sprime
#             Bprime = bprime
#             self.rdnoise = []
#             SN = signaltonoise
#
#             for i in inst_range:
#                 self.rednoise.append(inst_rdnoise)
#
#                 # PLOT SHIT HERE
#
#
# def s_integradeInterpolate(functions, interpolation_range):
#     for i, f in enumerate(functions):
#         if i == 0:
#             x = np.ones(len(interpolation_range))
#         x = f(interpolation_range) * x
#
#     return interpolate.InterpolatedUnivariateSpline(interpolation_range, (x * interpolation_range))
