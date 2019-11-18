#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from os.path import dirname, abspath
from astropy.io import ascii
from scipy.interpolate import InterpolatedUnivariateSpline


class Instrument:
    """This object represents the instrument used.

    The ARC 3.5 meter telescope has a 5 different instruments available. For a list of instrument names and
    specifications, check out their `user manual <https://www.apo.nmsu.edu/arc35m/>`_.

    Parameters
    -----------
    inst_name : str
        The name of the isntrument used.

    Attributes
    -----------
    inst_name : str
        The name of the instrument.

    readout_noise : float
        The readout noise of the instrument.

    gain : float
        The gain of the instrument.

    sensitivity : arr
        The quantum efficiency (per angstrom) of the instrument.

    wavelength : arr
        The wavelengths corresponding to the sensitivity attribute.

    plate_scale : float
        The plate scale of the instrument.

    """

    def __init__(self, inst_name):
        """The constructor for the Instrument class.

        """
        # Path to the directory containing instrument data:
        self.filter_range = None
        path_to_dir = dirname(abspath(__file__)) + '/data/APO/' + inst_name

        # set attributes
        # --------------

        # Name of instrument
        self.inst_name = inst_name  # i.e., 'Artcic'

        # Read out noise of instrument
        self.readout_noise = np.array(ascii.read(path_to_dir + '/readout_noise.dat',
                                                 format='no_header'
                                                 )
                                      )[0][0]  # electrons/pix

        # Gain of instrument
        self.gain = np.array(ascii.read(path_to_dir + '/gain.dat',
                                        format='no_header'
                                        )
                             )[0][0]  # gain
        # Quantum effiency of instrument
        self.sensitivity = ascii.read(path_to_dir + '/qe.dat')['col2'] / 100

        # Wavelength of quantum efficiency
        self.wavelength = ascii.read(path_to_dir + '/qe.dat')['col1'] * 10
        self.plate_scale = np.array(ascii.read(path_to_dir + '/plate_scale.dat')
                                    )[0][0]

    def filter(self, filter):
        """Method that returns the transmission of specified filter.

        The Astrophysical Research Consortium 3.5m telescope uses both Johnson-Cousins and SDSS prime filters.
        the `filter` parameter must equal one of the following::

            *'U'        *'gprime'
            *'B'        *'uprime'
            *'V'        *'rprime'
            *'R'        *'zprime'
            *'I'        *'iprime'

        Parameters
        ----------
        filter : str
            The name of the filter used (i.e., 'U','B','gprime', etc.).

        """

        # Define the path to the filter .dat files based on the type of filter
        # used.

        path_to_dir = dirname(abspath(__file__)) + '/data/APO/Filters/'

        # Get the transmission and wavelength from the .dat file
        filt_data = ascii.read(path_to_dir + filter + '.dat')

        # Get the transmission and wavelength from the .dat file

        filt_wavelength = filt_data[0][:]
        filt_transmission = filt_data[1][:]

        # Interpolate the filter
        filt = InterpolatedUnivariateSpline(filt_wavelength,
                                            filt_transmission,
                                            k=3
                                            )

        # Set the range of the filter as an object attribute.
        setattr(self, 'filter_range', (filt_wavelength[0], filt_wavelength[-1]))

        # Return the interpolated filter.
        return filt

    def interpolate_efficiency(self):
        """Method that interpolates the quantum efficiency.

        Return
        ------
        efficiency : Interpolated object
            The efficiency of the instrument interpolated over the appropriate wavelengths (in Angstroms).
        """
        # Interpolate the efficiency
        efficiency = InterpolatedUnivariateSpline(self.wavelength,
                                                  self.sensitivity,
                                                  k=3)

        return efficiency


class Telescope:
    """Object that represents the telescope used.

    For our purposes, this telescope really represents the Astrophysical Research Consortium 3.5m telescope
    but we hope to make it more versatile in the future. The ``throughput`` attribute is defined by the reflectivity
    of aluminum from 2000 A to 100000 A as determined by
    `laserbeamproducts.com <https://laserbeamproducts.wordpress.com/2014/06/19/reflectivity-of-aluminium-uv-visible-and-infrared/>`_.

    Parameters
    -----------
    obs_name : str, optional
        The name of the observatory. Default is 'ARC 3.5m'.
    mirror_diameter : float, optional
        The diameter of the primary mirror in cgs. Default is 350 cm.
    coord : tuple, optional
        The latitude and longitude of the observatory. Default set to APO

    Attributes
    -----------
    obs_name : str
        The name of the observatory. Default is 'ARC 3.5m'.
    area : float
        The light gathering area of the primary mirror.

    throughput : Intepolated object
        The throughput of the telescope which is based on the reflectivity of aluminum.

    latitude : str
        Latitude of the observatory. Defaults to 32.7803.

    longitude : str
        Longitude of the observatory. Defaults to -105.8203
    """

    def __init__(self,
                 obs_name='ARC 3.5m',
                 mirror_diameter=350,
                 coord=['32.7803', '-105.8203']
                 ):

        """The constructor of the Telescope class.

        """

        # Set the path to the aluminum data file.

        path_to_file = dirname(abspath(__file__)) + '/data/APO/al_reflectivity.dat'

        # Wavelength of reflectivity of Aluminum is the first column (in Angstroms)
        wavelength = ascii.read(path_to_file)[0]
        transmission = ascii.read(path_to_file)[1]

        # Name of observatory
        self.obs_name = obs_name

        # Aperature area
        self.area = np.pi * (mirror_diameter / 2) ** 2  # pi r^2 #cm

        # Throughput of telescope, based on reflectivity of aluminum.
        self.throughput = InterpolatedUnivariateSpline(wavelength, transmission)
        self.throughput_range = [wavelength[0], wavelength[-1]]

        # Latitude and longitude of the observatory
        self.latitude = coord[0]
        self.longitude = coord[1]
