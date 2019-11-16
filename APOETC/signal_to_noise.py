#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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