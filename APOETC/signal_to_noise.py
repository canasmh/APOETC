#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 10 23:32:35 2019

@author: Manuel H. Canas
"""


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