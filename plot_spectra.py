#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 14:21:00 2018

@author: charlee
"""

import numpy as np
from matplotlib import pyplot as plt

#Constants
h = 6.626e-27 #planck const in erg/s
c = 3e10 #speed of light cm/s 
kT = 0.35 #keV - 0.288 @ infinity, 0.35 @ surface
dist = 6.1713552e20 #200 pc in cm
Req = 1.2e6 #equatorial radius in cm 
M_over_R = 0.34466/2 #GM/Rc^2
redshift = 1/np.sqrt(1 - 2*M_over_R) #dimensionless 
solid_angle = np.pi*redshift**2*(Req/dist)**2

kev_to_erg = 1.6e-9
E_bin = 0.1*kev_to_erg

#Load spectral data from flux_master code
data = np.loadtxt('spectra.txt', skiprows = 1)

E_obs = data[:,0] #kev
Flux = data[:,1] #photons/cm^2
BB_only = data[:,2] #erg/cm^2

E_em = E_obs*redshift #observed photon energy in keV
E_obs_erg = E_obs*kev_to_erg #observed photon energy in ergs
E_em_erg = E_obs*kev_to_erg #emitted photon energy in ergs

#Compute a BB Spectrum at given energy kT
BB_em_1 = (2*(E_em_erg)**3/(c**2*h**2))
BB_obs_1 = (2*(redshift*E_obs_erg)**3/(c**2*h**2))
BB_em_2 = np.e**((E_em)/kT) - 1
BB_obs_2 = np.e**((E_obs*redshift)/kT) - 1

BB_em = (BB_em_1/BB_em_2)
BB_obs = BB_obs_1/BB_obs_2
Flux_obs = (BB_obs*solid_angle)/(redshift**3) #erg/cm^2/spectral unit

#Plot the spectra 
plt.figure(1, figsize = (8,6))
plt.xlabel('Energy (keV)', size = 14)
plt.ylabel(r'Flux ($\mathrm{erg/cm}^{2}$)', size = 14)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
#plt.plot(E_obs, BB_em, color = 'b', label = 'BB_em from python')
#plt.plot(E_obs, BB_obs, color = 'g', label = 'BB_obs from python')
#plt.plot(E_obs, BB_only, color = 'k', label = 'BB_obs from flux_master')
plt.plot(E_obs, Flux, color = 'r', label = 'Flux from flux_master (with redshift etc)')
plt.plot(E_obs, Flux_obs, color = 'k', label = 'Flux from python')
plt.legend()

plt.show()