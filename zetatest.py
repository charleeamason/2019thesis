#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 12:04:36 2019

@author: charlee
"""
"""For comparing BB flux emisison from fixed zeta angles (Newton)"""

import numpy as np

#Enter T value between 5.1 and 6.8
log_T = 5.3 #Kelvin
T = 10**log_T
#Enter g value between 13.7 and 14.7
log_g = 14.4 #g = GM/R^2 * redshift^2
inclination = 0

#List of E/kT values to calculate energies for
E_list = np.linspace(0.06,100,166)

#Define constants
h = 6.626e-27 #planck const in erg*s
c = 3e10 #speed of light cm/s 
k_B = 1.38e-16 #Boltzmann in erg/K

#Initial parameters:
M = 2.78e33 #1.4 solar masses in grams
R = 2.06e6 #cm
dist = 6.1713552e20 #200 pc in cm
solid_const = R**2/dist**2

T = 10**log_T
const_inte = (k_B*T)/h 
const = 2/(c**2*h**2) 

incl = np.radians(inclination)

E_dx = E_list[1] - E_list[0]


#theta_all = np.linspace(1e-06, np.pi,10)
#phi_all = np.linspace(-np.pi,np.pi, 10)

#theta_dx = theta_all[1] - theta_all[0]
#phi_dx = phi_all[1] - phi_all[0]

#theta_integrand = np.zeros(len(theta_all))
#phi_integrand = np.zeros(len(phi_all))

spectral_flux = np.zeros(len(E_list))
spectral_I = np.zeros(len(E_list))

#for i in range(len(theta_all)):
    #for j in range(len(phi_all)):
        #zeta = cos psi (= cos alpha for newton)
zeta = 1

Inu = const*((k_B*T*E_list)**3)/(np.e**(E_list) - 1)
F_integrand = np.zeros(len(Inu))

#exlude sections of rings we cannot see (negative values of eta)

for k in range(len(Inu)):
    F_integrand[k]  = zeta*Inu[k]*const_inte*E_dx
    spectral_flux[k] = spectral_flux[k] + F_integrand[k]
    spectral_I[k] = spectral_I[k] + Inu[k]*const_inte
       # phi_integrand[j] = sum(F_integrand)*E_dx
   # theta_integrand[i] = sum(phi_integrand)*phi_dx #total flux at each theta
#bolo_flux = solid_const*sum(theta_integrand)*theta_dx
bolo_flux = sum(F_integrand)*solid_const
    
print(bolo_flux)