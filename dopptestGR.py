#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  1 14:04:16 2019

@author: charlee
"""

"""
Created on Mon Jan  7 13:12:49 2019

@author: charlee

Calculates flux for a relativistic star with a BB atmosphere
in a spin-aligned system given inclination and T."""

import numpy as np
from Interp_g_T_zeta import lookup_alpha 

#Enter T value between 5.1 and 6.8
log_T = 5.34 #Kelvin
T = 10**log_T
#Enter g value between 13.7 and 14.7
log_g = 14.36 #g = GM/R^2 * redshift^2

#List of E/kT values to calculate energies for
E_list = np.linspace(0.05,100,166)

inclination = 90

#Define constants
h = 6.626e-27 #planck const in erg*s
c = 3e10 #speed of light cm/s 
k_B = 1.38e-16 #Boltzmann in erg/K
G = 6.674e-8 #grav constant in dyne cm^2 / g^2

#Initial parameters:
M = 2.78e33 #1.4 solar masses in grams
R = 2.06e6 #cm
M_over_R = 0.1 #GM/Rc^2
z = (1 - 2*M_over_R)**-0.5 - 1
redshift = 1/(1+z) #gravitational redshift factor
dist = 6.1713552e20 #200 pc in cm
solid_const = R**2/dist**2
spin_freq = 600 #Hz
v = spin_freq*2*np.pi*R #linear velocity 

T = 10**log_T
const = 2/(c**2*h**2) 
const_inte = (k_B*T)/h 

incl = np.radians(inclination)

theta_all = np.linspace(1e-06, np.pi, 10)
phi_all = np.linspace(-np.pi,np.pi, 10)

E_dx = E_list[1] - E_list[0]
theta_dx = theta_all[1] - theta_all[0]
phi_dx = phi_all[1] - phi_all[0]

theta_integrand = np.zeros(len(theta_all))
phi_integrand = np.zeros(len(phi_all))

spectral_flux = np.zeros(len(E_list))
spectral_I = np.zeros(len(E_list))
 
for i in range(len(theta_all)):
    for j in range(len(phi_all)):
        #zeta = cos psi 
        #print(i,j)
        zeta = np.cos(incl)*np.cos(theta_all[i])\
                         + np.sin(incl)*np.sin(theta_all[i])*np.cos(phi_all[j])
                         
        psi = np.arccos(zeta) #arccos only valid for 0 - 90
        print('theta =', theta_all[i], 'phi =', phi_all[j])
        print('psi input =', psi)
        
        alpha, dmu_dzeta = lookup_alpha(M_over_R, psi)
        mu = np.cos(alpha)
        
        #Calculate Doppler boost term
        cos_xi = (-np.sin(alpha)*np.sin(incl)*np.sin(phi_all[j]))/(np.sin(psi))
        boost = np.sqrt(1 - v**2/c**2)/(1 + (v/c)*cos_xi)
        
        print('alpha =', alpha, 'dmu_dzeta =', dmu_dzeta)
        Inu = const*((k_B*T*E_list)**3)/(np.e**(E_list) - 1)*boost**3
        F_integrand = np.zeros(len(Inu))
        
        #exlude sections of rings we cannot see (negative values of mu)
        if mu > 0:
            #F_integrand = np.zeros(len(Inu))
            for k in range(len(Inu)):
                F_integrand[k] = (np.sin(theta_all[i])*mu*Inu[k]\
                            *const_inte*dmu_dzeta)
                spectral_flux[k] = spectral_flux[k] + F_integrand[k]
                spectral_I[k] = spectral_I[k] + Inu[k]*const_inte*redshift**3
        phi_integrand[j] = sum(F_integrand)*E_dx
    theta_integrand[i] = sum(phi_integrand)*phi_dx #total flux at each theta
bolo_flux = solid_const*(1 - 2*M_over_R)*sum(theta_integrand)*theta_dx
    
print(bolo_flux)