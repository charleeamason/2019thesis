#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 14:06:25 2019

@author: charlee
"""

"""Calculates flux for a Newtonian star with an H atmosphere
in a spin-aligned system given inclination, T, and g."""

import numpy as np
import Interp_g_T_zeta as Ho

#Enter T value between 5.1 and 6.8
log_T = 5.34 #Kelvin
T = 10**log_T
#Enter g value between 13.7 and 14.7
log_g = 14.36 #g = GM/R^2 * redshift^2

#List of E/kT values to calculate energies for
E_list = np.linspace(0.05,100,166)

inclination = 0

#Define constants
h = 6.626e-27 #planck const in erg*s
c = 3e10 #speed of light cm/s 
k_B = 1.38e-16 #Boltzmann in erg/K

#Initial parameters:
M = 2.78e33 #1.4 solar masses in grams
R = 2.06e6 #cm
dist = 6.1713552e20 #200 pc in cm
solid_const = R**2/dist**2
spin_freq = 600 #Hz
#v = spin_freq*2*np.pi*R #linear velocity 

T = 10**log_T
const_inte = (k_B*T)/h 

incl = np.radians(inclination)
 
#Get H atm model for given temp and g
E, Inu, Inu_T_high, Inu_T_low, zeta_array,\
            zeta_g_high, zeta_g_low = Ho.interp_T_and_g(log_T, log_g)

E_dx = 0.2
theta_all = np.linspace(1e-06, np.pi,10)
phi_all = np.linspace(-np.pi,np.pi, 10)

theta_dx = theta_all[1] - theta_all[0]
phi_dx = phi_all[1] - phi_all[0]

phi_integrand = np.zeros(len(phi_all))
theta_integrand = np.zeros(len(theta_all))

#LENGTH OF HO FILES IS 166 - MAKE SURE THIS IS STILL TRUE
spectral_flux = np.zeros(166)
spectral_I = np.zeros(166)
for i in range(len(theta_all)):
    for j in range(len(phi_all)):
        #zeta = cos psi (= cos alpha for newton)
        zeta = np.cos(incl)*np.cos(theta_all[i])\
                        + np.sin(incl)*np.sin(theta_all[i])*np.cos(phi_all[j])
        #Doppler boosting term 
       #boost = np.sqrt(1 - v**2/c**2)/(1 + (v/c)*np.sin(incl)*np.sin(phi_all[j]))   
                
        E_temp, Inu_temp = Ho.interp_zeta(zeta, E, Inu, Inu_T_high,\
        Inu_T_low, zeta_array, zeta_g_high, zeta_g_low, log_g, log_T)
    
        Inu_final = Ho.interp_E(E_temp, Inu_temp, E_list)
        
        #exlude sections of rings we cannot see (negative values of mu)
        F_integrand = np.zeros(len(Inu_final))
        I_integrand = np.zeros(len(Inu_final))
        if zeta > 0:
            for k in range(len(Inu_final)):
                F_integrand[k] = zeta*np.sin(theta_all[i])*\
                                 Inu_final[k]*const_inte*T**3 #*boost**3
                I_integrand[k] = Inu_final[k]*const_inte*T**3 #*boost**3
                spectral_flux[k] = spectral_flux[k] + F_integrand[k]
                spectral_I[k] = spectral_I[k] + I_integrand[k]
        phi_integrand[j] = sum(F_integrand)*E_dx
    theta_integrand[i] = sum(phi_integrand)*phi_dx #total flux at each theta
bolo_flux = solid_const*sum(theta_integrand)*theta_dx
         
print(bolo_flux)