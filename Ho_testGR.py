#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 15:39:36 2019

@author: charlee
"""

import numpy as np
import Interp_g_T_zeta as Ho
from Interp_g_T_zeta import lookup_alpha

#Define constants
h = 6.626e-27 #planck const in erg*s
c = 3e10 #speed of light cm/s 
k_B = 1.38e-16 #Boltzmann in erg/K
G = 6.674e-8 #grav constant in dyne cm^2 / g^2

inclination = 0
log_T = 5.3
log_g = 14.61
E_list = np.linspace(0.06, 100, 166)

#Initial parameters:
#cross check M, R, dist, and M_over_R with main.py
M = 2.78e33 #1.4 solar masses in grams
R = 8.03544e5 #cm
M_over_R = 0.2564 #GM/Rc^2
z = (1 - 2*M_over_R)**-0.5 - 1
redshift = 1/(1+z)
dist = 6.1713552e20 #200 pc in cm
solid_const = R**2/dist**2
spin_freq = 600 #Hz
v = spin_freq*2*np.pi*R #linear velocity 

T = 10**log_T
const_inte = (k_B*T)/h 

incl = np.radians(inclination)
 
#Get H atm model for given temp and g
E, Inu, Inu_T_high, Inu_T_low, zeta_array,\
            zeta_g_high, zeta_g_low = Ho.interp_T_and_g(log_T, log_g)

theta_all = np.linspace(1e-06, np.pi,10)
phi_all = np.linspace(-np.pi,np.pi, 20) #since phi space is double, double stepsize

E_dx = E_list[1] - E_list[0]
theta_dx = theta_all[1] - theta_all[0]
phi_dx = phi_all[1] - phi_all[0]

phi_integrand = np.zeros(len(phi_all))
theta_integrand = np.zeros(len(theta_all))

#LENGTH OF HO FILES IS 166 - MAKE SURE THIS IS STILL TRUE
spectral_flux = np.zeros(166)
spectral_I = np.zeros(166)
for i in range(len(theta_all)):
    for j in range(len(phi_all)):
        #zeta = cos psi 
        zeta = np.cos(incl)*np.cos(theta_all[i])\
                        + np.sin(incl)*np.sin(theta_all[i])*np.cos(phi_all[j])
                        
        psi = np.arccos(zeta)
        alpha, dmu_dzeta = lookup_alpha(M_over_R, psi)
        mu = np.cos(alpha)
        
        #Calculate Doppler boost term
        cos_xi = (-np.sin(alpha)*np.sin(incl)*np.sin(phi_all[j]))/(np.sin(psi))
        boost = np.sqrt(1 - v**2/c**2)/(1 + (v/c)*cos_xi)
    
        #use mu instead of zeta in interpolation function 
        E_temp, Inu_temp = Ho.interp_zeta(mu, E, Inu, Inu_T_high,\
        Inu_T_low, zeta_array, zeta_g_high, zeta_g_low, log_g, log_T)
    
        Inu_final = Ho.interp_E(E_temp, Inu_temp, E_list)
        
        #exlude sections of rings we cannot see (negative values of mu)
        
        F_integrand = np.zeros(len(Inu_final))
        I_integrand = np.zeros(len(Inu_final))
            
        if mu > 0:
            
            for k in range(len(Inu_final)):
                
                F_integrand[k] = T**3*mu*np.sin(theta_all[i])*Inu_final[k]\
                                 *const_inte*dmu_dzeta*redshift**3*(1-2*M_over_R)**-0.5
                print(theta_all[i], phi_all[j], Inu_final[k], mu, np.sin(theta_all[i]), T**3)
                print(F_integrand[k])
                break 
            
                I_integrand[k] = T**3*(Inu_final[k]*const_inte*\
                                 redshift**3)/(1-2*M_over_R)**0.5
                spectral_flux[k] = spectral_flux[k] + F_integrand[k]
                spectral_I[k] = spectral_I[k] + I_integrand[k]
        phi_integrand[j] = sum(F_integrand)*E_dx
    theta_integrand[i] = sum(phi_integrand)*phi_dx #total flux at each theta
bolo_flux = solid_const*sum(theta_integrand)*theta_dx