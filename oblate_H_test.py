#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 11:45:40 2019

@author: charlee
"""

"""Calculates flux for a relativistic star with an H atmosphere
in a spin-aligned system given inclination, T, and g."""

import numpy as np
import Interp_g_T_zeta as Ho
from Interp_g_T_zeta import lookup_alpha_MR
from get_R_theta import get_R
from get_g import get_g

log_T = 5.3 #Kelvin
T = 10**log_T
#Initial paramters
M = 2.78e33 #1.4 solar masses in grams
R_eq = 12e5 #cm
dist = 6.1713552e20 #200 pc in cm
spin_freq = 0

#List of E/kT values to calculate energies for
E_list_obs = np.linspace(0.06,60,166) 
#Define constants
h = 6.626e-27 #planck const in erg*s
c = 3e10 #speed of light cm/s 
k_B = 1.38e-16 #Boltzmann in erg/K
G = 6.674e-8 #grav constant in dyne cm^2 / g^2

T = 10**log_T
const_inte = (k_B*T)/h 

incl = np.radians(90)
            
#flux-master tables must be changed when dtheta/dphi is altered!
theta_all = np.linspace(1e-06, np.pi,10)
phi_all = np.linspace(-np.pi,np.pi, 20) #since phi space is double, double stepsize

E_dx = E_list_obs[1] - E_list_obs[0]
theta_dx = theta_all[1] - theta_all[0]
phi_dx = phi_all[1] - phi_all[0]

phi_integrand = np.zeros(len(phi_all))
theta_integrand = np.zeros(len(theta_all))
#LENGTH OF HO FILES IS 166 - MAKE SURE THIS IS STILL TRUE
spectral_flux = np.zeros(166)
spectral_I = np.zeros(166)

dOmega_all = np.loadtxt('flux-master/dOmega.txt', skiprows = 1, unpack = True).T
boost_all = np.loadtxt('flux-master/boost.txt', skiprows = 1, unpack = True).T
cos_beta_all = np.loadtxt('flux-master/cosbeta.txt', skiprows = 1, unpack = True).T


M_over_R = (G*M)/(c**2*R_eq) 
g = get_g(theta_all[1], R_eq, M, M_over_R, spin_freq)
log_g = np.round(np.log10(g), 2)

#Get H atm model for given temp and g
E, Inu, Inu_T_high, Inu_T_low, zeta_array,\
            zeta_g_high, zeta_g_low = Ho.interp_T_and_g(log_T, log_g)
            
for i in range(len(theta_all)):
    for j in range(len(phi_all) - 1):
    
        #R = get_R(theta_all[i], M, R_eq, spin_freq)
        #M_over_R = (G*M)/(c**2*R) #GM/Rc^2 - acceptable MoverR range 0.1 to 0.284
        #print("M over R =", M_over_R)
        dOmega = dOmega_all[i][j]
        #print("dOmega = ", dOmega)
       
        #acceptable g values between 13.7 and 14.7
        #g = get_g(theta_all[i], R, M, M_over_R, spin_freq)
        #log_g = np.round(np.log10(g), 2)
        #print("log_g =", log_g) 
        
        R = R_eq
            
        z = (1 - 2*M_over_R)**-0.5 - 1
        redshift = 1/(1+z)
        solid_const = R**2/dist**2
        v = spin_freq*2*np.pi*R #linear velocity 
        boost = boost_all[i][j]

        mu = cos_beta_all[i][j]
        
        #Calculate Doppler boost term
        #cos_xi = (-np.sin(alpha)*np.sin(incl)*np.sin(phi_all[j]))/(np.sin(psi))
        #boost = np.sqrt(1 - v**2/c**2)/(1 + (v/c)*cos_xi)
    
        #use mu instead of zeta in interpolation function 
        E_temp, Inu_temp = Ho.interp_zeta(mu, E, Inu, Inu_T_high,\
        Inu_T_low, zeta_array, zeta_g_high, zeta_g_low, log_g, log_T)
        
        E_list_em = E_list_obs*(redshift*boost)**-1
        
        Inu_final = Ho.interp_E(E_temp, Inu_temp, E_list_em)
        
        F_integrand = np.zeros(len(Inu_final))
        I_integrand = np.zeros(len(Inu_final))
        
        #exlude sections of rings we cannot see (negative values of mu)
        if mu > 0:
        #redshift**4 if using E_list_obs 
            for k in range(len(Inu_final)):
                
                F_integrand[k] = solid_const*const_inte*T**3*boost**4\
                                *Inu_final[k]*dOmega*redshift**3
                I_integrand[k] = solid_const*const_inte*T**3*boost**4\
                                *Inu_final[k]*redshift**3
                spectral_flux[k] = spectral_flux[k] + F_integrand[k]
                spectral_I[k] = spectral_I[k] + I_integrand[k]
        phi_integrand[j] = sum(F_integrand)*E_dx
    theta_integrand[i] = sum(phi_integrand)*phi_dx #total flux at each theta
bolo_flux = sum(theta_integrand)*theta_dx
         
print(bolo_flux)