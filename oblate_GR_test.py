#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 11:33:51 2019

@author: charlee
"""

"""Calculates flux for a relativistic star with a BB atmosphere
in a spin-aligned system given inclination and T."""

import numpy as np
from Interp_g_T_zeta import lookup_alpha_MR 
from get_R_theta import get_R
from get_g import get_g

#Define constants
h = 6.626e-27 #planck const in erg*s
c = 3e10 #speed of light cm/s 
k_B = 1.38e-16 #Boltzmann in erg/K
G = 6.674e-8 #grav constant in dyne cm^2 / g^2
log_T = 5.3 #Kelvin
T = 10**log_T
#Initial paramters
M = 2.78e33 #1.4 solar masses in grams
R_eq = 12e5 #cm
dist = 6.1713552e20 #200 pc in cm
spin_freq = 0
inclination = 90

#List of E/kT values to calculate energies for
E_list_obs = np.linspace(0.06,60,166)
#Initial parameters:
T = 10**log_T
const = 2/(c**2*h**2) 
const_inte = (k_B*T)/h 

incl = np.radians(inclination)

#flux master tables must be altered when dtheta, dphi is changed!
theta_all = np.linspace(1e-06, np.pi, 30)
phi_all = np.linspace(-np.pi,np.pi, 60)

E_dx = E_list_obs[1] - E_list_obs[0]

spectral_flux = np.zeros(len(E_list_obs))
spectral_I = np.zeros(len(E_list_obs))

#theta_all, phi_all = np.loadtxt('flux-master/angles.txt',skiprows = 1, unpack = True)
dOmega_all = np.loadtxt('flux-master/dOmega.txt', skiprows = 1, unpack = True).T
boost_all = np.loadtxt('flux-master/boost.txt', skiprows = 1, unpack = True).T
cos_beta_all = np.loadtxt('flux-master/cosbeta.txt', skiprows = 1, unpack = True).T

theta_dx = theta_all[1] - theta_all[0]
phi_dx = phi_all[1] - phi_all[0]

theta_integrand = np.zeros(len(theta_all))
theta_integrand2 = np.zeros(len(theta_all))
phi_integrand = np.zeros(len(phi_all))
phi_integrand2 = np.zeros(len(phi_all))

for i in range(len(theta_all)):
    for j in range(len(phi_all) - 1):
        #zeta = cos psi 
        #print(i,j)
        R = get_R(theta_all[i], M, R_eq, spin_freq)
        M_over_R = (G*M)/(c**2*R) #GM/Rc^2 - acceptable MoverR range 0.1 to 0.284
        #print("M over R =", M_over_R)
        dOmega = dOmega_all[i][j]
        #print("dOmega = ", dOmega)
        
        #acceptable g values between 13.7 and 14.7
        g = get_g(theta_all[i], R_eq, M, M_over_R, spin_freq)
        log_g = np.round(np.log10(g), 2)
        #print("log_g =", log_g)
        
        z = (1 - 2*M_over_R)**-0.5 - 1
        redshift = 1/(1+z) #gravitational redshift factor
        dist = 6.1713552e20 #200 pc in cm
        solid_const = R**2/dist**2
        v = spin_freq*2*np.pi*R #linear velocity 
        boost = boost_all[i][j]
        
        
        #zeta = np.cos(incl)*np.cos(theta_all[i])\
        #                 + np.sin(incl)*np.sin(theta_all[i])*np.cos(phi_all[j])
       
        #psi = np.arccos(zeta) #arccos only valid for 0 - 90
       # print('theta =', theta_all[i], 'phi =', phi_all[j])
        #print('psi input =', psi)
        
        #alpha, dmu_dzeta = lookup_alpha_MR(M_over_R, psi)
        
        #dOmega_2 = (np.sin(theta_all[i])*np.cos(alpha)*dmu_dzeta)*(1-2*M_over_R)**-1
        #print("dOmega 2 =", dOmega_2)
        
        #mu_alpha = np.cos(alpha)
        #print("theta =", theta_all[i])
        #print("phi =", phi_all[j])
        #print("mu_alpha =", mu_alpha)
        mu = cos_beta_all[i][j]
        #print("mu_beta =", mu)
        
        #Calculate Doppler boost term
        #cos_xi = (-np.sin(alpha)*np.sin(incl)*np.sin(phi_all[j]))/(np.sin(psi))
        #boost2 = np.sqrt(1 - v**2/c**2)/(1 + (v/c)*cos_xi)
        #print('boost calculated =', boost2)
        #print('boost from flux-master =', boost)
        
        #print('alpha =', alpha, 'dmu_dzeta =', dmu_dzeta)
       
        #need redshifted/observed energy list for spinning star 
        E_list_em = E_list_obs*(redshift*boost)**-1
        
        Inu = const*((k_B*T*E_list_em)**3)/(np.e**(E_list_em) - 1)
        #Inu = const*((k_B*T*E_list_obs)**3)/(np.e**(E_list_obs) - 1)
        #Inu_obs = Inu*(E_list_obs/E_list)**3
        
        F_integrand = np.zeros(len(Inu))
        #F_integrand2 = np.zeros(len(Inu))
        Inu_integrand = np.zeros(len(Inu))
        
        #exlude sections of rings we cannot see (negative values of mu)
        if mu > 0:
            #F_integrand = np.zeros(len(Inu))
            for k in range(len(Inu)):
                #redshift**4 if using E_list_obs 
                F_integrand[k] = solid_const*const_inte*boost**4*Inu[k]*dOmega*redshift**3
                
                #F_integrand2[k] = solid_const*boost2**4*(np.sin(theta_all[i])*mu_alpha*Inu[k]\
                                #*const_inte*dmu_dzeta*redshift**3)*(1-2*M_over_R)**-1
                #Inu_integrand[k] = boost**4*(Inu[k]*const_inte*dmu_dzeta*\
                             #redshift**3)*(1-2*M_over_R)**-1
                spectral_flux[k] = spectral_flux[k] + F_integrand[k]
                #spectral_I[k] = spectral_I[k] + Inu_integrand[k]
        phi_integrand[j] = sum(F_integrand)*E_dx
        #phi_integrand2[j] = sum(F_integrand2)*E_dx
    theta_integrand[i] = sum(phi_integrand)*phi_dx #total flux at each theta
    #theta_integrand2[i] = sum(phi_integrand2)*phi_dx
bolo_flux = sum(theta_integrand)*theta_dx
#bolo_flux_2 = sum(theta_integrand2)*theta_dx
    
print(bolo_flux)
