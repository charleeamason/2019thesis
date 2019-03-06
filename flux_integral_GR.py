#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  7 13:12:49 2019

@author: charlee

A program which finds the flux emitted by a spherical, nonrotating neutron star
using general relativity in an SPIN-aligned coordinate system. 
A hydrogen atmosphere model from Ho et al. is compared to a simple blackbody.
"""

def get_flux_BB(inclination, log_T, E_list):
    
    """Calculates flux for a relativistic star with a BB atmosphere
    in a spin-aligned system given inclination and T."""

    import numpy as np
    from Interp_g_T_eta import lookup_alpha 
    
    #Define constants
    h = 6.626e-27 #planck const in erg*s
    c = 3e10 #speed of light cm/s 
    k_B = 1.38e-16 #Boltzmann in erg/K
    G = 6.674e-8 #grav constant in dyne cm^2 / g^2
    
    #Initial parameters:
    spin_freq = 0 #1/s
    M = 2.78e33 #1.4 solar masses in grams
    R = 2.06e6 #cm
    M_over_R = 0.1 #GM/Rc^2
    z = (1 - 2*M_over_R)**-0.5 - 1
    redshift = 1/(1+z) #gravitational redshift factor
    dist = 6.1713552e20 #200 pc in cm
    solid_const = R**2/dist**2
    
    T = 10**log_T
    const = 2/(c**2*h**2) 
    const_inte = (k_B*T)/h 
    
    incl = np.radians(inclination)
    
    theta_all = np.linspace(1e-06, np.pi, 40)
    phi_all = np.linspace(-np.pi,np.pi, 40)
    
    E_dx = E_list[1] - E_list[0]
    theta_dx = theta_all[1] - theta_all[0]
    phi_dx = phi_all[1] - phi_all[0]
    
    theta_integrand = np.zeros(len(theta_all))
    phi_integrand = np.zeros(len(phi_all))

    spectral_flux = np.zeros(len(E_list))
    spectral_I = np.zeros(len(E_list))
     
    for i in range(len(theta_all)):
        for j in range(len(phi_all)):
            #eta = cos psi 
            #print(i,j)
            eta = np.cos(incl)*np.cos(theta_all[i])\
                             + np.sin(incl)*np.sin(theta_all[i])*np.cos(phi_all[j])
           
            psi = np.arccos(eta)
            print('theta =', theta_all[i], 'phi =', phi_all[j])
            print('psi input =', psi)
            
            alpha, dmu_deta = lookup_alpha(M_over_R, psi)
            mu = np.cos(alpha)
            
            print('alpha =', alpha, 'dmu_deta =', dmu_deta)
            Inu = const*((k_B*T*E_list)**3)/(np.e**(E_list) - 1)
            F_integrand = np.zeros(len(Inu))
            
            #exlude sections of rings we cannot see (negative values of mu)
            if mu > 0:
                #F_integrand = np.zeros(len(Inu))
                for k in range(len(Inu)):
                    F_integrand[k] = (np.sin(theta_all[i])*mu*Inu[k]\
                                *const_inte*dmu_deta)
                    spectral_flux[k] = spectral_flux[k] + F_integrand[k]
                    spectral_I[k] = spectral_I[k] + Inu[k]*const_inte*redshift**3
            phi_integrand[j] = sum(F_integrand)*E_dx
        theta_integrand[i] = sum(phi_integrand)*phi_dx #total flux at each theta
    bolo_flux = solid_const*(1 - 2*M_over_R)*sum(theta_integrand)*theta_dx
        
    return bolo_flux, spectral_flux, spectral_I
    
def get_flux_H(inclination, log_T, log_g, E_list):
    """Calculates flux for a relativistic star with an H atmosphere
    in a spin-aligned system given inclination, T, and g."""

    import numpy as np
    import Interp_g_T_eta as Ho
    from Interp_g_T_eta import lookup_alpha
    
    #Define constants
    h = 6.626e-27 #planck const in erg*s
    c = 3e10 #speed of light cm/s 
    k_B = 1.38e-16 #Boltzmann in erg/K
    G = 6.674e-8 #grav constant in dyne cm^2 / g^2
    
    #Initial parameters:
    spin_freq = 0 #1/s
    M = 2.78e33 #1.4 solar masses in grams
    R = 2.06e6  #cm
    M_over_R = 0.1 #GM/Rc^2
    z = (1 - 2*M_over_R)**-0.5 - 1
    redshift = (1/(1+z)**3)
    dist = 6.1713552e20 #200 pc in cm
    solid_const = R**2/dist**2
    
    T = 10**log_T
    const_inte = (k_B*T)/h 
    
    incl = np.radians(inclination)
 
    #Get H atm model for given temp and g
    E, Inu, Inu_T_high, Inu_T_low, eta_array,\
                eta_g_high, eta_g_low = Ho.interp_T_and_g(log_T, log_g)
    
    theta_all = np.linspace(1e-06, np.pi,10)
    phi_all = np.linspace(-np.pi,np.pi, 10)
    
    E_dx = 0.2
    theta_dx = theta_all[1] - theta_all[0]
    phi_dx = phi_all[1] - phi_all[0]
    
    phi_integrand = np.zeros(len(phi_all))
    theta_integrand = np.zeros(len(theta_all))
    
    #LENGTH OF HO FILES IS 166 - MAKE SURE THIS IS STILL TRUE
    spectral_flux = np.zeros(166)
    spectral_I = np.zeros(166)
    for i in range(len(theta_all)):
        for j in range(len(phi_all)):
            #eta = cos psi 
            eta = np.cos(incl)*np.cos(theta_all[i])\
                            + np.sin(incl)*np.sin(theta_all[i])*np.cos(phi_all[j])
                            
            psi, alpha, dmu_deta = lookup_alpha(M_over_R, np.arccos(eta))
            mu = np.cos(alpha)
            
            #use mu instead of eta in interpolation function 
            E_temp, Inu_temp = Ho.interp_eta(mu, E, Inu, Inu_T_high,\
            Inu_T_low, eta_array, eta_g_high,eta_g_low, log_g, log_T)
        
            Inu_final = Ho.interp_E(E_temp, Inu_temp, E_list)
            
            #exlude sections of rings we cannot see (negative values of mu)
            if mu > 0:
                
                F_integrand = np.zeros(len(Inu_final))
                I_integrand = np.zeros(len(Inu_final))
                
                for k in range(len(Inu_final)):
                    F_integrand[k] = (mu*np.sin(theta_all[i])*\
                                     Inu_final[k]*const_inte*redshift)/(1-2*M_over_R)
                    I_integrand[k] = (Inu_final[k]*const_inte*redshift)/(1-2*M_over_R)
                    spectral_flux[k] = spectral_flux[k] + F_integrand[k]
                    spectral_I[k] = spectral_I[k] + I_integrand[k]
            phi_integrand[j] = sum(F_integrand)*E_dx
        theta_integrand[i] = sum(phi_integrand)*phi_dx #total flux at each theta
    bolo_flux = T**3*solid_const*sum(theta_integrand)*theta_dx
             
    return bolo_flux, spectral_flux, spectral_I