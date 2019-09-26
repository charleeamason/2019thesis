#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 12:41:21 2018

@author: charlee

A program which finds the flux emitted by a spherical, rotating (optional) neutron star
using Newtonian gravity only in an SPIN-aligned coordinate system. 
A hydrogen atmosphere model from Ho et al. is compared to a simple blackbody.
"""

def get_flux_BB(inclination, log_T, E_list, M, R, dist, spin_freq, ntheta, nphi, M_over_R):
    
    """Calculates flux for a Newtonian star with a BB atmosphere
    in a spin-aligned system given inclination and T."""

    import numpy as np
    
    #Define constants
    h = 6.626e-27 #planck const in erg*s
    c = 3e10 #speed of light cm/s 
    k_B = 1.38e-16 #Boltzmann in erg/K
    
    #Initial parameters:
    solid_const = R**2/dist**2
    v = spin_freq*2*np.pi*R #linear velocity 
    
    T = 10**log_T
    const = 2/(c**2*h**2) 
    const_inte = (k_B*T)/h 
    
    incl = np.radians(inclination)
    
    E_dx = E_list[1] - E_list[0]
    theta_all = np.linspace(1e-06, np.pi,ntheta)
    phi_all = np.linspace(-np.pi,np.pi, nphi)
    
    dtheta = theta_all[1] - theta_all[0]
    dphi = phi_all[1] - phi_all[0]
    
    print("dphi*dtheta =", dphi*dtheta)

    theta_integrand = np.zeros(len(theta_all))
    phi_integrand = np.zeros(len(phi_all))

    spectral_flux = np.zeros(len(E_list))

   # spectral_I = np.zeros(len(E_list))
    
    for i in range(len(theta_all)):
        for j in range(len(phi_all)):
            #zeta = cos psi (= cos alpha for newton)
            zeta = np.cos(incl)*np.cos(theta_all[i])\
                             + np.sin(incl)*np.sin(theta_all[i])*np.cos(phi_all[j])
            #Doppler boosting term
            #boost = np.sqrt(1 - v**2/c**2)/(1 + (v/c)*np.sin(incl)*np.sin(phi_all[j]))   
              
            #print(eta)
            Inu = const*((k_B*T*E_list)**3)/(np.e**(E_list) - 1) #*boost**3
            
            F_integrand = np.zeros(len(Inu))
            
            #exlude sections of rings we cannot see (negative values of eta)
            
            if zeta > 0:
                #F_integrand = np.zeros(len(Inu))
                for k in range(len(Inu)):
                    F_integrand[k]  = solid_const*zeta*np.sin(theta_all[i])*\
                                      Inu[k]*(dtheta*dphi)
                    spectral_flux[k] = spectral_flux[k] + F_integrand[k]
            phi_integrand[j] = sum(F_integrand)*E_dx
        theta_integrand[i] = sum(phi_integrand)*dphi #total flux at each theta
    bolo_flux = sum(theta_integrand)*dtheta*const_inte/(dtheta*dphi)
    
    return bolo_flux, spectral_flux
    
def get_flux_H(inclination, log_T, E_list, M, R, dist, spin_freq, ntheta, nphi, M_over_R):
    """Calculates flux for a Newtonian star with an H atmosphere
    in a spin-aligned system given inclination, T, and g."""

    import numpy as np
    import Interp_g_T_zeta as Ho
    from get_g import get_g
    #Define constants
    h = 6.626e-27 #planck const in erg*s
    c = 3e10 #speed of light cm/s 
    k_B = 1.38e-16 #Boltzmann in erg/K
    G = 6.67e-8
    
    #Initial parameters:
    solid_const = R**2/dist**2
    v = spin_freq*2*np.pi*R #linear velocity 
    M_over_R = (G*M)/(c**2*R)
     
    T = 10**log_T
    const_inte = (k_B*T)/h 
    
    incl = np.radians(inclination)
    
    E_dx = E_list[1] - E_list[0]
    theta_all = np.linspace(1e-06, np.pi,ntheta)
    phi_all = np.linspace(-np.pi,np.pi, nphi)
    
    dtheta = theta_all[1] - theta_all[0]
    dphi= phi_all[1] - phi_all[0]
    
    phi_integrand = np.zeros(len(phi_all))
    theta_integrand = np.zeros(len(theta_all))
    
    #g same for all values of theta for spherical star
    g = get_g(theta_all[0], R, M, M_over_R, spin_freq)
    log_g = np.round(np.log10(g), 2)
    #print("g =", log_g)
    
    #Get H atm model for given temp and g
    E, Inu, Inu_T_high, Inu_T_low, zeta_array,\
                zeta_g_high, zeta_g_low = Ho.interp_T_and_g(log_T, log_g)

    #LENGTH OF HO FILES IS 166 - MAKE SURE THIS IS STILL TRUE
    spectral_flux = np.zeros(166)
    spectral_I = np.zeros(166)
    for i in range(len(theta_all)):
        for j in range(len(phi_all)):
            #zeta = cos psi (= cos alpha for newton)
            zeta = np.cos(incl)*np.cos(theta_all[i])\
                            + np.sin(incl)*np.sin(theta_all[i])*np.cos(phi_all[j])
            #Doppler boost term
           # boost = np.sqrt(1 - v**2/c**2)/(1 + (v/c)*np.sin(incl)*np.sin(phi_all[j]))
                    
            E_temp, Inu_temp = Ho.interp_zeta(zeta, E, Inu, Inu_T_high,\
            Inu_T_low, zeta_array, zeta_g_high,zeta_g_low, log_g, log_T)
        
            Inu_final = Ho.interp_E(E_temp, Inu_temp, E_list)
            
            #exlude sections of rings we cannot see (negative values of zeta)
            F_integrand = np.zeros(len(Inu_final))
            
            if zeta > 0:
                
                for k in range(len(Inu_final)):
                    
                    F_integrand[k] = solid_const*T**3*zeta*np.sin(theta_all[i])*\
                                     Inu_final[k]*dtheta*dphi #*boost**3
                    spectral_flux[k] = spectral_flux[k] + F_integrand[k]
            phi_integrand[j] = sum(F_integrand)*E_dx
        theta_integrand[i] = sum(phi_integrand)*dphi #total flux at each theta
    bolo_flux = sum(theta_integrand)*dtheta*const_inte/(dtheta*dphi)
             
    return bolo_flux, spectral_flux