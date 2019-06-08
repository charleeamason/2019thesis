#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 11:25:51 2019

@author: charlee

A program which finds the flux emitted by an oblate, rotating neutron star
using general relativity in an SPIN-aligned coordinate system. 
A hydrogen atmosphere model from Ho et al. is compared to a simple blackbody.
"""

def get_flux_BB(inclination, log_T, E_list_obs, M, R_eq, dist, spin_freq):
    
    """Calculates flux for an oblate relativistic star with a BB atmosphere
    in a spin-aligned system given inclination and T."""
    
    import numpy as np
    from get_R_theta import get_R
    from get_g import get_g
    
    #Define constants
    h = 6.626e-27 #planck const in erg*s
    c = 3e10 #speed of light cm/s 
    k_B = 1.38e-16 #Boltzmann in erg/K
    G = 6.674e-8 #grav constant in dyne cm^2 / g^2

    #Initial parameters:
    T = 10**log_T
    const = 2/(c**2*h**2) 
    const_inte = (k_B*T)/h 

    E_dx = E_list_obs[1] - E_list_obs[0]
    
    dtheta = np.loadtxt('flux-master/angles_sph_spin.txt',skiprows = 1, usecols = -1, unpack = True)[0]
    dOmega_all = np.loadtxt('flux-master/dOmega_sph_spin.txt', skiprows = 1, unpack = True).T
    boost_all = np.loadtxt('flux-master/boost_sph_spin.txt', skiprows = 1, unpack = True).T
    cos_beta_all = np.loadtxt('flux-master/cosbeta_sph_spin.txt', skiprows = 1, unpack = True).T
    
    #flux master tables must be altered when ntheta, nphi is changed!
    ntheta = 30 #number of theta divisions
    nphi = 60 #number of phi divisions 
    dphi = 2*np.pi/nphi

    theta_all = np.zeros(ntheta)
    phi_all = np.zeros(nphi)

    #make theta list
    for m in range(2): #loop through hemispheres
        for x in range(int(ntheta/2)): #loop through latitudes 
            #print(x)
            if m == 0:
                #northern hemisphere
                theta_all[x] = x*dtheta + 1e-06 + 0.5*dtheta
            else:
                #southern hemisphere
                theta_all[x + int(ntheta/2)] = np.pi/2 + ((x + 1)*dtheta + 1e-06) - 0.5*dtheta
           
    #make phi list 
    for y in range(nphi):
        phi_all[y] = 1e-06 + y*dphi
        if y == nphi - 1:
            #endpoints identical
            phi_all[y] = 1e-06
            
    spectral_flux = np.zeros(len(E_list_obs))
    spectral_I = np.zeros(len(E_list_obs))
    theta_integrand = np.zeros(len(theta_all))
    theta_integrand_SA = np.zeros(len(theta_all))
    phi_integrand = np.zeros(len(phi_all))
    phi_integrand_SA = np.zeros(len(phi_all))
         
    for i in range(len(theta_all)):
        for j in range(len(phi_all) - 1):
            #print(i,j)
            R = get_R(theta_all[i], M, R_eq, spin_freq)
            M_over_R = (G*M)/(c**2*R) #GM/Rc^2 - acceptable MoverR range 0.1 to 0.284
            #print("M over R =", M_over_R)
            dOmega = dOmega_all[i][j]
            #print("dOmega = ", dOmega)
            
            #acceptable g values between 13.7 and 14.7
            g = get_g(theta_all[i], R, M, M_over_R, spin_freq)
            log_g = np.round(np.log10(g), 2)
            #print("log_g =", log_g)
            
            z = (1 - 2*M_over_R)**-0.5 - 1
            redshift = 1/(1+z) #gravitational redshift factor
            solid_const = R**2/dist**2
            v = spin_freq*2*np.pi*R #linear velocity 
            boost = boost_all[i][j]
            
            
            mu = cos_beta_all[i][j]
            #print("mu_beta =", mu)
           
            #need redshifted/observed energy list for spinning star
           # boost = 1
            #redshift = 1
            E_list_em = E_list_obs*(redshift*boost)**-1
            E_dx_em = E_dx*(redshift*boost)**-1
            #print(E_list_em)
            Inu = const*((k_B*T*E_list_em)**3)/(np.e**(E_list_em) - 1)
            #Inu = const*((k_B*T*E_list_obs)**3)/(np.e**(E_list_obs) - 1)
            #Inu_obs = Inu*(E_list_obs/E_list)**3
            
            F_integrand = np.zeros(len(Inu))
            Inu_integrand = np.zeros(len(Inu))
            SA_integrand = np.zeros(len(Inu))
            SA_integrand2 = np.zeros(len(Inu))    
            
            #redshift = 1
            #exlude sections of rings we cannot see (negative values of mu)
            if mu > 0:
                #SA_integrand2 = (dOmega*solid_const)
                #print("Solid angle 2 =", SA_integrand2)
                #print("E_initial =", E_list_em[0], "E_final =", E_list_em[-1])
                for k in range(len(Inu)):
                    #redshift**4 if using E_list_obs 
                    #Inu[k] = 1
                    F_integrand[k] = solid_const*const_inte*boost**4*redshift**3*Inu[k]*dOmega
                    Inu_integrand[k] = solid_const*const_inte*boost**4*redshift**3*Inu[k]
                    SA_integrand[k] = SA_integrand[k] + (dOmega*solid_const)\
                    /(E_list_em[-1] - E_list_em[0] + E_dx_em)
                    spectral_flux[k] = spectral_flux[k] + F_integrand[k]/(dphi*dtheta)
                    spectral_I[k] = spectral_I[k] + Inu_integrand[k]
            phi_integrand[j] = sum(F_integrand)*E_dx_em
            phi_integrand_SA[j] = sum(SA_integrand)*E_dx_em
            #print("Solid angle integrand =", sum(SA_integrand)*E_dx_em)
        theta_integrand[i] = sum(phi_integrand)*dphi #total flux at each theta
        theta_integrand_SA[i] = sum(phi_integrand_SA) 
    bolo_flux = sum(theta_integrand)*dtheta/(dphi*dtheta) #dOmega already has dphi dtheta in it
    SA = sum(theta_integrand_SA)
    print("Solid angle =", SA)
    
    return bolo_flux, spectral_flux, spectral_I, log_g, M_over_R
    
def get_flux_H(inclination, log_T, E_list_obs, M, R_eq, dist, spin_freq):
    """Calculates flux for a relativistic star with an H atmosphere
    in a spin-aligned system given inclination, T, and g."""

    import numpy as np
    import Interp_g_T_zeta as Ho
    #from Interp_g_T_zeta import lookup_alpha_MR
    from get_R_theta import get_R
    from get_g import get_g
  
    #Define constants
    h = 6.626e-27 #planck const in erg*s
    c = 3e10 #speed of light cm/s 
    k_B = 1.38e-16 #Boltzmann in erg/K
    G = 6.674e-8 #grav constant in dyne cm^2 / g^2
    
    T = 10**log_T
    const_inte = (k_B*T)/h 
                
    dtheta = np.loadtxt('flux-master/angles.txt',skiprows = 1, usecols = -1, unpack = True)[0]
    dOmega_all = np.loadtxt('flux-master/dOmega.txt', skiprows = 1, unpack = True).T
    boost_all = np.loadtxt('flux-master/boost.txt', skiprows = 1, unpack = True).T
    cos_beta_all = np.loadtxt('flux-master/cosbeta.txt', skiprows = 1, unpack = True).T
    
    #flux master tables must be altered when ntheta, nphi is changed!
    ntheta = 30 #number of theta divisions
    nphi = 60 #number of phi divisions 
    dphi = 2*np.pi/nphi

    theta_all = np.zeros(ntheta)
    phi_all = np.zeros(nphi)

    #make theta list
    for m in range(2): #loop through hemispheres
        for x in range(int(ntheta/2)): #loop through latitudes 
            #print(x)
            if m == 0:
                #northern hemisphere
                theta_all[x] = x*dtheta + 1e-06 + 0.5*dtheta
            else:
                #southern hemisphere
                theta_all[x + int(ntheta/2)] = np.pi/2 + ((x + 1)*dtheta + 1e-06) - 0.5*dtheta
           
    #make phi list 
    for y in range(nphi):
        phi_all[y] = 1e-06 + y*dphi
        if y == nphi - 1:
            #endpoints identical
            phi_all[y] = 1e-06
            
    theta_integrand = np.zeros(len(theta_all))
    phi_integrand = np.zeros(len(phi_all))
    
    E_dx = E_list_obs[1] - E_list_obs[0]
    #LENGTH OF HO FILES IS 166 - MAKE SURE THIS IS STILL TRUE
    spectral_flux = np.zeros(166)
    spectral_I = np.zeros(166)

    for i in range(len(theta_all)):
        for j in range(len(phi_all) - 1):
            
            R = get_R(theta_all[i], M, R_eq, spin_freq)
            M_over_R = (G*M)/(c**2*R) #GM/Rc^2 - acceptable MoverR range 0.1 to 0.284
            #print("M over R =", M_over_R)
            dOmega = dOmega_all[i][j]
            #print("dOmega = ", dOmega)
           
            #acceptable g values between 13.7 and 14.7
            g = get_g(theta_all[i], R, M, M_over_R, spin_freq)
            log_g = np.round(np.log10(g), 2)
            #print("log_g =", log_g)
            
            #Get H atm model for given temp and g
            E, Inu, Inu_T_high, Inu_T_low, zeta_array,\
                zeta_g_high, zeta_g_low = Ho.interp_T_and_g(log_T, log_g)
                
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
            
            #boost = 1
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
                    spectral_flux[k] = spectral_flux[k] + F_integrand[k]/(dtheta*dphi)
                    spectral_I[k] = spectral_I[k] + I_integrand[k]
            phi_integrand[j] = sum(F_integrand)*E_dx
        theta_integrand[i] = sum(phi_integrand)*dphi #total flux at each theta
    bolo_flux = sum(theta_integrand)*dtheta/(dtheta*dphi)
             
    return bolo_flux, spectral_flux, spectral_I, log_g, M_over_R