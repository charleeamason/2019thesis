"""Calculates flux for a Newtonian star with an H atmosphere
in a spin-aligned system given inclination, T, and g."""

import numpy as np
import Interp_g_T_zeta as Ho
#Define constants
h = 6.626e-27 #planck const in erg*s
c = 3e10 #speed of light cm/s 
k_B = 1.38e-16 #Boltzmann in erg/K
inclination = 0
log_T = 5.3
log_g = 14.61
E_list = np.linspace(0.06, 100, 166)

#Initial parameters:
#cross check M, R, and dist with main.py
R = 8.03544e5 #cm
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

E_dx = E_list[1] - E_list[0]
theta_all = np.linspace(1e-06, np.pi,10)
phi_all = np.linspace(-np.pi,np.pi, 20)

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
        #Doppler boost term
        boost = np.sqrt(1 - v**2/c**2)/(1 + (v/c)*np.sin(incl)*np.sin(phi_all[j]))
                
        E_temp, Inu_temp = Ho.interp_zeta(zeta, E, Inu, Inu_T_high,\
        Inu_T_low, zeta_array, zeta_g_high,zeta_g_low, log_g, log_T)
    
        Inu_final = Ho.interp_E(E_temp, Inu_temp, E_list)
        
        #exlude sections of rings we cannot see (negative values of zeta)
        F_integrand = np.zeros(len(Inu_final))
        I_integrand = np.zeros(len(Inu_final))
        
        if zeta > 0:
            
            for k in range(len(Inu_final)):
                
                F_integrand[k] = T**3*zeta*np.sin(theta_all[i])*\
                                 Inu_final[k]*const_inte #*boost**3
                print(theta_all[i], phi_all[j])
                print(F_integrand[k])
                break
                I_integrand[k] = T**3*Inu_final[k]*const_inte #*boost**3
                spectral_flux[k] = spectral_flux[k] + F_integrand[k]
                spectral_I[k] = spectral_I[k] + I_integrand[k]
        phi_integrand[j] = sum(F_integrand)*E_dx
    theta_integrand[i] = sum(phi_integrand)*phi_dx #total flux at each theta
bolo_flux = solid_const*sum(theta_integrand)*theta_dx