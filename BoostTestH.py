#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 15:21:50 2019

@author: charlee
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 14:24:28 2019

@author: charlee
"""
import time
import matplotlib.pyplot as plt
import numpy as np
from Interp_g_T_zeta import lookup_alpha 
import Interp_g_T_zeta as Ho

"""Calculates flux for a relativistic star with a BB atmosphere
in a spin-aligned system given inclination and T."""
start_time = time.time()

log_T = 5.3 #Kelvin
T = 10**log_T

#Initial paramters
M = 2.78e33 #1.4 solar masses in grams
#R = 8.03544e5 #cm
R = 7.25453e5
G = 6.674e-8 #grav constant in dyne cm^2 / g^2
c = 3e10
#acceptable MoverR range 0.1 to 0.284
M_over_R_actual = (G*M)/(c**2*R) #GM/Rc^2
#M_over_R = 0.2564
M_over_R = 0.284
print("M over R =", M_over_R_actual)
print("M over R approx =", M_over_R)

#acceptable g values between 13.7 and 14.7
g = (G*M)/R**2*(1 - 2*M_over_R)**-0.5
#log_g = np.round(np.log10(g), 2)
log_g = 14.7
print("log_g =", log_g)

#List of E/kT values to calculate energies for
#E_list = np.linspace(0.06,100,166)
E_list_obs = np.linspace(0.06,60,166)

#Relevant angles (radians):
#***avoid |i| > 90 (quadrant issues)***
inclination = 90

#Define constants
h = 6.626e-27 #planck const in erg*s
c = 3e10 #speed of light cm/s 
k_B = 1.38e-16 #Boltzmann in erg/K
G = 6.674e-8 #grav constant in dyne cm^2 / g^2

#Initial parameters:
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

#Get H atm model for given temp and g
E, Inu, Inu_T_high, Inu_T_low, zeta_array,\
            zeta_g_high, zeta_g_low = Ho.interp_T_and_g(log_T, log_g)

theta_all = np.linspace(1e-06, np.pi,4)
phi_all = np.linspace(-np.pi,np.pi, 8) #since phi space is double, double stepsize

#E_dx = E_list[1] - E_list[0]
E_dx_obs = E_list_obs[1] - E_list_obs[0]
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
        
        E_list_em = E_list_obs*(redshift*boost)**-1

        Inu_final = Ho.interp_E(E_temp, Inu_temp, E_list_em)
        
        #exlude sections of rings we cannot see (negative values of mu)
        
        F_integrand = np.zeros(len(Inu_final))
        I_integrand = np.zeros(len(Inu_final))
            
        if mu > 0:
            
            for k in range(len(Inu_final)):
                
                F_integrand[k] = boost**4*T**3*(mu*np.sin(theta_all[i])*Inu_final[k]\
                                 *const_inte*dmu_dzeta*redshift**3)*(1-2*M_over_R)**-1
                           
                I_integrand[k] = boost**4*T**3*(Inu_final[k]*const_inte*dmu_dzeta*\
                                 redshift**3)*(1-2*M_over_R)**-1
                spectral_flux[k] = spectral_flux[k] + F_integrand[k]
                spectral_I[k] = spectral_I[k] + I_integrand[k]
        phi_integrand[j] = sum(F_integrand)*E_dx_obs
    theta_integrand[i] = sum(phi_integrand)*phi_dx #total flux at each theta
bolo_flux = solid_const*sum(theta_integrand)*theta_dx
             
print("bolo_flux =", bolo_flux)
print("--- %s seconds ---" % (time.time() - start_time))

#Plot the spectra
plt.figure(1, figsize = (8,6))
plt.title('log T = {0} log g = {1} incl = {2}'.format(log_T,log_g, inclination), size = 16)
plt.xlabel(r'$E/kT$ [dimensionless]', size = 16)
plt.ylabel(r'Spectral Flux [$\mathrm{erg/cm}^2\mathrm{/s/ster}$]', size = 16)
#plt.ylabel(r'$I_\nu \mathrm{(erg/cm}^2\mathrm{/ster}\mathrm{)}$', size = 16)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
#plt.ylim(0,2e17)
#plt.plot(E_list, spectral_I_BB, label = 'Pure BB')
plt.plot(E_list_em, spectral_flux, marker = 'o', label = 'Emitted energy' )
plt.plot(E_list_obs, spectral_flux, label = 'Observed energy' )
plt.legend()
plt.show()
#plt.savefig('FluxPlot_incl{0}_newt'.format(inclination))
#plt.clf()          

