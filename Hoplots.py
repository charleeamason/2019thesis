#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 12:29:21 2019

@author: charlee
"""

"""Ho curves to compare to V.E. Zavlin (1996) H atmosphere models (fixed zeta)"""

import numpy as np
import Interp_g_T_zeta as Ho
import matplotlib.pyplot as plt
#Enter a value between 5.10 and 6.50
log_T = 5.3 #Kelvin
T = 10**log_T
#Enter g value between 13.7 and 14.7
log_g = 14.4 #g = GM/R^2 * redshift^2. 
inclination = 0

#logE_all, mu_all, logInu_over_T3_all, logT_all,logg_all =\
#        np.loadtxt("nsx_H_v171019.txt", usecols = (0,1,2,3,4), unpack = True)

#Define constants
h = 6.626e-27 #planck const in erg*s
c = 3e10 #speed of light cm/s 
k_B = 1.38e-16 #Boltzmann in erg/K

#Define Energy range
E_list = np.linspace(0.06, 100, 166)
E_list_Zavlin = np.linspace(0.058,581,166) #E/kT in cgs
erg_to_kev = 6.242e8
E_list_kev = E_list*erg_to_kev*k_B*T #E in keV
log_E_list_kev = np.log10(E_list_kev)
#print(E_list_kev)

#Initial parameters:
M = 2.78e33 #1.4 solar masses in grams
R = 2.06e6 #cm
dist = 6.1713552e20 #200 pc in cm
solid_const = R**2/dist**2

T = 10**log_T
const_inte = (k_B*T)/h 
const = 2/(c**2*h**2) 

incl = np.radians(inclination)
 
#Get H atm model for given temp and g
E, Inu, Inu_T_high, Inu_T_low, zeta_array,\
            zeta_g_high, zeta_g_low = Ho.interp_T_and_g(log_T, log_g)

E_dx = E_list[1] - E_list[0]
#theta_all = np.linspace(1e-06, np.pi,10)
#phi_all = np.linspace(-np.pi,np.pi, 10)

#theta_dx = theta_all[1] - theta_all[0]
#phi_dx = phi_all[1] - phi_all[0]

#phi_integrand = np.zeros(len(phi_all))
#theta_integrand = np.zeros(len(theta_all))

#LENGTH OF HO FILES IS 166 - MAKE SURE THIS IS STILL TRUE
spectral_flux = np.zeros(166)
spectral_flux_BB = np.zeros(166)
spectral_I = np.zeros(166)
#for i in range(len(theta_all)):
    #for j in range(len(phi_all)):
        #zeta = cos psi (= cos alpha for newton)
zeta = 0.17 #zeta = 1 normal to surface, zeta = 0 tangent to surface
      
E_temp, Inu_temp = Ho.interp_zeta(zeta, E, Inu, Inu_T_high,\
Inu_T_low, zeta_array, zeta_g_high, zeta_g_low, log_g, log_T)

Inu_final = Ho.interp_E(E_temp, Inu_temp, E_list)
Inu_BB = const*((k_B*T*E_list)**3)/(np.e**(E_list) - 1)

#exlude sections of rings we cannot see (negative values of mu)
F_integrand = np.zeros(len(Inu_final))
F_integrand_BB = np.zeros(len(Inu_BB))
I_integrand = np.zeros(len(Inu_final))
for k in range(len(Inu_final)):
    F_integrand[k] = zeta*Inu_final[k]*const_inte*T**3*E_dx
    F_integrand_BB[k] = zeta*Inu_BB[k]*const_inte*E_dx
                                  
    I_integrand[k] = Inu_final[k]*const_inte*T**3*E_dx
    
    spectral_flux[k] = spectral_flux[k] + F_integrand[k]
    
    spectral_flux_BB[k] = spectral_flux_BB[k] + F_integrand_BB[k]
    spectral_I[k] = spectral_I[k] + I_integrand[k]
        #phi_integrand[j] = sum(F_integrand)*E_dx
    #theta_integrand[i] = sum(phi_integrand)*phi_dx #total flux at each theta
#bolo_flux = solid_const*sum(theta_integrand)*theta_dx
bolo_flux = sum(F_integrand)*solid_const      
bolo_flux_BB = sum(F_integrand_BB)*solid_const 
print('Hydrogen =', bolo_flux)
print('BB =', bolo_flux_BB)
log_spectral_flux = np.log10(spectral_flux/(k_B*T*erg_to_kev))
log_spectral_flux_BB = np.log10(spectral_flux_BB/(k_B*T*erg_to_kev))

#Plot the spectra
plt.figure(1, figsize = (8,6))
#plt.title('log T = {0} log g ={1}'.format(log_T,log_g), size = 16)
plt.xlabel(r'$\mathrm{Log\ E}$ (keV)', size = 16)
plt.ylabel(r'$\mathrm{Log\ Spectral\ Flux}\mathrm{(erg/cm}^2\mathrm{/s}\mathrm{/keV}\mathrm{)}$', size = 16)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.ylim(13,20)
plt.xlim(-3,1)
plt.plot(log_E_list_kev, log_spectral_flux,\
         label = 'log T = {0} log g ={1}'.format(log_T,log_g))
plt.plot(log_E_list_kev, log_spectral_flux_BB,\
         label = 'BB log T = {0} log g ={1}'.format(log_T,log_g))
#plt.plot(E_list, spectral_flux_BB, marker = 'o', label = 'Pure BB flux' )
#plt.plot(E_list, spectral_flux*T**3, label = 'W.Ho' )
plt.legend()
plt.show() 