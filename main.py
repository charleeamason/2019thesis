#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 13:45:13 2018

@author: charlee

A program which finds the flux emitted by a given NS.
A hydrogen atmosphere model from Ho et al. is compared to a simple blackbody.
"""
import numpy as np
import time
import sys
import matplotlib.pyplot as plt

start_time = time.time()
#T value should be between 5.1 and 6.8
log_T = 6.0 #Kelvin
T = 10**log_T

#Initial paramters
M = 2.78e33 #1.4 solar masses in grams
#R_eq = 12e5 #cm - for M/R = 0.17
R_eq = 9.37e5 #cm - for M/R = 0.22
#R_eq = 7.64e5 #cm = for M/R = 0.27
dist = 6.1713552e20 #200 pc in cm
spin_freq = 0 

c = 3e10 #speed of light cm/s 
G = 6.674e-8 #grav constant in dyne cm^2 / g^2
h = 6.626e-27 #planck const, erg s
k = 1.38e-16 #blotzmann const erg/K

print((G*M)/(c**2*R_eq))

M_over_R = np.round((G*M)/(c**2*R_eq), 2)
print("M over Req (rounded to 2 decimal places)= ", M_over_R)


#List of E/kT values to calculate energies for
min_E = 0.06
max_E = 60
E_list = np.linspace(min_E,max_E,166)

#flux master tables must be altered when ntheta, nphi is changed!
ntheta = 30 #number theta bins
nphi = 2*ntheta #number phi bins
print('ntheta =', ntheta)

#inclination angles:
#***avoid |i| > 90 (quadrant issues)***
incl_array = [0,30,45,60,90]
#incl_array = np.arange(0,95,5)
#incl_array = np.insert(incl_array, 1, [1,3])

newt = True #Newtonian gravity? 
oblate = False #Is star oblate?

if newt == True:
    #Newtonian no spin
    fname = "newt"
    import flux_Newton as f
    if spin_freq != 0:
        print("spin turned on with Newtonian gravity")
        sys.exit()
        
elif oblate == True:
    #Oblate with spin
    fname = "obl_spin" + str(spin_freq)
    import flux_GR_oblate as f
    if spin_freq == 0:
        print("oblateness turned on without spin")
        sys.exit()
        
else: 
    import flux_GR_sphere as f
    fname = "sph_spin" + str(spin_freq)

        
print(fname)

bolo_flux_array = []
H_bolo_flux_array = []
for inclination in incl_array:
    print("inclination = ", inclination)
    if newt == True: 
    #Newtonian 
        bolo_flux_BB, spectral_flux_BB =\
           f.get_flux_BB(inclination, log_T, E_list, M, R_eq, dist, spin_freq, ntheta, nphi, M_over_R)
        
        bolo_flux_Ho, spectral_flux =\
          f.get_flux_H(inclination, log_T, E_list, M, R_eq, dist, spin_freq, ntheta, nphi, M_over_R)
    else:                      
    #GR
        bolo_flux_BB, spectral_flux_BB, log_g, M_over_R =\
            f.get_flux_BB(inclination, log_T, E_list, M, R_eq, dist, spin_freq, fname, ntheta, nphi, M_over_R)
        bolo_flux_Ho, spectral_flux, log_g, M_over_R =\
           f.get_flux_H(inclination, log_T, E_list, M, R_eq, dist, spin_freq, fname, ntheta, nphi, M_over_R)                        
       
    bolo_flux_array.append(bolo_flux_BB)
    H_bolo_flux_array.append(bolo_flux_Ho)

    np.savez('./fluxdata/ntheta' + str(ntheta) + '/BB_spectra_' + fname + '_MR{0}_incl{1}'\
             .format(int(100*np.round(M_over_R,2)), inclination),\
             E_list, spectral_flux_BB)
    
    np.savez('./fluxdata/ntheta' + str(ntheta) + '/H_spectra_' + fname + '_MR{0}_incl{1}'\
            .format(int(100*np.round(M_over_R,2)), inclination),\
              E_list, spectral_flux)
  
    np.savez('./fluxdata/ntheta' + str(ntheta) + '/BB_bolo_' + fname +'_MR{0}'\
             .format(int(100*np.round(M_over_R,2)), inclination),\
          incl_array, bolo_flux_array)
          
    np.savez('./fluxdata/ntheta' + str(ntheta) + '/H_bolo_' + fname +'_MR{0}'\
            .format(int(100*np.round(M_over_R,2)), inclination),\
         incl_array, H_bolo_flux_array)


    print("--- %s seconds ---" % (time.time() - start_time))

print("--- %s seconds ---" % (time.time() - start_time))
print('flux Ho =',  H_bolo_flux_array)
print('BB flux =', bolo_flux_array)

spectral_flux_peak = np.max(spectral_flux_BB)
spectral_flux_peak_analytical = (R_eq/dist)**2*(2*np.pi/(c**2*h**2))*(3*k*T)**3*(np.e**3 - 1)**(-1)
print("analytical spectral flux peak =", spectral_flux_peak_analytical)
spectral_flux_peak_E = E_list[np.where(spectral_flux_BB == spectral_flux_peak)]
print("spectral flux peak = ", spectral_flux_peak)
print("spectral flux peak energy (should be 3 for BB) = ", spectral_flux_peak_E)

plt.title('E/kT vs. Spectral Flux')
plt.plot(E_list,spectral_flux_BB)
plt.xlim(0, 15)
plt.show()