#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 14:14:13 2018

@author: charlee

A program which finds the flux emitted by a spherical, nonrotating neutron star
using Newtonian gravity only in an observer-aligned coordinate system. 
A hydrogen atmosphere model from Ho et al. is compared to a simple blackbody.
"""

import numpy as np
import glob
import Interp_g_T_mu as Ho
import scipy.integrate as inte
import time
import matplotlib.pyplot as plt
import os
import shutil

#Create a directory to store files
directory = './Interp_Ho_data/'
if os.path.exists(directory):
    shutil.rmtree(directory)
os.makedirs(directory)
    
start_time = time.time()

#Define constants
h = 6.626e-27 #planck const in erg*s
c = 3e10 #speed of light cm/s 
k = 1.38e-16 #Boltzmann in erg/K

#Initial parameters:
M = 2.78e33 #1.4 solar masses in grams
R = 1.2e6 #cm
dist = 6.1713552e20 #200 pc in cm

#Enter T value between 5.1 and 6.8
log_T = 5.34 #Kelvin
T = 10**log_T
#Enter g value between 13.7 and 14.7
log_g = 14.36 #g = GM/R^2 * redshift^2
                  
#Define range of mu points on star to integrate 
alpha = np.linspace(0, np.pi/2, 166) #size 166 needed if plotting spectra
mu_for_inte = np.linspace(1e-06, 1, len(alpha)) #cos alpha
solid_const = R**2/dist**2
solid_angle = solid_const*(np.sin(alpha)*np.cos(alpha))
E_kT = np.linspace(0.05, 20, len(mu_for_inte)) #E/kT
E_kT_obs = E_kT
"""
Jmu = [] #Function of Inu dependent on mu
#Get H atm model for given temp and g
E, Inu, Inu_T_high, Inu_T_low, mu_array,\
                mu_g_high, mu_g_low = Ho.interp_T_and_g(log_T, log_g) 
                
#print("--- %s seconds ---" % (time.time() - start_time))
for mu in mu_for_inte:
              
    #Modify H model given mu       
    Ho.interp_mu(mu, E, Inu, Inu_T_high, Inu_T_low, mu_array, mu_g_high,\
                    mu_g_low, log_g, log_T)
    #print("--- %s seconds ---" % (time.time() - start_time))
    
#Load files from H atm model
print("Loading data.")
filenames = glob.glob('./Interp_Ho_data/output_*.txt') #2D array E, Inu
#print("--- %s seconds ---" % (time.time() - start_time))
filenames.sort()

#Load data into arrays
E_temp = []
Inu_temp = []
for i in range(len(filenames)):
    E_temp, Inu_temp = np.loadtxt(filenames[i], unpack = True)
    #pick a couple mu values for plotting spectra
    if mu_for_inte[i] == 1.0: 
        Inu_1 = Inu_temp
    elif mu_for_inte[i] == 1e-6:
        Inu_0 = Inu_temp
    #Get Inu_obs
    #Inu_obs = Inu*(E_obs/E)**3

    #Integrate Inu for fixed mu and convert to correct units
    E_inte = inte.simps(Inu_temp, E_temp)*(k*T)*T**3/h 
    Jmu.append(E_inte)
    #print("--- %s seconds ---" % (time.time() - start_time))
    #print('Inu integral =', E_inte) 
    
#Get Ho flux    
Jmu_inte = 2*np.pi*inte.simps(Jmu, alpha)  
solid_angle_inte = inte.simps(solid_angle, alpha)
flux_Ho = Jmu_inte*solid_angle_inte 
print('flux Ho =', flux_Ho)
#print("--- %s seconds ---" % (time.time() - start_time))
"""
#Get BB model
const = 2/(c**2*h**2) 
Inu_BB = const*((k*T*E_kT)**3)/(np.e**(E_kT) - 1) 
Inu_BB_plot = const*((k*E_kT)**3)/(np.e**(E_kT) - 1) #want to plot Inu/T^3 not Inu
Inu_BB_obs = Inu_BB*(E_kT_obs/E_kT)**3

#Get BB flux
print("calculating flux_BB")
const_inte = (k*T)/h #Want to integrate over nu (freq), so need this conv. factor
E_inte_BB = const_inte*inte.simps(Inu_BB, E_kT) #E approx 4e16
#print("E_inte_BB =", E_inte_BB)
alpha_inte_BB = 2*np.pi*inte.simps(solid_angle, alpha) #phi and alpha approx. 1e-29
print(alpha_inte_BB)
#print("alpha_inte_BB =", alpha_inte_BB)
flux_BB = E_inte_BB*alpha_inte_BB
print('flux BB =', flux_BB)
#flux_BB_array.append(flux_BB)

print("--- %s seconds ---" % (time.time() - start_time))
"""
#Plot the spectra
plt.figure(1, figsize = (8,6))
plt.title('log T = {0} log g ={1}'.format(log_T,log_g), size = 16)
plt.xlabel(r'$E/kT$ (dimensionless)', size = 16)
plt.ylabel(r'$I_\nu/T^{3} \mathrm{(erg/cm}^2\mathrm{/ster/K}^3\mathrm{)}$', size = 16)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
#plt.xlim(0,30)
plt.plot(E_kT, Inu_BB_plot, label = 'Pure BB')
plt.plot(E_kT, Inu_0, label = 'W.Ho, mu = 0')
plt.plot(E_kT, Inu_1, label = 'W.Ho, mu = 1')
plt.legend()
plt.show() """