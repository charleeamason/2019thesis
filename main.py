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
import matplotlib.pyplot as plt
import flux_integral_GR as f
#import os
#import shutil

"""
#Create a directory to store files
directory = './Interp_Ho_data/'
if os.path.exists(directory):
    shutil.rmtree(directory)
os.makedirs(directory)
"""
start_time = time.time()
#Enter T value between 5.1 and 6.8
log_T = 5.34 #Kelvin
T = 10**log_T
#Enter g value between 13.7 and 14.7
log_g = 14.36 #g = GM/R^2 * redshift^2

#List of E/kT values to calculate energies for
E_list = np.linspace(0.05,100,166)
#Relevant angles (radians):
#***avoid |i| > 90 (quadrant issues)***
incl_array = [0,30,60,90]
bolo_flux_array = []

for inclination in incl_array:
#inclination = 30 #inclination in degrees (from spin axis)

    bolo_flux_BB, spectral_flux_BB, spectral_I_BB = f.get_flux_BB(inclination, log_T, E_list)
    #bolo_flux_Ho, spectral_flux, spectral_I = f.get_flux_H(inclination, log_T, log_g, E_list)
    bolo_flux_array.append(bolo_flux_BB)
np.savetxt('GRflux_40pts.txt', bolo_flux_array, comments = '#flux for i = 0,30,60,90' )
print("--- %s seconds ---" % (time.time() - start_time))


#print('flux Ho =', bolo_flux_Ho)
print('BB flux =', bolo_flux_array)

#Plot the spectra
plt.figure(1, figsize = (8,6))
plt.title('log T = {0} log g ={1}'.format(log_T,log_g), size = 16)
plt.xlabel(r'$E/kT$ (dimensionless)', size = 16)
plt.ylabel(r'$I_\nu \mathrm{(erg/cm}^2\mathrm{/ster}\mathrm{)}$', size = 16)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
#plt.xlim(0,20)
plt.plot(E_list, spectral_I_BB, label = 'Pure BB')
#plt.plot(E_list, spectral_flux_BB, marker = 'o', label = 'Pure BB flux' )
#plt.plot(E_list, spectral_flux*T**3, label = 'W.Ho' )
plt.legend()
plt.show()