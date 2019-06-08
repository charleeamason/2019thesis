#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 31 15:57:34 2019

@author: charlee
"""

import numpy as np
import matplotlib.pyplot as plt 
theta_all, phi_all, dOmega_sph_all = np.loadtxt('flux-master/angles_sph.txt', skiprows = 1, unpack = True)
dOmega_obl_all = np.loadtxt('flux-master/angles_obl.txt', skiprows = 1, usecols = 2, unpack = True)

i = 9 #value of phi 
j = 9 #value of theta
theta_filter = np.where(phi_all == phi_all[i]) #pick a value of phi
theta = theta_all[theta_filter]
dOmega_sph_theta = dOmega_sph_all[theta_filter]
dOmega_obl_theta = dOmega_obl_all[theta_filter]

plt.figure(1)
plt.plot(theta, dOmega_sph_theta, label = 'sph')
plt.plot(theta, dOmega_obl_theta, label = 'obl')
plt.legend()

phi_filter = np.where(theta_all == theta_all[j])
phi = phi_all[phi_filter]
dOmega_sph_phi = dOmega_sph_all[phi_filter]
dOmega_obl_phi = dOmega_sph_all[phi_filter]

plt.figure(2)
plt.plot(phi, dOmega_sph_phi, label = 'sph')
plt.plot(phi, dOmega_obl_phi, label = 'obl')
plt.legend()
plt.show(1, 2)