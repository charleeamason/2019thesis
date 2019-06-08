#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 12:17:41 2019

@author: charlee
"""
import numpy as np
import matplotlib.pyplot as plt 
"""Inclination = 90, M = 2.78e33g, R = 12e5 cm all fixed"""

npzfile_bb = np.load("./fluxdata/SphfluxBB_30pts_spectra_incl90.npz")
npzfile_bb_newt = np.load("./fluxdata/NewtfluxBB_30pts_spectra_incl90.npz")
npzfile_bb_boost = np.load("./fluxdata/SphfluxBB_boost_30pts_spectra_incl90.npz")
npzfile_bb_obl_boost = np.load("./fluxdata/OblfluxBB_30pts_spectra_incl90.npz")

#npzfile_ho = np.load("./fluxdata/SphfluxHo_30pts_spectra_incl90.npz")
npzfile_ho_newt = np.load("./fluxdata/NewtfluxHo_30pts_spectra_incl90.npz")
#npzfile_ho_boost = np.load("./fluxdata/SphfluxHo_boost_30pts_spectra_incl90.npz")
#npzfile_ho_obl_boost = np.load("./fluxdata/OblfluxHo_30pts_spectra_incl90.npz")

E_bb_newt = npzfile_bb_newt['arr_0']
F_bb_newt = npzfile_bb_newt['arr_1']

E_bb = npzfile_bb['arr_0']
F_bb = npzfile_bb['arr_1']

E_bb_boost = npzfile_bb_boost['arr_0']
F_bb_boost = npzfile_bb_boost['arr_1']

E_bb_obl_boost = npzfile_bb_obl_boost['arr_0']
F_bb_obl_boost = npzfile_bb_obl_boost['arr_1']


E_H_newt = npzfile_ho_newt['arr_0']
F_H_newt = npzfile_ho_newt['arr_1']

#E_H_boost = npzfile_ho_boost['arr_0']
#F_H_boost = npzfile_ho_boost['arr_1']

#E_H_obl_boost = npzfile_ho_obl_boost['arr_0']
#F_H_obl_boost = npzfile_ho_obl_boost['arr_1']

#Plot the spectra
plt.figure(1, figsize = (8,6))

plt.title('Blackbody in GR', size = 16)
plt.xlabel(r'$E/kT$ [dimensionless]', size = 16)
plt.ylabel(r'Spectral Flux [$\mathrm{erg/cm}^2\mathrm{/s/ster}$]', size = 16)
#plt.ylabel(r'$I_\nu \mathrm{(erg/cm}^2\mathrm{/ster}\mathrm{)}$', size = 16)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.xlim(0,15)
plt.plot(E_bb,F_bb, label = 'No Doppler boost', color = 'black')
plt.plot(E_bb_boost, F_bb_boost, label = 'With Doppler boost', color = 'cornflowerblue')
plt.legend()

plt.figure(2, figsize = (8,6))
plt.title('Hydrogen vs. Blackbody', size = 16)
plt.xlabel(r'$E/kT$ [dimensionless]', size = 16)
plt.ylabel(r'Spectral Flux [$\mathrm{erg/cm}^2\mathrm{/s/ster}$]', size = 16)
#plt.ylabel(r'$I_\nu \mathrm{(erg/cm}^2\mathrm{/ster}\mathrm{)}$', size = 16)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.xlim(0,30)
plt.plot(E_bb_newt,F_bb_newt, label = 'Blackbody', color = 'black')
plt.plot(E_H_newt, F_H_newt, label = 'Hydrogen Atmosphere', color = 'red')
plt.legend()
    
plt.figure(3, figsize = (8,6))
plt.title('Spherical vs. Oblate Blackbody', size = 16)
plt.xlabel(r'$E/kT$ [dimensionless]', size = 16)
plt.ylabel(r'Spectral Flux [$\mathrm{erg/cm}^2\mathrm{/s/ster}$]', size = 16)
#plt.ylabel(r'$I_\nu \mathrm{(erg/cm}^2\mathrm{/ster}\mathrm{)}$', size = 16)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.xlim(0,15)
plt.plot(E_bb,F_bb, label = 'Spherical No Doppler boost', color = 'cyan')
plt.plot(E_bb_boost,F_bb_boost, label = 'Spherical', color = 'black')
plt.plot(E_bb_obl_boost, F_bb_obl_boost, label = 'Oblate', color = 'green')
plt.legend()

plt.figure(4, figsize = (8,6))
plt.title('Newton vs. GR (Stationary)', size = 16)
plt.xlabel(r'$E/kT$ [dimensionless]', size = 16)
plt.ylabel(r'Spectral Flux [$\mathrm{erg/cm}^2\mathrm{/s/ster}$]', size = 16)
#plt.ylabel(r'$I_\nu \mathrm{(erg/cm}^2\mathrm{/ster}\mathrm{)}$', size = 16)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.xlim(0,15)
plt.plot(E_bb_newt,F_bb_newt, label = 'Newtonian Gravity', color = 'black')
plt.plot(E_bb, F_bb, label = 'GR', color = 'red')
plt.legend()
