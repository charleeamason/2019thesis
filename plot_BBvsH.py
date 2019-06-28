#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 15:24:55 2019

@author: charlee
"""
import numpy as np
import matplotlib.pyplot as plt

"""Plot comparing BB and H atmospheres For a spherical NS. Fixed M/R = 0.22, incl = 45, spin = 0. """

ntheta = 30

npzfile_bb = np.load('./fluxdata/ntheta' + str(ntheta) + '/BB_spectra_sph_spin0_MR22_incl45.npz')
npzfile_H = np.load('./fluxdata/ntheta' + str(ntheta) + '/H_spectra_sph_spin0_MR22_incl45.npz')

E_bb = npzfile_bb['arr_0']
F_bb = npzfile_bb['arr_1']

E_H = npzfile_H['arr_0']
F_H = npzfile_H['arr_1']

plt.figure(1, figsize = (8,6))
#plt.title('Blackbody vs. Hydrogen Atmosphere', size = 36)
plt.xlabel(r'$E/kT$ [dimensionless]', size = 30)
plt.ylabel(r'Spectral Flux [$\mathrm{10^{-12}\ erg/cm}^2\mathrm{/s/ster}$]', size = 30)
#plt.ylabel(r'$I_\nu \mathrm{(erg/cm}^2\mathrm{/ster}\mathrm{)}$', size = 16)
plt.xticks(fontsize = 30)
plt.yticks(fontsize = 30)
plt.xlim(0,15)
plt.ylim(0, 7)
plt.plot(E_bb,F_bb/1e-12, label = 'Blackbody',linewidth = 3)
plt.plot(E_H,F_H/1e-12, label = 'Hydrogen Atmosphere',linewidth = 3)
plt.legend(prop={'size': 30})
plt.show()
#plt.savefig('./fluxdata/ntheta' + str(ntheta) + '/incl_compare_BBvsH.png')


