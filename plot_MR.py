#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 14:52:44 2019

@author: charlee
"""
import numpy as np
import matplotlib.pyplot as plt

"""Plot to compare 3 M/R values. Inclination fixed @ 45 deg. Spin fixed @ 300 Hz."""

ntheta = 30

npzfile_MR17 = np.load('./fluxdata/ntheta' + str(ntheta) + '/BB_spectra_sph_spin300_MR17_incl45.npz')
npzfile_MR22 = np.load('./fluxdata/ntheta' + str(ntheta) + '/BB_spectra_sph_spin300_MR22_incl45.npz')
npzfile_MR27 = np.load('./fluxdata/ntheta' + str(ntheta) + '/BB_spectra_sph_spin300_MR27_incl45.npz')

E_MR17 = npzfile_MR17['arr_0']
F_MR17 = npzfile_MR17['arr_1']

E_MR22 = npzfile_MR22['arr_0']
F_MR22 = npzfile_MR22['arr_1']

E_MR27 = npzfile_MR27['arr_0']
F_MR27 = npzfile_MR27['arr_1']

plt.figure(1, figsize = (8,6))
#plt.title('M/R values', size = 36)
plt.xlabel(r'$E/kT$ [dimensionless]', size = 30)
plt.ylabel(r'Spectral Flux [$\mathrm{10^{-12}\ erg/cm}^2\mathrm{/s/ster}$]', size = 30)
#plt.ylabel(r'$I_\nu \mathrm{(erg/cm}^2\mathrm{/ster}\mathrm{)}$', size = 16)
plt.xticks(fontsize = 30)
plt.yticks(fontsize = 30)
plt.xlim(0,15)
plt.ylim(0, 7)
plt.plot(E_MR17,F_MR17/1e-12, label = 'M/R = 0.17', color = 'black',\
             linewidth = 3)
plt.plot(E_MR22,F_MR22/1e-12, label = 'M/R = 0.22', color = 'red',\
             linewidth = 3)
plt.plot(E_MR27,F_MR27/1e-12, label = 'M/R = 0.27', color = 'blue',\
             linewidth = 3)
plt.legend(prop={'size': 30})
plt.show()
#plt.savefig('./fluxdata/ntheta' + str(ntheta) + '/MR_comparison.png')