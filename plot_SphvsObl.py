#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 15:29:33 2019

@author: charlee
"""

import numpy as np
import matplotlib.pyplot as plt

"""Plot for comparing a spherical vs. an oblate NS. Blackbody. Fixed M/R = 0.22, spin = 600. All inclinations."""
ntheta = 30 

npzfile_obl_incl0 = np.load('./fluxdata/ntheta' + str(ntheta) + '/BB_spectra_obl_spin600_MR22_incl0.npz')
npzfile_obl_incl30 = np.load('./fluxdata/ntheta' + str(ntheta) + '/BB_spectra_obl_spin600_MR22_incl30.npz')
npzfile_obl_incl45 = np.load('./fluxdata/ntheta' + str(ntheta) + '/BB_spectra_obl_spin600_MR22_incl45.npz')
npzfile_obl_incl60 = np.load('./fluxdata/ntheta' + str(ntheta) + '/BB_spectra_obl_spin600_MR22_incl60.npz')
npzfile_obl_incl90 = np.load('./fluxdata/ntheta' + str(ntheta) + '/BB_spectra_obl_spin600_MR22_incl90.npz')

E_obl_incl0 = npzfile_obl_incl0['arr_0']
F_obl_incl0 = npzfile_obl_incl0['arr_1']

E_obl_incl30 = npzfile_obl_incl30['arr_0']
F_obl_incl30 = npzfile_obl_incl30['arr_1']

E_obl_incl45 = npzfile_obl_incl45['arr_0']
F_obl_incl45 = npzfile_obl_incl45['arr_1']

E_obl_incl60 = npzfile_obl_incl60['arr_0']
F_obl_incl60 = npzfile_obl_incl60['arr_1']

E_obl_incl90 = npzfile_obl_incl90['arr_0']
F_obl_incl90 = npzfile_obl_incl90['arr_1']


npzfile_sph_incl0 = np.load('./fluxdata/ntheta' + str(ntheta) + '/BB_spectra_sph_spin600_MR22_incl0.npz')
npzfile_sph_incl30 = np.load('./fluxdata/ntheta' + str(ntheta) + '/BB_spectra_sph_spin600_MR22_incl30.npz')
npzfile_sph_incl45 = np.load('./fluxdata/ntheta' + str(ntheta) + '/BB_spectra_sph_spin600_MR22_incl45.npz')
npzfile_sph_incl60 = np.load('./fluxdata/ntheta' + str(ntheta) + '/BB_spectra_sph_spin600_MR22_incl60.npz')
npzfile_sph_incl90 = np.load('./fluxdata/ntheta' + str(ntheta) + '/BB_spectra_sph_spin600_MR22_incl90.npz')

E_sph_incl0 = npzfile_sph_incl0['arr_0']
F_sph_incl0 = npzfile_sph_incl0['arr_1']

E_sph_incl30 = npzfile_sph_incl30['arr_0']
F_sph_incl30 = npzfile_sph_incl30['arr_1']

E_sph_incl45 = npzfile_sph_incl45['arr_0']
F_sph_incl45 = npzfile_sph_incl45['arr_1']

E_sph_incl60 = npzfile_sph_incl60['arr_0']
F_sph_incl60 = npzfile_sph_incl60['arr_1']

E_sph_incl90 = npzfile_sph_incl90['arr_0']
F_sph_incl90 = npzfile_sph_incl90['arr_1']

plt.figure(1, figsize = (8,6))
#plt.title('spin = 300, bb', size = 36)
plt.xlabel(r'$E/kT$ [dimensionless]', size = 30)
plt.ylabel(r'Spectral Flux [$\mathrm{10^{-12}\ erg/cm}^2\mathrm{/s/ster}$]', size = 30)
#plt.ylabel(r'$I_\nu \mathrm{(erg/cm}^2\mathrm{/ster}\mathrm{)}$', size = 16)
plt.xticks(fontsize = 30)
plt.yticks(fontsize = 30)
plt.xlim(0,15)
plt.ylim(0, 7)
plt.plot(E_obl_incl0,F_obl_incl0/1e-12, label = 'obl, incl = 0',linewidth = 3)
plt.plot(E_obl_incl30,F_obl_incl30/1e-12, label = 'obl, incl = 30',linewidth = 3)
plt.plot(E_obl_incl45,F_obl_incl45/1e-12, label = 'obl, incl = 45',linewidth = 3)
plt.plot(E_obl_incl60,F_obl_incl60/1e-12, label = 'obl, incl = 60',linewidth = 3)
plt.plot(E_obl_incl90,F_obl_incl90/1e-12, label = 'obl, incl = 90',linewidth = 3)

plt.plot(E_sph_incl0,F_sph_incl0/1e-12, label = 'sph, incl = 0',linewidth = 3)
plt.plot(E_sph_incl30,F_sph_incl30/1e-12, label = 'sph, incl = 30',linewidth = 3)
plt.plot(E_sph_incl45,F_sph_incl45/1e-12, label = 'sph, incl = 45',linewidth = 3)
plt.plot(E_sph_incl60,F_sph_incl60/1e-12, label = 'sph, incl = 60',linewidth = 3)
plt.plot(E_sph_incl90,F_sph_incl90/1e-12, label = 'sph, incl = 90',linewidth = 3)
plt.legend(prop={'size': 30})
plt.show()
#plt.savefig('./fluxdata/ntheta' + str(ntheta) + '/incl_compare_SphvsObl.png')