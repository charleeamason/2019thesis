#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 15:39:26 2019

@author: charlee
"""

import numpy as np
import matplotlib.pyplot as plt

"""Plot for comparing a GR blackbody vs a Newtonian blackbody. Spherical NS. Fixed M/R = 0.22, incl = 45, spin = 0."""
ntheta = 30 

npzfile_GR = np.load('./fluxdata/ntheta' + str(ntheta) + '/BB_spectra_sph_spin0_MR22_incl45.npz')
npzfile_newt =  np.load('./fluxdata/ntheta' + str(ntheta) + '/BB_spectra_newt_MR22_incl45.npz')

E_GR = npzfile_GR['arr_0']
F_GR = npzfile_GR['arr_1']

E_newt = npzfile_newt['arr_0']
F_newt = npzfile_newt['arr_1']

plt.figure(1, figsize = (8,6))
#plt.title('GR vs. Newtonian Gravity', size = 36)
plt.xlabel(r'$E/kT$ [dimensionless]', size = 30)
plt.ylabel(r'Spectral Flux [$\mathrm{10^{-12}\ erg/cm}^2\mathrm{/s/ster}$]', size = 30)
#plt.ylabel(r'$I_\nu \mathrm{(erg/cm}^2\mathrm{/ster}\mathrm{)}$', size = 16)
plt.xticks(fontsize = 30)
plt.yticks(fontsize = 30)
plt.xlim(0,15)
plt.ylim(0, 7)
plt.plot(E_GR,F_GR/1e-12, label = 'General Relativity',linewidth = 3)
plt.plot(E_newt,F_newt/1e-12, label = 'Newtonian Gravity',linewidth = 3)
plt.legend(prop={'size': 30})
plt.show()
#plt.savefig('./fluxdata/ntheta' + str(ntheta) + '/incl_compare_GRvsNewt.png')