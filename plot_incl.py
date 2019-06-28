#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 15:06:58 2019

@author: charlee
"""

import numpy as np
import matplotlib.pyplot as plt

"""Comparision of inclinations for oblate NS. M/R fixed @ 0.22. Four plots, for each value
of spin (=300 Hz, 600 Hz) and atmosphere type (BB vs H).""" 

ntheta = 30

# spin = 300, bb
npzfile_bb_spin300_incl0 = np.load('./fluxdata/ntheta' + str(ntheta) + '/BB_spectra_obl_spin300_MR22_incl0.npz')
npzfile_bb_spin300_incl30 = np.load('./fluxdata/ntheta' + str(ntheta) + '/BB_spectra_obl_spin300_MR22_incl30.npz')
npzfile_bb_spin300_incl45 = np.load('./fluxdata/ntheta' + str(ntheta) + '/BB_spectra_obl_spin300_MR22_incl45.npz')
npzfile_bb_spin300_incl60 = np.load('./fluxdata/ntheta' + str(ntheta) + '/BB_spectra_obl_spin300_MR22_incl60.npz')
npzfile_bb_spin300_incl90 = np.load('./fluxdata/ntheta' + str(ntheta) + '/BB_spectra_obl_spin300_MR22_incl90.npz')

E_bb_spin300_incl0 = npzfile_bb_spin300_incl0['arr_0']
F_bb_spin300_incl0 = npzfile_bb_spin300_incl0['arr_1']

E_bb_spin300_incl30 = npzfile_bb_spin300_incl30['arr_0']
F_bb_spin300_incl30 = npzfile_bb_spin300_incl30['arr_1']

E_bb_spin300_incl45 = npzfile_bb_spin300_incl45['arr_0']
F_bb_spin300_incl45 = npzfile_bb_spin300_incl45['arr_1']

E_bb_spin300_incl60 = npzfile_bb_spin300_incl60['arr_0']
F_bb_spin300_incl60 = npzfile_bb_spin300_incl60['arr_1']

E_bb_spin300_incl90 = npzfile_bb_spin300_incl90['arr_0']
F_bb_spin300_incl90 = npzfile_bb_spin300_incl90['arr_1']

#spin = 300, H
npzfile_H_spin300_incl0 = np.load('./fluxdata/ntheta' + str(ntheta) + '/H_spectra_obl_spin300_MR22_incl0.npz')
npzfile_H_spin300_incl30 = np.load('./fluxdata/ntheta' + str(ntheta) + '/H_spectra_obl_spin300_MR22_incl30.npz')
npzfile_H_spin300_incl45 = np.load('./fluxdata/ntheta' + str(ntheta) + '/H_spectra_obl_spin300_MR22_incl45.npz')
npzfile_H_spin300_incl60 = np.load('./fluxdata/ntheta' + str(ntheta) + '/H_spectra_obl_spin300_MR22_incl60.npz')
npzfile_H_spin300_incl90 = np.load('./fluxdata/ntheta' + str(ntheta) + '/H_spectra_obl_spin300_MR22_incl90.npz')

E_H_spin300_incl0 = npzfile_H_spin300_incl0['arr_0']
F_H_spin300_incl0 = npzfile_H_spin300_incl0['arr_1']

E_H_spin300_incl30 = npzfile_H_spin300_incl30['arr_0']
F_H_spin300_incl30 = npzfile_H_spin300_incl30['arr_1']

E_H_spin300_incl45 = npzfile_H_spin300_incl45['arr_0']
F_H_spin300_incl45 = npzfile_H_spin300_incl45['arr_1']

E_H_spin300_incl60 = npzfile_H_spin300_incl60['arr_0']
F_H_spin300_incl60 = npzfile_H_spin300_incl60['arr_1']

E_H_spin300_incl90 = npzfile_H_spin300_incl90['arr_0']
F_H_spin300_incl90 = npzfile_H_spin300_incl90['arr_1']

#spin = 600, bb
npzfile_bb_spin600_incl0 = np.load('./fluxdata/ntheta' + str(ntheta) + '/BB_spectra_obl_spin600_MR22_incl0.npz')
npzfile_bb_spin600_incl30 = np.load('./fluxdata/ntheta' + str(ntheta) + '/BB_spectra_obl_spin600_MR22_incl30.npz')
npzfile_bb_spin600_incl45 = np.load('./fluxdata/ntheta' + str(ntheta) + '/BB_spectra_obl_spin600_MR22_incl45.npz')
npzfile_bb_spin600_incl60 = np.load('./fluxdata/ntheta' + str(ntheta) + '/BB_spectra_obl_spin600_MR22_incl60.npz')
npzfile_bb_spin600_incl90 = np.load('./fluxdata/ntheta' + str(ntheta) + '/BB_spectra_obl_spin600_MR22_incl90.npz')

E_bb_spin600_incl0 = npzfile_bb_spin600_incl0['arr_0']
F_bb_spin600_incl0 = npzfile_bb_spin600_incl0['arr_1']

E_bb_spin600_incl30 = npzfile_bb_spin600_incl30['arr_0']
F_bb_spin600_incl30 = npzfile_bb_spin600_incl30['arr_1']

E_bb_spin600_incl45 = npzfile_bb_spin600_incl45['arr_0']
F_bb_spin600_incl45 = npzfile_bb_spin600_incl45['arr_1']

E_bb_spin600_incl60 = npzfile_bb_spin600_incl60['arr_0']
F_bb_spin600_incl60 = npzfile_bb_spin600_incl60['arr_1']

E_bb_spin600_incl90 = npzfile_bb_spin600_incl90['arr_0']
F_bb_spin600_incl90 = npzfile_bb_spin600_incl90['arr_1']

#spin = 600, H
npzfile_H_spin600_incl0 = np.load('./fluxdata/ntheta' + str(ntheta) + '/H_spectra_obl_spin600_MR22_incl0.npz')
npzfile_H_spin600_incl30 = np.load('./fluxdata/ntheta' + str(ntheta) + '/H_spectra_obl_spin600_MR22_incl30.npz')
npzfile_H_spin600_incl45 = np.load('./fluxdata/ntheta' + str(ntheta) + '/H_spectra_obl_spin600_MR22_incl45.npz')
npzfile_H_spin600_incl60 = np.load('./fluxdata/ntheta' + str(ntheta) + '/H_spectra_obl_spin600_MR22_incl60.npz')
npzfile_H_spin600_incl90 = np.load('./fluxdata/ntheta' + str(ntheta) + '/H_spectra_obl_spin600_MR22_incl90.npz')

E_H_spin600_incl0 = npzfile_H_spin600_incl0['arr_0']
F_H_spin600_incl0 = npzfile_H_spin600_incl0['arr_1']

E_H_spin600_incl30 = npzfile_H_spin600_incl30['arr_0']
F_H_spin600_incl30 = npzfile_H_spin600_incl30['arr_1']

E_H_spin600_incl45 = npzfile_H_spin600_incl45['arr_0']
F_H_spin600_incl45 = npzfile_H_spin600_incl45['arr_1']

E_H_spin600_incl60 = npzfile_H_spin600_incl60['arr_0']
F_H_spin600_incl60 = npzfile_H_spin600_incl60['arr_1']

E_H_spin600_incl90 = npzfile_H_spin600_incl90['arr_0']
F_H_spin600_incl90 = npzfile_H_spin600_incl90['arr_1']

#spin = 300, bb
plt.figure(1, figsize = (8,6))
#plt.title('spin = 300, bb', size = 36)
plt.xlabel(r'$E/kT$ [dimensionless]', size = 30)
plt.ylabel(r'Spectral Flux [$\mathrm{10^{-12}\ erg/cm}^2\mathrm{/s/ster}$]', size = 30)
#plt.ylabel(r'$I_\nu \mathrm{(erg/cm}^2\mathrm{/ster}\mathrm{)}$', size = 16)
plt.xticks(fontsize = 30)
plt.yticks(fontsize = 30)
plt.xlim(0,15)
plt.ylim(0, 7)
plt.plot(E_bb_spin300_incl0,F_bb_spin300_incl0/1e-12, label = 'incl = 0',linewidth = 3)
plt.plot(E_bb_spin300_incl30,F_bb_spin300_incl30/1e-12, label = 'incl = 30',linewidth = 3)
plt.plot(E_bb_spin300_incl45,F_bb_spin300_incl45/1e-12, label = 'incl = 45',linewidth = 3)
plt.plot(E_bb_spin300_incl60,F_bb_spin300_incl60/1e-12, label = 'incl = 60',linewidth = 3)
plt.plot(E_bb_spin300_incl90,F_bb_spin300_incl90/1e-12, label = 'incl = 90',linewidth = 3)
plt.legend(prop={'size': 30})
plt.show()
#plt.savefig('./fluxdata/ntheta' + str(ntheta) + '/incl_compare_bb_spin300.png')

#spin = 300, H
plt.figure(2, figsize = (8,6))
#plt.title('spin = 300, H', size = 36)
plt.xlabel(r'$E/kT$ [dimensionless]', size = 30)
plt.ylabel(r'Spectral Flux [$\mathrm{10^{-12}\ erg/cm}^2\mathrm{/s/ster}$]', size = 30)
#plt.ylabel(r'$I_\nu \mathrm{(erg/cm}^2\mathrm{/ster}\mathrm{)}$', size = 16)
plt.xticks(fontsize = 30)
plt.yticks(fontsize = 30)
plt.xlim(0,15)
plt.ylim(0, 7)
plt.plot(E_H_spin300_incl0,F_H_spin300_incl0/1e-12, label = 'incl = 0',linewidth = 3)
plt.plot(E_H_spin300_incl30,F_H_spin300_incl30/1e-12, label = 'incl = 30',linewidth = 3)
plt.plot(E_H_spin300_incl45,F_H_spin300_incl45/1e-12, label = 'incl = 45',linewidth = 3)
plt.plot(E_H_spin300_incl60,F_H_spin300_incl60/1e-12, label = 'incl = 60',linewidth = 3)
plt.plot(E_H_spin300_incl90,F_H_spin300_incl90/1e-12, label = 'incl = 90',linewidth = 3)
plt.legend(prop={'size': 30})
plt.show()
#plt.savefig('./fluxdata/ntheta' + str(ntheta) + '/incl_compare_H_spin300.png')

#spin = 600, bb
plt.figure(3, figsize = (8,6))
#plt.title('spin = 600, bb', size = 36)
plt.xlabel(r'$E/kT$ [dimensionless]', size = 30)
plt.ylabel(r'Spectral Flux [$\mathrm{10^{-12}\ erg/cm}^2\mathrm{/s/ster}$]', size = 30)
#plt.ylabel(r'$I_\nu \mathrm{(erg/cm}^2\mathrm{/ster}\mathrm{)}$', size = 16)
plt.xticks(fontsize = 30)
plt.yticks(fontsize = 30)
plt.xlim(0,15)
plt.ylim(0, 7)
plt.plot(E_bb_spin600_incl0,F_bb_spin600_incl0/1e-12, label = 'incl = 0',linewidth = 3)
plt.plot(E_bb_spin600_incl30,F_bb_spin600_incl30/1e-12, label = 'incl = 30',linewidth = 3)
plt.plot(E_bb_spin600_incl45,F_bb_spin600_incl45/1e-12, label = 'incl = 45',linewidth = 3)
plt.plot(E_bb_spin600_incl60,F_bb_spin600_incl60/1e-12, label = 'incl = 60',linewidth = 3)
plt.plot(E_bb_spin600_incl90,F_bb_spin600_incl90/1e-12, label = 'incl = 90',linewidth = 3)
plt.legend(prop={'size': 30})
plt.show()
#plt.savefig('./fluxdata/ntheta' + str(ntheta) + '/incl_compare_bb_spin600.png')

#spin = 600, H
plt.figure(4, figsize = (8,6))
#plt.title('spin = 600, H', size = 36)
plt.xlabel(r'$E/kT$ [dimensionless]', size = 30)
plt.ylabel(r'Spectral Flux [$\mathrm{10^{-12}\ erg/cm}^2\mathrm{/s/ster}$]', size = 30)
#plt.ylabel(r'$I_\nu \mathrm{(erg/cm}^2\mathrm{/ster}\mathrm{)}$', size = 16)
plt.xticks(fontsize = 30)
plt.yticks(fontsize = 30)
plt.xlim(0,15)
plt.ylim(0, 7)
plt.plot(E_H_spin600_incl0,F_H_spin600_incl0/1e-12, label = 'incl = 0',linewidth = 3)
plt.plot(E_H_spin600_incl30,F_H_spin600_incl30/1e-12, label = 'incl = 30',linewidth = 3)
plt.plot(E_H_spin600_incl45,F_H_spin600_incl45/1e-12, label = 'incl = 45',linewidth = 3)
plt.plot(E_H_spin600_incl60,F_H_spin600_incl60/1e-12, label = 'incl = 60',linewidth = 3)
plt.plot(E_H_spin600_incl90,F_H_spin600_incl90/1e-12, label = 'incl = 90',linewidth = 3)
plt.legend(prop={'size': 30})
plt.show()
#plt.savefig('./fluxdata/ntheta' + str(ntheta) + '/incl_compare_H_spin600.png')


