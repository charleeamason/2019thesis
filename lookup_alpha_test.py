#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  1 14:30:06 2019

@author: charlee
"""
  
psi_val = 2.0943954357265286
M_over_R_val = 0.1

import numpy as np
from Interp_g_T_zeta import get_nearest_uneven

   # print(psi_val)
alpha_all, psi_all, dzeta_dmu_all, M_over_R_all =\
    np.loadtxt("lookup_alpha.txt", usecols = (0,2,3,5), unpack = True)
MR_filter = np.where(M_over_R_val == M_over_R_all)[0]

alpha = alpha_all[MR_filter]
dzeta_dmu = dzeta_dmu_all[MR_filter]
psi = psi_all[MR_filter]
psi_rounded = np.zeros(len(psi))

for i in range(len(psi)): 
    psi_rounded[i] = np.round(psi[i], 4)
#print(psi)
max_psi = np.max(psi)
min_psi = np.min(psi)

if psi_val in psi:
    print("Exact psi is in table (within 4 decimal places).")
    psi_idx_temp = np.where(np.round(psi_val, 4) == psi)[0]
    psi_idx = int(psi_idx_temp)
    alpha_final = alpha[psi_idx][0]
    dzeta_dmu_final= dzeta_dmu[psi_idx][0]
    psi_final = psi_val
    print(psi_final, alpha_final, dzeta_dmu_final)
    
elif psi_val > max_psi:
    print("Psi value is larger than largest in table.")
    psi_idx_temp = np.where(np.round(max_psi, 4) == psi_rounded)[0]
    psi_idx = int(psi_idx_temp)
    print(max_psi)
    print(psi_idx)
    alpha_final = alpha[psi_idx]
    dzeta_dmu_final= dzeta_dmu[psi_idx]
    psi_final = max_psi
    
elif psi_val < min_psi: 
    print("Psi value is smaller than smallest in table.")
    psi_idx= np.where(np.round(min_psi, 4) == psi)[0]
    alpha_final = alpha[psi_idx]
    dzeta_dmu_final= dzeta_dmu[psi_idx]
    psi_final = min_psi
    
else: 
    print("Interpolating.")
    psi_high, psi_low = get_nearest_uneven(psi, psi_val)
    psi_high_idx_temp = np.where(psi_rounded == np.round(psi_high,4))[0]
    psi_low_idx_temp = np.where(psi_rounded == np.round(psi_low,4))[0]
    
    psi_high_idx = int(psi_high_idx_temp)
    psi_low_idx = int(psi_low_idx_temp)
    
    #Interpolate corresponding alpha values
    alpha_high = alpha[psi_high_idx]
    alpha_low = alpha[psi_low_idx]
    
    if alpha_high != alpha_low:
       # alpha_final = (alpha_high + alpha_low)/2
        alpha_final = ((alpha_high - alpha_low)/(psi_high - psi_low))\
            *(psi_val - psi_low) + alpha_low
       # print(alpha_final)
    else: 
        alpha_final = alpha_high
        #print(alpha_final) 
               
    #Interpolate corresponding dzeta_dmu values
    dzeta_dmu_high = dzeta_dmu[psi_high_idx]
    dzeta_dmu_low = dzeta_dmu[psi_low_idx]
    
    if dzeta_dmu_high != dzeta_dmu_low:
        #dzeta_dmu_final = (dzeta_dmu_high + deta_dmu_low)/2
        dzeta_dmu_final = ((dzeta_dmu_high - dzeta_dmu_low)/(psi_high - psi_low))\
            *(psi_val - psi_low) + dzeta_dmu_low
        #print(deta_dmu_final)
    else: 
        dzeta_dmu_final = dzeta_dmu_high
        #print(deta_dmu_final)
        
print(alpha_final, dzeta_dmu_final)