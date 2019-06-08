#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 13:19:09 2019

@author: charlee

This program takes a given value of M/R and a given value of Psi, then interpolates
using the table in lookup_alpha.txt to obtain values of alpha and dmu/dzeta.

Spherical star only. 
"""


import numpy as np
from Interp_g_T_zeta import get_nearest_uneven

M_over_R_val = 0.2007
psi_val = 0.01 

alpha_all, psi_all, dmu_dzeta_all, M_over_R_all =\
    np.loadtxt("lookup_alpha.txt", usecols = (0,2,3,5), unpack = True)

M_over_R_rounded = np.zeros(len(M_over_R_all))

for i in range(len(M_over_R_all)): 
    M_over_R_rounded[i] = np.round(M_over_R_all[i], 6)
#print(psi)
max_M_over_R = np.max(M_over_R_all)
min_M_over_R = np.min(M_over_R_all)

if np.round(M_over_R_val,6) in M_over_R_rounded: 
    print("Case 1")
    print("Exact M/R is in table (within 6 decimal places).")
    MR_filter = np.where(M_over_R_val == M_over_R_all)[0]
    
    alpha = alpha_all[MR_filter]
    dmu_dzeta = dmu_dzeta_all[MR_filter]
    psi = psi_all[MR_filter]
    psi_rounded = np.zeros(len(psi))
    
    for i in range(len(psi)): 
        psi_rounded[i] = np.round(psi[i], 6)
    #print(psi)
    max_psi = np.max(psi)
    min_psi = np.min(psi)
    
    if psi_val in psi:
        print("Exact psi is in table (within 6 decimal places).")
        psi_idx_temp = np.where(np.round(psi_val, 6) == psi_rounded)[0]
        psi_idx = int(psi_idx_temp)
        alpha_final = alpha[psi_idx]
        dmu_dzeta_final= dmu_dzeta[psi_idx]
        psi_final = psi_val
        print(psi_final, alpha_final, dmu_dzeta_final)
        
    elif psi_val > max_psi:
        print("Psi value is larger than largest in table.")
        psi_idx = np.where(np.round(max_psi, 6) == psi_rounded)[0]
        alpha_final = float(alpha[psi_idx])
        dmu_dzeta_final= float(dmu_dzeta[psi_idx])
        psi_final = float(max_psi)
        
    elif psi_val < min_psi: 
        print("Psi value is smaller than smallest in table.")
        psi_idx= np.where(np.round(min_psi, 6) == psi_rounded)[0]
        alpha_final = float(alpha[psi_idx])
        dmu_dzeta_final= float(dmu_dzeta[psi_idx])
        psi_final = float(min_psi)
        
    else: 
        print("Interpolating.")
        psi_high, psi_low = get_nearest_uneven(psi, psi_val)
        psi_high_idx_temp = np.where(psi_rounded == np.round(psi_high,6))[0]
        psi_low_idx_temp = np.where(psi_rounded == np.round(psi_low,6))[0]
        
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
                   
        #Interpolate corresponding dmu_dzeta values
        dmu_dzeta_high = dmu_dzeta[psi_high_idx]
        dmu_dzeta_low = dmu_dzeta[psi_low_idx]
        
        if dmu_dzeta_high != dmu_dzeta_low:
            #dmu_dzeta_final = (dmu_dzeta_high + deta_dmu_low)/2
            dmu_dzeta_final = ((dmu_dzeta_high - dmu_dzeta_low)/(psi_high - psi_low))\
                *(psi_val - psi_low) + dmu_dzeta_low
            #print(deta_dmu_final)
        else: 
            dmu_dzeta_final = dmu_dzeta_high
            #print(deta_dmu_final)
            
    print(alpha_final, dmu_dzeta_final)

elif M_over_R_val > max_M_over_R: 
    print("Case 2")
    print("M/R value is larger than largest in table." )
    MR_filter = np.where(np.round(max_M_over_R, 6) == M_over_R_rounded)[0]
    
    alpha = alpha_all[MR_filter]
    dmu_dzeta = dmu_dzeta_all[MR_filter]
    psi = psi_all[MR_filter]
    psi_rounded = np.zeros(len(psi))
    
    for i in range(len(psi)): 
        psi_rounded[i] = np.round(psi[i], 6)
    #print(psi)
    max_psi = np.max(psi)
    min_psi = np.min(psi)
    
    if psi_val in psi:
        print("Exact psi is in table (within 6 decimal places).")
        psi_idx_temp = np.where(np.round(psi_val, 6) == psi_rounded)[0]
        psi_idx = int(psi_idx_temp)
        alpha_final = alpha[psi_idx][0]
        dmu_dzeta_final= dmu_dzeta[psi_idx][0]
        psi_final = psi_val
        print(psi_final, alpha_final, dmu_dzeta_final)
        
    elif psi_val > max_psi:
        print("Psi value is larger than largest in table.")
        psi_idx = np.where(np.round(max_psi, 6) == psi_rounded)[0]
        alpha_final = alpha[psi_idx]
        dmu_dzeta_final= dmu_dzeta[psi_idx]
        psi_final = max_psi
        
    elif psi_val < min_psi: 
        print("Psi value is smaller than smallest in table.")
        psi_idx= np.where(np.round(min_psi, 6) == psi_rounded)[0]
        alpha_final = alpha[psi_idx]
        dmu_dzeta_final= dmu_dzeta[psi_idx]
        psi_final = min_psi
        
    else: 
        print("Interpolating.")
        psi_high, psi_low = get_nearest_uneven(psi, psi_val)
        psi_high_idx_temp = np.where(psi_rounded == np.round(psi_high,6))[0]
        psi_low_idx_temp = np.where(psi_rounded == np.round(psi_low,6))[0]
        
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
                   
        #Interpolate corresponding dmu_dzeta values
        dmu_dzeta_high = dmu_dzeta[psi_high_idx]
        dmu_dzeta_low = dmu_dzeta[psi_low_idx]
        
        if dmu_dzeta_high != dmu_dzeta_low:
            #dmu_dzeta_final = (dmu_dzeta_high + deta_dmu_low)/2
            dmu_dzeta_final = ((dmu_dzeta_high - dmu_dzeta_low)/(psi_high - psi_low))\
                *(psi_val - psi_low) + dmu_dzeta_low
            #print(deta_dmu_final)
        else: 
            dmu_dzeta_final = dmu_dzeta_high
            #print(deta_dmu_final)
            
    print(alpha_final, dmu_dzeta_final)
    #print(M_over_R_final)
elif M_over_R_val < min_M_over_R: 
    print("Case 3")
    print("M/R value is smaller than smallest in table." )
    
    MR_filter = np.where(np.round(min_M_over_R, 6) == M_over_R_rounded)[0]
        
    alpha = alpha_all[MR_filter]
    dmu_dzeta = dmu_dzeta_all[MR_filter]
    psi = psi_all[MR_filter]
    psi_rounded = np.zeros(len(psi))
    
    for i in range(len(psi)): 
        psi_rounded[i] = np.round(psi[i], 6)
    #print(psi)
    max_psi = np.max(psi)
    min_psi = np.min(psi)
    
    if psi_val in psi:
        print("Exact psi is in table (within 6 decimal places).")
        psi_idx_temp = np.where(np.round(psi_val, 6) == psi_rounded)[0]
        psi_idx = int(psi_idx_temp)
        alpha_final = alpha[psi_idx][0]
        dmu_dzeta_final= dmu_dzeta[psi_idx][0]
        psi_final = psi_val
        print(psi_final, alpha_final, dmu_dzeta_final)
        
    elif psi_val > max_psi:
        print("Psi value is larger than largest in table.")
        psi_idx = np.where(np.round(max_psi, 6) == psi_rounded)[0]
        alpha_final = alpha[psi_idx]
        dmu_dzeta_final= dmu_dzeta[psi_idx]
        psi_final = max_psi
        
    elif psi_val < min_psi: 
        print("Psi value is smaller than smallest in table.")
        psi_idx= np.where(np.round(min_psi, 6) == psi_rounded)[0]
        alpha_final = alpha[psi_idx]
        dmu_dzeta_final= dmu_dzeta[psi_idx]
        psi_final = min_psi
        
    else: 
        print("Interpolating.")
        psi_high, psi_low = get_nearest_uneven(psi, psi_val)
        psi_high_idx_temp = np.where(psi_rounded == np.round(psi_high,6))[0]
        psi_low_idx_temp = np.where(psi_rounded == np.round(psi_low,6))[0]
        
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
                   
        #Interpolate corresponding dmu_dzeta values
        dmu_dzeta_high = dmu_dzeta[psi_high_idx]
        dmu_dzeta_low = dmu_dzeta[psi_low_idx]
        
        if dmu_dzeta_high != dmu_dzeta_low:
            #dmu_dzeta_final = (dmu_dzeta_high + deta_dmu_low)/2
            dmu_dzeta_final = ((dmu_dzeta_high - dmu_dzeta_low)/(psi_high - psi_low))\
                *(psi_val - psi_low) + dmu_dzeta_low
            #print(deta_dmu_final)
        else: 
            dmu_dzeta_final = dmu_dzeta_high
            #print(deta_dmu_final)
            
    print(alpha_final, dmu_dzeta_final)
    
else: 
    print("Interpolating M/R.")
    #Catch erroneous psi values 
    min_psi = min(psi_all)
    max_psi = max(psi_all)
    psi_rounded = np.zeros(len(psi_all))
    
    for i in range(len(psi_all)): 
        psi_rounded[i] = np.round(psi_all[i], 6)
  
    print("Case 4 - General")
    print("Interpolating psi.")
    #Catch erroneous psi values 
    if psi_val < min_psi: 
        print("Psi value smaller than smallest in table.")
        psi_val = min_psi
    elif psi_val > max_psi: 
        print("Psi value larger than largest in table.")
        psi_val = max_psi 
        
    #Find closest pair to given M/R
    M_over_R_high, M_over_R_low = get_nearest_uneven(M_over_R_all, M_over_R_val)
    M_over_R_high_idx = np.where(M_over_R_rounded == np.round(M_over_R_high,6))[0]
    M_over_R_low_idx = np.where(M_over_R_rounded == np.round(M_over_R_low,6))[0] 

    #Find closest pair to given psi
    psi_all_MR_high = psi_all[M_over_R_high_idx]
    psi_all_MR_low = psi_all[M_over_R_low_idx]
    alpha_all_MR_high = alpha_all[M_over_R_high_idx]
    alpha_all_MR_low = alpha_all[M_over_R_low_idx]
    dmu_dzeta_all_MR_high = dmu_dzeta_all[M_over_R_high_idx]
    dmu_dzeta_all_MR_low = dmu_dzeta_all[M_over_R_low_idx]
    
    psi_all_MR_high_rounded = np.zeros(len(psi_all_MR_high))
    psi_all_MR_low_rounded = np.zeros(len(psi_all_MR_low))    
    for i in range(len(psi_all_MR_high)):
        psi_all_MR_high_rounded[i] = np.round(psi_all_MR_high[i],6)
        psi_all_MR_low_rounded[i] = np.round(psi_all_MR_low[i],6)
    
    #Catch erroneous psi values 
    new_max_psi = max(psi_all_MR_high)
    new_min_psi = min(psi_all_MR_high)
    if psi_val < new_min_psi: 
        print("Psi value smaller than smallest in table.")
        psi_val = new_min_psi
    elif psi_val > new_max_psi: 
        print("Psi value larger than largest in table.")
        psi_val = new_max_psi   
        
    psi_high_MR_high, psi_low_MR_high = \
                get_nearest_uneven(psi_all_MR_high, psi_val)
    psi_high_MR_low, psi_low_MR_low = \
                get_nearest_uneven(psi_all_MR_low, psi_val)
    
    psi_high_MR_high_idx = np.where(np.round(psi_high_MR_high, 6) == psi_all_MR_high_rounded)[0]
    psi_low_MR_high_idx = np.where(np.round(psi_low_MR_high, 6) == psi_all_MR_high_rounded)[0]
    psi_high_MR_low_idx = np.where(np.round(psi_high_MR_low, 6) == psi_all_MR_low_rounded)[0]
    psi_low_MR_low_idx = np.where(np.round(psi_low_MR_low, 6) == psi_all_MR_low_rounded)[0]
    
    alpha_psi_high_MR_high = alpha_all_MR_high[psi_high_MR_high_idx]
    alpha_psi_low_MR_high = alpha_all_MR_high[psi_low_MR_high_idx]
    alpha_psi_high_MR_low = alpha_all_MR_low[psi_high_MR_low_idx]
    alpha_psi_low_MR_low = alpha_all_MR_low[psi_low_MR_low_idx]
    
    dmu_dzeta_psi_high_MR_high = dmu_dzeta_all_MR_high[psi_high_MR_high_idx]
    dmu_dzeta_psi_low_MR_high = dmu_dzeta_all_MR_high[psi_low_MR_high_idx]
    dmu_dzeta_psi_high_MR_low = dmu_dzeta_all_MR_low[psi_high_MR_low_idx]
    dmu_dzeta_psi_low_MR_low = dmu_dzeta_all_MR_low[psi_low_MR_low_idx]

    #Interpolate between psi values
    if alpha_psi_high_MR_low != alpha_psi_low_MR_low:
        alpha_MR_low = float(((alpha_psi_high_MR_low - alpha_psi_low_MR_low)\
            /(psi_high_MR_low - psi_low_MR_low))*\
            (psi_val - psi_low_MR_low) + alpha_psi_low_MR_low)
        print("low alpha = ", alpha_MR_low)   
    else: 
        alpha_MR_low = float(alpha_psi_high_MR_low)
        print("low alpha = ", alpha_MR_low) 
    
    if alpha_psi_high_MR_high != alpha_psi_low_MR_high:
        alpha_MR_high = float(((alpha_psi_high_MR_high - alpha_psi_low_MR_high)\
            /(psi_high_MR_high - psi_low_MR_high))*\
            (psi_val - psi_low_MR_high) + alpha_psi_low_MR_high) 

        print("high alpha = ", alpha_MR_high)
    else: 
        alpha_MR_high = float(alpha_psi_high_MR_high)
        print("high alpha = ", alpha_MR_high)
    
    if dmu_dzeta_psi_high_MR_low != dmu_dzeta_psi_low_MR_low: 
        
        dmu_dzeta_MR_low = float(((dmu_dzeta_psi_high_MR_low - dmu_dzeta_psi_low_MR_low)\
            /(psi_high_MR_low - psi_low_MR_low))*\
            (psi_val - psi_low_MR_low) + dmu_dzeta_psi_low_MR_low)
        print("low dmu/dzeta = ", dmu_dzeta_MR_low)      
    else: 
        dmu_dzeta_MR_low = float(dmu_dzeta_psi_high_MR_low)
        print("low dmu/dzeta = ", dmu_dzeta_MR_low) 
        
    if dmu_dzeta_psi_high_MR_high != dmu_dzeta_psi_low_MR_high:
        dmu_dzeta_MR_high = float(((dmu_dzeta_psi_high_MR_high - dmu_dzeta_psi_low_MR_high)\
            /(psi_high_MR_high - psi_low_MR_high))*\
            (psi_val - psi_low_MR_high) + dmu_dzeta_psi_low_MR_high)
        print("high dmu/dzeta = ", dmu_dzeta_MR_high)
    else: 
        dmu_dzeta_MR_high = float(dmu_dzeta_psi_high_MR_high)
        print("high dmu/dzeta = ", dmu_dzeta_MR_high)
        
    #Interpolate bewteen M/R values
    if alpha_MR_high != alpha_MR_low: 
        alpha_final = float(((alpha_MR_high - alpha_MR_low)\
            /(M_over_R_high - M_over_R_low))*\
            (M_over_R_val - M_over_R_low) + alpha_MR_low)
        print("alpha final =", alpha_final)
    else: 
        alpha_final = alpha_MR_high
        print("alpha final =", alpha_final)

    if dmu_dzeta_MR_high != dmu_dzeta_MR_low:
        dmu_dzeta_final = float(((dmu_dzeta_MR_high - dmu_dzeta_MR_low)\
            /(M_over_R_high - M_over_R_low))*\
            (M_over_R_val - M_over_R_low) + dmu_dzeta_MR_low)
        print("dmu/dzeta final =", dmu_dzeta_final)
    else: 
        dmu_dzeta_final = dmu_dzeta_MR_high
        print("dmu/dzeta final =", dmu_dzeta_final)
   
        

