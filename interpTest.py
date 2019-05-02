#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 16:26:33 2019

@author: charlee
"""
import numpy as np
import glob 
import sys
from Interp_g_T_zeta import get_nearest_pair

log_T= 5.3
log_g = 14.39 

#Get filenames
log_T_name = './logT_files/logT_' + str(int(round(log_T*100))) + '.txt'
filenames = glob.glob('./logT_files/logT_*')

#Convert filenames into numbers for comparision to log_T
filenumber_temp = [file[-7:] for file in filenames]
filenumber = [int(file[:3])/100 for file in filenumber_temp]
filenumber.sort()

#Pick correct files for given log_T
for i in range(len(filenumber)): 
    if log_T > max(filenumber) or log_T < min(filenumber):
        print("log_T is out of range.")
        print("Terminating program.")
        sys.exit()
    elif log_T in filenumber:
        #"high" and "low" will be identical in this case
        print("No temperature interpolation necessary.")
        data_high = np.loadtxt(log_T_name)
        data_low = np.loadtxt(log_T_name)
        break
    else:    
        name_1 = int(round(get_nearest_pair(filenumber,0.05,log_T)[0]*100, 0))
        name_2 = int(round(get_nearest_pair(filenumber,0.05,log_T)[1]*100, 0))
        
        if name_1 > name_2: 
            data_high = np.loadtxt('./logT_files/logT_' + str(name_1) + '.txt')
            data_low = np.loadtxt('./logT_files/logT_' + str(name_2) + '.txt')
        else:
            data_high = np.loadtxt('./logT_files/logT_' + str(name_2) + '.txt')
            data_low = np.loadtxt('./logT_files/logT_' + str(name_1) + '.txt')
            
#Get log_g values from files 
log_g_array = data_high[:,4] #g identical in data_high and data_low

#Catch erroneous log_g values 
if log_g > max(log_g_array) or log_g < min(log_g_array):
    print("log_g is out of range.")
    print("Terminating program.")
    sys.exit()
else:
    #Find closest pair to given log_g
    log_g_high, log_g_low = get_nearest_pair(log_g_array,0.1,log_g)
       
#E and zeta same across all files    
log_E_over_kT = data_low[:,0]  #emitted photon energy/kT
zeta_array = data_low[:,1] #cos theta; theta = zenith angle

#Get specific intensity and temperatures from files
log_Inu_over_T_low = data_low[:,2] #specific intensity/T^3
log_Inu_over_T_high = data_high[:,2] #specific intensity/T^3

log_T_low = data_low[:,3][0] #temp in Kelvin
log_T_high = data_high[:,3][0] #temp in Kelvin

#convert to linear space
E_over_kT = 10**log_E_over_kT
Inu_over_T_low = 10**log_Inu_over_T_low
Inu_over_T_high = 10**log_Inu_over_T_high

#define fixed g ranges
fixed_g_high = np.where(log_g_array == log_g_high)
Inu_g_high_high = Inu_over_T_high[fixed_g_high] 
Inu_g_high_low = Inu_over_T_low[fixed_g_high] 
zeta_g_high = zeta_array[fixed_g_high]

fixed_g_low = np.where(log_g_array == log_g_low)
Inu_g_low_high = Inu_over_T_high[fixed_g_low] 
Inu_g_low_low = Inu_over_T_low[fixed_g_low] 
zeta_g_low = zeta_array[fixed_g_low]
    
E = E_over_kT[fixed_g_low]
if log_T_high != log_T_low: #T value not in array
    if log_g_high != log_g_low:
        #Interpolate 
        print("Interpolating T for log_g_high and log_g_low.")
        Inu_g_high = [None]*len(E)
        for i in range(len(E)):  
            Inu_g_high[i] = ((Inu_g_high_high[i] - Inu_g_high_low[i])/(log_T_high - log_T_low))\
            *(log_T - log_T_low) + Inu_g_high_low[i]  
            
        Inu_g_low = [None]*len(E)
        for i in range(len(E)):  
            Inu_g_low[i] = ((Inu_g_low_high[i] - Inu_g_low_low[i])/(log_T_high - log_T_low))\
            *(log_T - log_T_low) + Inu_g_low_low[i]
            
        print("Interpolating g.")
        Inu = [None]*len(E)
        for i in range(len(E)):  
           # Inu[i] = ((Inu_g_high[i] - Inu_g_low[i])/(log_T_high - log_T_low))\
           # *(log_T - log_T_low) + Inu_g_low[i]
             Inu[i] = ((Inu_g_high[i] - Inu_g_low[i])/(log_g_high - log_g_low))\
            *(log_g - log_g_low) + Inu_g_low[i]
            
        #Define upper and lower Inu     
        Inu_T_high = Inu_over_T_high[fixed_g_high]
        Inu_T_low = Inu_over_T_low[fixed_g_low]
        
    else: 
        print("No g interpolation necessary.")
        Inu = [None]*len(E)
        for i in range(len(E)):
            Inu[i] = ((Inu_over_T_high[i] - Inu_over_T_low[i])/(log_T_high - log_T_low))\
            *(log_T - log_T_low) + Inu_over_T_low[i]  
        #Define upper and lower Inu  
        Inu_T_high = Inu_over_T_high
        Inu_T_low = Inu_over_T_low
        
else: #T value is in array!
    if log_g_high != log_g_low:
        #Interpolate   
        print("Interpolating g.")
        Inu = [None]*len(E)
        for i in range(len(E)):  
           # Inu[i] = ((Inu_g_high[i] - Inu_g_low[i])/(log_T_high - log_T_low))\
           # *(log_T - log_T_low) + Inu_g_low[i]
            Inu[i] = ((Inu_g_high_high[i] - Inu_g_low_high[i])/(log_g_high - log_g_low))\
            *(log_g - log_g_low) + Inu_g_low_high[i]
            
        #Define upper and lower Inu     
        Inu_T_high = Inu_over_T_high[fixed_g_high]
        Inu_T_low = Inu_over_T_low[fixed_g_low]
        
    else: 
        print("No g interpolation necessary.")
        Inu = [None]*len(E)
        for i in range(len(E)):
            Inu[i] = ((Inu_over_T_high[i] - Inu_over_T_low[i])/(log_T_high - log_T_low))\
            *(log_T - log_T_low) + Inu_over_T_low[i]  
        #Define upper and lower Inu  
        Inu_T_high = Inu_over_T_high
        Inu_T_low = Inu_over_T_low
        
print(E) #Inu,Inu_T_high, Inu_T_low zeta_array, zeta_g_high,zeta_g_low)