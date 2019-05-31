#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 31 12:10:42 2019

@author: charlee
"""
#import Interp_g_T_zeta as Ho
from Interp_g_T_zeta import get_nearest_uneven
from Interp_g_T_zeta import get_nearest_pair
import sys
import glob
import numpy as np

log_T = 5.3
log_g = 14.27
#log_g = 14.03
zeta = 0.163541
#Get H atm model for given temp and g
#E, Inu, Inu_T_high, Inu_T_low, zeta_array,\
#                zeta_g_high, zeta_g_low = Ho.interp_T_and_g(log_T, log_g)

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
    log_g_high = np.round(log_g_high,2)
    log_g_low = np.round(log_g_low, 2)
    
#E and zeta same across all files    
log_E_over_kT = data_low[:,0]  #emitted photon energy/kT
zeta_array = data_low[:,1] #cos psi

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
Inu_g_high_high = Inu_over_T_high[fixed_g_high] #high T, high g
Inu_g_high_low = Inu_over_T_low[fixed_g_high] #low T, high g
zeta_g_high = zeta_array[fixed_g_high]

fixed_g_low = np.where(log_g_array == log_g_low)
Inu_g_low_high = Inu_over_T_high[fixed_g_low] #high T, low g
Inu_g_low_low = Inu_over_T_low[fixed_g_low] #low T, low g
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
        Inu = Inu_over_T_high[fixed_g_high] 
        #Define upper and lower Inu  
        Inu_T_high = Inu_over_T_high
        Inu_T_low = Inu_over_T_low
        
print(len(Inu))

"""               
 #Let 1e-6 be effectively 0 
if zeta < 1e-06: 
    zeta = 1e-06
    
#Catch erroneous eta values 
if zeta > max(zeta_array) or zeta < min(zeta_array):
    #let 0 = 1e-6
    print("zeta is out of range.")
    print("Terminating program.")
    sys.exit()
elif zeta in zeta_array: 
    zeta_high = zeta
    zeta_low = zeta
else:
    #Find closest pair to given eta
    zeta_high, zeta_low = get_nearest_uneven(zeta_array, zeta)
#print(eta_high, eta_low)

#zeta_g_high and zeta_g_low are identical
fixed_zeta_high = np.where(zeta_g_high == zeta_high)[0]
fixed_zeta_low = np.where(zeta_g_high == zeta_low)[0]

Inu_np = np.asarray(Inu) #convert Inu to np array

#Interpolate zeta values
if zeta_high != zeta_low:
    Inu_zeta_high = Inu_np[fixed_zeta_high]
    Inu_zeta_low = Inu_np[fixed_zeta_low]
    #print("Interpolating eta.")
    
    Inu_zeta = [None]*len(Inu_zeta_high)
    for i in range(len(Inu_zeta_high)):
        Inu_zeta[i] = ((Inu_zeta_high[i] - Inu_zeta_low[i])/(zeta_high - zeta_low))\
        *(zeta - zeta_low) + Inu_zeta_low[i]
        
    Inu_zeta = np.asarray(Inu_zeta)
else:
    #print("No eta interpolation necessary.")
    fixed_zeta = np.where(zeta_g_high == zeta) #theta = pi/2
    Inu_zeta = Inu_np[fixed_zeta]
    
#Inu_high = Inu_T_high[fixed_eta_high] #highest T, highest eta, highest g
#Inu_low = Inu_T_low[fixed_eta_low] #lowest T, lowest eta, lowest g
E_zeta = E[fixed_zeta_high]

print(E_zeta, Inu_zeta)"""
