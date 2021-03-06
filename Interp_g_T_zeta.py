"""
Created on Wed Aug  1 13:13:09 2018

@author: charlee

A program which, given a particular value for log T, log g, and zeta (cos psi), 
interpolates using discrete values from Ho et al. to produce an interpolated 
intensity profile. 
"""
import numpy as np
import glob
import sys
import os
import shutil

#Define functions to get nearest values for g, T, and eta

def find_closest(val1, val2, target):
    #returns closest given two test values and a target value
    return val2 if target - val1 >= val2 - target else val1

#to be used for T and g values
def get_nearest_pair(arr, interval, target):
    #returns two integers from list nearest to target
    #interval is spacing between integers in list
    n = len(arr)
    left = 0
    right = n - 1
    mid = 0
    
    while left < right:
        mid = (left + right)//2  # find the mid
        if target < arr[mid]:
            right = mid
        elif target > arr[mid]:
            left = mid + 1
        else:
            return (arr[mid], arr[mid])
    
    if target < arr[mid]:
        temp_1 = find_closest(arr[mid - 1], arr[mid], target)
        if target < temp_1: 
            temp_2 = temp_1 - interval
        else: 
            temp_2 = temp_1 + interval          
        return (temp_1, temp_2)
    else:
        temp_1 = find_closest(arr[mid + 1], arr[mid], target)
        if target < temp_1: 
            temp_2 = round(temp_1 - interval, 2)
        else: 
            temp_2 = round(temp_1 + interval, 2)         
        return (temp_1, temp_2) 
    
#E and eta values are not evenly spaced; requires different approach    
def get_nearest_uneven(array, value):
    array = np.asarray(array)
    array.sort()
    diff = array - value
    if value > max(array):
        print("Value is greater than any in array.")
        max_closest = max(array)
        min_closest = max(array)
    elif value < min(array):
        print("Value is smaller than any in array.")
        max_closest = min(array)
        min_closest = min(array)
    elif value in array:
        max_closest = value
        min_closest = value
    else:        
        neg_diff = []
        pos_diff = []      
        for i in range(len(diff)):
            if diff[i] < 0:
                neg_diff.append(array[i])
            elif diff[i] > 0: 
                pos_diff.append(array[i])  

        max_closest = pos_diff[0]       
        min_closest = neg_diff[-1]
            
    return max_closest, min_closest

def interp_T_and_g(log_T, log_g):

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
            Inu = Inu_over_T_high[fixed_g_high] 
            #Define upper and lower Inu  
            Inu_T_high = Inu_over_T_high
            Inu_T_low = Inu_over_T_low
        
    return E, Inu, Inu_T_high, Inu_T_low, zeta_array, zeta_g_high,\
            zeta_g_low
    
def interp_zeta(zeta, E, Inu, Inu_T_high, Inu_T_low, zeta_array, zeta_g_high,\
              zeta_g_low, log_g, log_T): 

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
    
    #eta_g_high and eta_g_low are identical
    fixed_zeta_high = np.where(zeta_g_high == zeta_high)
    fixed_zeta_low = np.where(zeta_g_high == zeta_low)
    
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
    
    """
    #Save spectra to a file
    head = 'E \t \t I_nu'
    final_data = np.array([E_eta, Inu_eta]).T
    np.savetxt('./Interp_Ho_Data/output_eta_' + str(format(eta, '.18f')) + '_g_' + str(log_g) + '_T_' + str(log_T) + '.txt', \
               final_data, fmt = '%.5e %.5e', header = head) """
    return E_zeta, Inu_zeta

def interp_E(E_zeta, Inu_zeta, E_list):
    #E_zeta and Inu_zeta are lists produced by zeta interp function
    #E_list is list of desired E/kT values
    Inu_final = np.zeros(len(E_list))
    E_np = np.around(E_zeta, 6)
    for i in range(len(E_list)):
         #Catch erroneous E values 
        if E_list[i] > max(E_zeta): 
            #print("value out of range--rounding down", i)
            E_high = max(E_zeta)
            E_low = max(E_zeta)
        elif E_list[i] < min(E_zeta):   
            #print("value out of range--rounding up", i)
            E_high = min(E_zeta)
            E_low = min(E_zeta)
        elif E_list[i] in E_zeta: 
            #print("No E interp necessary")
            E_high = E_list[i]
            E_low = E_list[i]         
        else:
            #print("Choosing nearest pair to E val in E_list",i)
            #Find closest pair to given E
            E_high, E_low = get_nearest_uneven(E_zeta, E_list[i])

        #Interpolate E values
        if E_high != E_low:
            #print("Interpolating E", i)
            E_high_round = np.around(E_high, 6) 
            E_low_round = np.around(E_low, 6)   
            E_high_idx = np.where(E_np == E_high_round)[0]
            E_low_idx = np.where(E_np == E_low_round)[0] 
             
            #returns single value for Inu corresponding to highest and lowest E
            Inu_high = Inu_zeta[E_high_idx]
            Inu_low = Inu_zeta[E_low_idx]
            #print("Interpolating E.")
        
            Inu_final[i] = ((Inu_high - Inu_low)/(E_high - E_low))\
            *(E_list[i] - E_low) + Inu_low
                
        else:
            #print("No E interpolation necessary.",i)
            Inu_final[i] = Inu_zeta[i]
    
    return Inu_final

def lookup_alpha(M_over_R_val, psi_val):
  
    import numpy as np
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
        #print("Exact psi is in table (within 4 decimal places).")
        psi_idx_temp = np.where(np.round(psi_val, 4) == psi_rounded)[0]
        psi_idx = int(psi_idx_temp)
        alpha_final = alpha[psi_idx][0]
        dzeta_dmu_final= dzeta_dmu[psi_idx][0]
        psi_final = psi_val
        print(psi_final, alpha_final, dzeta_dmu_final)
        
    elif psi_val > max_psi:
        #print("Psi value is larger than largest in table.")
        psi_idx = np.where(np.round(max_psi, 4) == psi_rounded)[0]
        alpha_final = alpha[psi_idx]
        dzeta_dmu_final= dzeta_dmu[psi_idx]
        psi_final = max_psi
        
    elif psi_val < min_psi: 
        #print("Psi value is smaller than smallest in table.")
        psi_idx= np.where(np.round(min_psi, 4) == psi_rounded)[0]
        alpha_final = alpha[psi_idx]
        dzeta_dmu_final= dzeta_dmu[psi_idx]
        psi_final = min_psi
        
    else: 
        #print("Interpolating.")
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
            
    return alpha_final, dzeta_dmu_final