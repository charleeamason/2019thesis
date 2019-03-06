#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 15:20:30 2019

@author: charlee
"""
import numpy as np 

value = 0.123213
array = np.random.uniform(low=0, high=1.8, size=(1000,))
array = np.asarray(array)
array.sort()
diff = array - value

import time
start_time = time.time()

#print(diff)
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
    print("--- %s seconds to get pair ---" % (time.time() - start_time))
    
    """   
    
    i = 0
    while diff[0]*diff[i] > 0 and i < len(diff):
        i = i + 1
    print(i, diff[i]) 
    
max_closest = array[i]
min_closest = array[i - 1] 
print("--- %s seconds to get pair ---" % (time.time() - start_time)) """

print(max_closest, min_closest)
            