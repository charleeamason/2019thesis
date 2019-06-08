#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 14:43:22 2019

@author: charlee

This program computes the radius of an oblate NS at a particular latitude
given the equatorial radius (R), spin frequency, and theta
"""
def get_R(theta, M, R, spin_freq):
    
    import numpy as np
    
    #Define constants:
    G = 6.67e-8
    c = 3e10
    
    #print("spin =", spin_freq)
    Omega = spin_freq*np.pi*2 #angular velocity at equator
    #print("Omega =", Omega)
    Omega_bar = Omega*(R**3/(G*M))**0.5 #dimensionless Omega
    #print("Omega bar =", Omega_bar)
    x = (G*M)/(R*c**2)
    
    sigma_2 = Omega_bar**2*(-0.788 + 1.030*x)
    R_theta = R*(1 + sigma_2*np.cos(theta)**2)
    
    return R_theta