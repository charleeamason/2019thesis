#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 16:07:49 2019

@author: charlee

This function determines the gravitational acceleration g of a NS given mass (M), 
equatorial radius (R), spin frequency, and theta
"""
def get_g(theta, R, M, M_over_R, spin_freq):
    import numpy as np
    
    #Define constants
    c = 3e10 #speed of light cm/s 
    G = 6.67e-8
    
    #Initial parameters:
    z = (1 - 2*M_over_R)**-0.5 - 1
    #print("spin =", spin_freq)
    
    Omega = spin_freq*np.pi*2 #angular velocity at equator
    #print("Omega =", Omega)
    Omega_bar = Omega*(R**3/(G*M))**0.5 #dimensionless Omega
    #print("Omega bar =", Omega_bar)
    
    g_eq = (M*G/R**2)*(z + 1)
    #print("log geq =", np.log10(g_eq))
    x = (G*M)/(R*c**2)
    
    #fixed quantities from Mohammed AlGendy & Sharon 2014
    c_e = (-0.791  + 0.766*x)*Omega_bar**2
    c_p = (1.138 - 1.431*x)*Omega_bar**2 
    d_e = (-1.315 + 2.431*x**2)*Omega_bar**4
    f_e = -1.172*x*Omega_bar**6
    d_p = (0.653 - 2.864*x**2)*Omega_bar**4
    f_p = 0.975*x*Omega_bar**6
    d_60 = (13.47*x - 27.13*x**2)*Omega_bar**4
    #f_60 = 1.69*Omega_bar**6
    
    if Omega_bar**2 <= 0.1:
        #slow rotation limit
        g = (1 + c_e*np.sin(theta)**2 + c_p*np.cos(theta)**2)*g_eq
    
        #print("g =", g)
    
    else: 
        #rapid rotation 
        g = (1 + (c_e + d_e + f_e)*np.sin(theta)**2 +\
                   (c_p + d_p + f_p + d_60)*np.cos(theta)**2\
                   + d_60*np.cos(theta))*g_eq
    #print("g =", g)
    
    return g