#!/usr/bin/env python
#-*- coding: utf-8 -*-

import numpy as np
import copy


##----------------------------------
#  Methode de Euler explicite :
def Euler(dt,Times,y0,f):
    sol = [y0];
    p0 = y0*1
    p1 = y0*0
    for t in Times:
        # Completer ICI :
        p1 = p0+dt*f(t,p0)
        sol.append(p1)
        p0 = p1*1.0
    sol.pop()
    return sol
##----------------------------------    

    
    
##----------------------------------
#  Methode de RK :
def RK(dt,Times,y0,f):
    sol = [y0];
    p0 = y0*1
    p1 = y0*0
    # Completer ICI :
    a = 1/2
    c1 = 0
    c2 = 1
    for t in Times:
        # Completer ICI :
        p1 = p0+dt*(c1*f(t,p0)+c2*f(t+dt*a,p0+dt*a*f(t,p0)))
        sol.append(p1)
        p0 = p1*1.0
    sol.pop()
    return sol
##----------------------------------
    
##----------------------------------
#  Methode de Multi-pas :
def MultiPas(dt,Times,y0,f):
    sol = [y0];
    
    # Initialisation :
    soltmp = RK(dt,np.array([0,dt]),y0,f)
    p0 = y0*1
    p1 = soltmp[1]
    p2 = y0*0
    
    f0 = f(Times[0],p0)
    for t in Times:
        # Completer ICI :
        p2 = p1 + dt*1/2*(3*f(t,p1)-f0)
        f0 = f(t,p1)
        sol.append(p1)
        p0 = p1*1.0
        p1 = p2*1.0
        
    sol.pop()
    return sol
##----------------------------------     
    
    
