#!/usr/bin/env python
#-*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import data

from pylab import *        


def afficher(list_x,list_y):
    for j,absic in enumerate(list_x):
        plt.plot(absic,list_y[j], '-',linewidth=4.0)
    
    plt.show()
    
def afficherConv(list_x,list_y):
    for j,absic in enumerate(list_x):
        plt.plot(absic,list_y[j], 'x-',linewidth=4.0)
    
    plt.legend(["log10(err) vs log10(1/dt)"])
    plt.show()    
    
    
    
def afficherMultiInfo(list_x,list_y1,list_y2,list_y3,list_y4):
    plt.subplot(2,2,1)
    for j,absic in enumerate(list_x):
        plt.plot(absic,list_y1[j], '-',linewidth=4.0)
    plt.legend(["N","B","A","S"])
    plt.title('Population non immunisee')
    
    plt.subplot(2,2,2)
    for j,absic in enumerate(list_x):
        plt.plot(absic,list_y2[j], '-',linewidth=4.0)
    plt.legend(["Ni","Bi","Ai","Si"])
    plt.title('Population infectee')
    
    plt.subplot(2,2,3)
    for j,absic in enumerate(list_x):
        plt.plot(absic,list_y3[j], '-',linewidth=4.0)
    plt.legend(["Nr","Br","Ar","Sr"])
    plt.title('Population immunisee')
    
    plt.subplot(2,2,4)
    dt = absic[1] - absic[0]
    for j,absic in enumerate(list_x):
        ll = len(list_y4[j])
        diff = (list_y4[j][1:ll] - list_y4[j][0:ll-1]) / dt 
        plt.plot(absic[0:ll-1],diff, '-',linewidth=4.0)
    plt.legend(["Nx","Bx","Ax","Sx"])
    plt.title('Population decedee (cause maladie)')
    
    plt.show()
    
    
def afficherPopTot(list_x,list_y1,list_y2,list_y3,list_y4):
    SainNonImu = list_y1[0]+list_y1[1]+list_y1[2]+list_y1[3]
    Inf = list_y2[0]+list_y2[1]+list_y2[2]+list_y2[3]
    Imu = list_y3[0]+list_y3[1]+list_y3[2]+list_y3[3]
    
    absic = list_x[0]
    plt.plot(absic,SainNonImu+Inf+Imu,'k-',linewidth=4.0)
    plt.plot(absic,SainNonImu,'b--',linewidth=2.0)
    plt.plot(absic,Inf,'r--',linewidth=2.0)
    plt.plot(absic,Imu,'g--',linewidth=2.0)
    
    plt.legend(["Total","Saine NI","Infec.","Imu."])
    plt.title('Population Totale')
    
    plt.show()   
