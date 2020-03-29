#!/usr/bin/env python
#-*- coding: utf-8 -*-

import numpy as np
import affichage

import algo
import data

import time



"""
    Etude de modèles d'évolution de population  
"""

##----------------------------------
# But du programme :
# 1. Etude du modele de Verhulst.
# 2. Analyse convergence.
# 3. Etude du modele en epidemiologie
butPrgm = 1
print("Que souhaitez vous faire ?")
print(" 1. Simuler le modèle de Verhulst.")
print(" 2. Analyse de convergence.")
print(" 3. Simuler le modèle d'épidémiologie.")
butPrgm = input("Entrer le choix : ")

print("")
# Choix de la methode :
print("Choisir la méthode :")
print(" 1. Euler")
print(" 2. RK")
print(" 3. Multi-Pas")
choixMeth = input("Entrer votre choix : ")
##----------------------------------


##----------------------------------
# Partie I :
if (butPrgm=="1"):
    y0 = 0.2
    SOL = []
    
    ##----------------------------------
    # Paramètres 
    # Temps simulation
    T = 10
    # Discretisation
    dt = 0.1 
    Times = np.arange(0.0,T,dt)
    ##----------------------------------
    
    ##----------------------------------
    # Calcul solution approchée :
    start = time.time()
    if (choixMeth == "1"):
        SOL = algo.Euler(dt,Times,y0,data.FVerhulst)
    if (choixMeth == "2"):
        SOL = algo.RK(dt,Times,y0,data.FVerhulst)
    if (choixMeth == "3"):
        SOL = algo.MultiPas(dt,Times,y0,data.FVerhulst)
    print("Temps calculs : ",time.time() - start,"s")    
    ##----------------------------------
         
    ##----------------------------------
    # Affichage :     
    affichage.afficher([Times],[SOL])
    ##----------------------------------
##----------------------------------

##----------------------------------
# Partie II :
if (butPrgm=="2"):
    list_dt = [0.1,0.01,1e-3,1e-4,1e-5]
    list_err = []
    
    T = 3.0
    for dt in list_dt:
        Times = np.arange(0.0,T,dt)
        
        ll = 1.0
        y0 = 0.2
        SOL = []
        
        ##----------------------------------
        # Calcul solution approchée :
        start = time.time()
        if (choixMeth == "1"):
            SOL = algo.Euler(dt,Times,y0,data.FVerhulst)
        if (choixMeth == "2"):
            SOL = algo.RK(dt,Times,y0,data.FVerhulst)
        if (choixMeth == "3"):
            SOL = algo.MultiPas(dt,Times,y0,data.FVerhulst)
        print("Temps calculs : ",time.time() - start,"s")     
        ##----------------------------------
        
        ##----------------------------------
        # Calcul de la solution exacte :
        # Completer ICI :
        SOL_e = y0*1/(y0+(1-y0)*np.exp(-ll*Times))  
        ##----------------------------------
        
        ##----------------------------------
        # Analyse erreur :
        maxerr = max([np.sqrt(np.dot(e-SOL_e[i], e-SOL_e[i])) for i,e in enumerate(SOL)])
        print("dt = ",dt," et err =",maxerr)
        list_err.append(maxerr) 
        ##----------------------------------    
    
    ## Tracer courbe erreur :
    poly = np.polyfit(np.log10(list_dt)*(-1),np.log10(list_err),1)
    print("")
    print("Droite approchee : y=",poly[0],"x - ",np.abs(poly[1]))     
    affichage.afficherConv([np.log10(list_dt)*(-1)], [np.log10(list_err)])
##----------------------------------


##----------------------------------
# Partie III :
if (butPrgm=="3"):
    # Population initiale dans chaque categorie :
    # y = [Ns, Bs, As, Ss, Ni, Bi, Ai, Si, Nr, Br, Ar, Sr, Nx, Bx, Ax, Sx]
    y0 = np.array([0.5,1.0,10.0,1.0,0.0,0.0,0.001,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    #y0 = np.array([0.0,0.0,20.0,0.0,0.0,0.0,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    
    ##----------------------------------
    # Paramètres 
    # Temps simulation
    T = 100.0
    # Discretisation
    dt = 0.01
    Times = np.arange(0.0,T,dt)
    ##----------------------------------
    
    ##----------------------------------
    # Calcul solution approchée :
    SOL = []
    start = time.time()
    if (choixMeth == "1"):
        SOL = algo.Euler(dt,Times,y0,data.Fcomp)
    if (choixMeth == "2"):
        SOL = algo.RK(dt,Times,y0,data.Fcomp)
    if (choixMeth == "3"):
        SOL = algo.MultiPas(dt,Times,y0,data.Fcomp)
    print("Temps calculs : ",time.time() - start,"s") 
    ##----------------------------------
    
    ##----------------------------------
    # Population non immunisée
    N = np.array([e[0] for e in SOL])
    B = np.array([e[1] for e in SOL])
    A = np.array([e[2] for e in SOL])
    S = np.array([e[3] for e in SOL])
    # Population infectée 
    Ni = np.array([e[4] for e in SOL])
    Bi = np.array([e[5] for e in SOL])
    Ai = np.array([e[6] for e in SOL])
    Si = np.array([e[7] for e in SOL])
    # Population immunisée 
    Nr = np.array([e[8] for e in SOL])
    Br = np.array([e[9] for e in SOL])
    Ar = np.array([e[10] for e in SOL])
    Sr = np.array([e[11] for e in SOL])
    # Population decedée pour cause de la maladie (cumul)
    Nx = np.array([e[12] for e in SOL])
    Bx = np.array([e[13] for e in SOL])
    Ax = np.array([e[14] for e in SOL])
    Sx = np.array([e[15] for e in SOL])
    ##----------------------------------

    ##----------------------------------
    # Affichage des résultats 
    affichage.afficherMultiInfo([Times,Times,Times,Times], [N,B,A,S],[Ni,Bi,Ai,Si],[Nr,Br,Ar,Sr],[Nx,Bx,Ax,Sx])
    affichage.afficherPopTot([Times,Times,Times,Times], [N,B,A,S],[Ni,Bi,Ai,Si],[Nr,Br,Ar,Sr],[Nx,Bx,Ax,Sx])
    ##----------------------------------

##----------------------------------

print("FIN \n")
########################################
