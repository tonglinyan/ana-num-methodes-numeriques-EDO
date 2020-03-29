#!/usr/bin/env python
#-*- coding: utf-8 -*-

import numpy as np
import copy


##----------------------------------
# Cas SIRX :
def FVerhulst(t,y):
    # Lambda :
    ll = 1.0
    #Completer ICI :
    return (ll*y*(1-y))
##----------------------------------
    

##----------------------------------
# Cas comp vaccin :
def Fcomp(t,y):
    # y = [N B A S Ni Bi Ai Si Nr Br Ar Sr Nx Bx Ax Sx]
    res = y*0.0
    [Ns,Bs,As,Ss,Ni,Bi,Ai,Si,Nr,Br,Ar,Sr,Nx,Bx,Ax,Sx] = y
    
    # Seuil a partir du quel on considere qu'il n'y a plus d'individus dans la categorie :
    seuil = 1e-7
    for ite,e in enumerate(y):
        if(e < seuil):
            y[ite] = 0.0
    
    ##----------------------------------
    ## Parametres :
    # Vieillisement naturel :
    c = 0.04
    cnb = 0.5
    cba = 0.2
    cas = 0.02
    m = 0.2
    
    Leslie = np.array([[-cnb,0,0,0],[cnb,-cba,0,0],[0,cba,-cas,0],[0,0,cas,-m]])
    
    # Infection :
    p = 0.2
    mni = 0.99
    mbi = 0.1
    mai = 0.3
    msi = 0.98
    rni = 0.01
    rbi = 0.9
    rai = 0.7
    rsi = 0.02
    
    # Politique de vaccination
    Tf = 2.0
    vns =1.0 if t <= Tf else 0.8
    vbs =0.8 if t <= Tf else 0.2
    vas =0.0 if t <= Tf else 0.2
    vss =1.0 if t <= Tf else 0.2
    ##----------------------------------
    
    ##----------------------------------
    # Vieillisement naturel :
    res[0:4] = np.dot(Leslie,[Ns,Bs,As,Ss])
    ## Completer ICI :
    res[4:8] = np.dot(Leslie,[Ni,Bi,Ai,Si])
    res[8:12] = np.dot(Leslie,[Nr,Br,Ar,Sr])
    # Croissance :
    res[0] = res[0] + c*(As+Ai+Ar)
    ##----------------------------------
    
    ##----------------------------------
    # Effet de l'infection :
    #Population totale infectee :
    PopI = (Ni + Bi + Ai + Si)
    # Sain qui tombe malade
    res[0:4] = res[0:4] - p * np.array([Ns,Bs,As,Ss])*PopI
    ## Completer ICI :
    res[4:8] = res[4:8] + p * np.array([Ns,Bs,As,Ss])*PopI
    # Malade qui deviennent soit remis, soit decedes
    res[4:8] = res[4:8] - np.array([Ni,Bi,Ai,Si]) * np.array([rni+mni,rbi+mbi,rai+mai,rsi+msi])
    ## Completer ICI :
    res[8:12] = res[8:12] + np.array([Ni,Bi,Ai,Si]) * np.array([rni,rbi,rai,rsi])
    res[12:16] = res[12:16] + np.array([Ni,Bi,Ai,Si]) * np.array([mni,mbi,mai,msi])
    ##----------------------------------
    
    ##----------------------------------
    # Politique de vaccination :
    res[0:4] = res[0:4] - np.array([Ns,Bs,As,Ss]) * np.array([vns,vbs,vas,vss])
    ## Completer ICI :
    res[8:12] = res[8:12] + np.array([Ns,Bs,As,Ss]) * np.array([vns,vbs,vas,vss])
    ##----------------------------------
    
    return res;
##----------------------------------




##----------------------------------
# Cas comp vaccin :
def FcompSav(t,y):
    ## Modele :
    # N' = c * (A + Ai + Ar) - p * N * (Ni*0.25 + Bi*0.25 + Ai*0.5 + Si*0.25) - nb * N - vn * N
    # B' = nb * N - p * B * (Ni*0.25 + Bi*0.25 + Ai*0.5 + Si*0.25) - ba * B  - vb * B
    # A' = ba * B - p * A * (Ni*0.25 + Bi*0.25 + Ai*0.5 + Si*0.25) - aas * A - va * A 
    # S' = aas * A - p * S * (Ni*0.25 + Bi*0.25 + Ai*0.5 + Si*0.25) - vs * S - m * S
    # Ni' = p * N * (Ni*0.25 + Bi*0.25 + Ai*0.5 + Si*0.25) - (rn+mn) * Ni 
    # Bi' = p * B * (Ni*0.25 + Bi*0.25 + Ai*0.5 + Si*0.25) - (rb+mb) * Bi
    # Ai' = p * A * (Ni*0.25 + Bi*0.25 + Ai*0.5 + Si*0.25) - (ra+ma) * Ai
    # Si' = p * S * (Ni*0.25 + Bi*0.25 + Ai*0.5 + Si*0.25) - (rs+ms+m) * Si
    # Nr' = vn * N + rn * Ni - nb * Nr
    # Br' = vb * B + rb * Bi - ba * Br
    # Ar' = va * A + ra * Ai - aas * Ar
    # Sr' = vs * S + rs * Si - m * Sr 
    # Nx' = mn * Ni
    # Bx' = mb * Bi
    # Ax' = ma * Ai
    # Sx' = (ms + m) * Si + m * (S + Sr)
    
    # Parametres :
    c = 0.04
    p = 0.2
    nb = 0.5
    ba = 0.2
    aas = 0.02
    m = 0.2
    mn = 0.99
    mb = 0.1
    ma = 0.3
    ms = 0.98
    rn = 0.01
    rb = 0.9
    ra = 0.7
    rs = 0.02
    
    # Politique de vaccination
    Tf = 10.0
    vn = 1.0 * (t <= Tf)
    vb = 0.8 * (t <= Tf)
    va = 0.2 * (t <= Tf)
    vs = 0.0 * (t <= Tf)
    
    # Comportement periodique 
    # c = 1.0
#     K = 0.0
#     p = 0.1
#     v = 0.3
#     m = 0.8
#     g = 0.002
    
    # y = [N B A S Ni Bi Ai Si Nr Br Ar Sr Nx Bx Ax Sx]
    res = y*1.0
    [Ns,Bs,As,Ss,Ni,Bi,Ai,Si,Nr,Br,Ar,Sr,Nx,Bx,Ax,Sx] = y
    
    # Seuil a partir du quel on considere qu'il n'y a plus d'individus infectes dans la categorie :
    seuil = 1e-7
    if(Ni < seuil):
        Ni = 0.0
    if(Bi < seuil):
        Bi = 0.0
    if(Ai < seuil):
        Ai = 0.0
    if(Si < seuil):
        Si = 0.0
    
    #PopI = (Ni*0.125 + Bi*0.125 + Ai*0.5 + Si*0.25)
    PopI = (Ni + Bi + Ai + Si)
    
    res[0] = c * (As + Ar + Ai) - p * Ns * PopI - nb * Ns - vn * Ns
    res[1] = nb * Ns - p * Bs * PopI - ba * Bs  - vb * Bs
    res[2] = ba * Bs - p * As * PopI - aas * As - va * As
    res[3] = aas * As - p * Ss * PopI - vs * Ss - m * Ss
    res[4] = p * Ns * PopI - (rn+mn) * Ni - nb * Ni
    res[5] = p * Bs * PopI - (rb+mb) * Bi + nb * Ni - ba * Bi
    res[6] = p * As * PopI - (ra+ma) * Ai + ba * Bi - aas * Ai
    res[7] = p * Ss * PopI - (rs+ms+m) * Si + aas * Ai
    res[8] = vn * Ns - nb * Nr + rn * Ni
    res[9] = vb * Bs + nb * Nr + rb * Bi - ba * Br
    res[10] = va * As + ba * Br + ra * Ai - aas * Ar
    res[11] = vs * Ss + aas * Ar + rs * Si - m * Sr
    res[12] = mn * Ni
    res[13] = mb * Bi
    res[14] = ma * Ai
    res[15] = ms * Si #(ms + m) * Si + m * (S + Sr)
    
    return res;
##----------------------------------



##----------------------------------
# Cas SIRX :
def Fsirx(t,y):
    ## Modele :
    # S' = cS S -K S^2 - p SI - v S + cR R
    # I' = p SI - (m + g) I 
    # R' = g I + v S  
    # X' = m I
    
    # Parametres :
    cS = 0.2
    cR = 0.2
    K = 0.0
    p = 0.2
    v = 0.0
    m = 0.8
    g = 0.2
    
    # Comportement periodique 
    # c = 1.0
#     K = 0.0
#     p = 0.1
#     v = 0.3
#     m = 0.8
#     g = 0.002
    
    # y = [S I R X]
    res = y*1.0
    
    res[0] = cS*100 - K * y[0] * y[0] - p * y[0] * y[1] - v * y[0] #+ cS * y[0]+ cR * y[2]
    res[1] = p * y[0] * y[1] - (m+g) * y[1]
    res[2] = g * y[1] + v * y[0]
    res[3] = m * y[1]
    
    return res;
##----------------------------------

