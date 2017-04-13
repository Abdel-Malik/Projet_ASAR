#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 18:00:01 2017

@author: ledentp
"""

# Caclul du vecteur contenant les coefficients d'autocorrelation
# Yn est le vecteur des donnes
# R est un vecteur contenant les coefficitens d'autocorrelation
# Complexite O(n)
def coeffsAuto(Yn):
    for i in range(0, len(Yn)-1):
        R[i] = Yn[i] * Yn[i+1];
    return R;

# Application de l'algorithme de Levinson-Durbin pour le modele AR(p) a partir des coefficients d'autocorrelation
# Entrees : p est l'ordre du modele
#           R est le vecteur des coefficients d'autocorellation
# Sortie  : A est le vecteur contenant les coefficients du modele AR
# Complexite : O(n^2)
def RecLevinsonDurbin(p,R):
    # Variables locales :
    # E est un vecteur d'inconnues
    # lambd est un nombre dans une combinaison lineaire
    # Initialisation
    Ak = zeros(p+1)
    Ak[0] = 1;
    Ek   = R[0]

    # Recurtion de Levinson-Durbin 
    for k in range(p):
        # calcul de lambda
        lambd = 0
        for j in range(k+1):
            lambd -= Ak[j]*R[k+1-j]
        lambd /= Ek
        
        # Mise a jour de Ak
        for n in range(1+int((k+1)/2)):
            temp = Ak[k+1-n]+lambd*Ak[n]
            Ak[n]=Ak[n]+lambd*Ak[k+1-n]
            Ak[k+1-n] = temp
        
        # Mise a jour de Ek
        Ek *= 1-lambd*lambd
    
    return Ak

    
# Application de l'algorithme de Levinson-Durbin pour le modele AR(p) a partir des donnes
def LevinsonDurbin(Yn,p):
    return RecLevinsonDurbin(p,coeffsAuto(Yn));
