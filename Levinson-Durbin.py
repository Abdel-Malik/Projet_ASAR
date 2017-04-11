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
    # somme est un nombre
    # Initialisation
    A[0] = 1;
    A[1] = -R[1] / R[0];
    E[1] = R[0] + A[1] * R[1];
    
    # Recurrence
    for k in range(1, p):
        # Calcul de la somme de j = 0 a k des A[j] * R[k+1-j]
        somme = 0;
        for j in range(0,k):
            somme = somme + A[j] * R[k+1-j];
        # calcul de lambda
        lambd = -somme / E(k);
        # calcul de A[k+1] = U[k+1] + lambda * V[k+1]
        # U = [A,0]
        # V contient les coefficients de U dans l'ordre inverse
        for x in range(0,k):
            U[x]   = A[x];
            V[x+1] = A[k-x];
        U[k+1] = 0;
        V[0]   = A[1];
        A[k+1] = U[k+1] + lambd * V[k+1]
        # Mise a jours de E
        E[k+1] = (1 - lambd**2) * E[k] # = E[k] + lambda * somme
    return A;

    
# Application de l'algorithme de Levinson-Durbin pour le modele AR(p) a partir des donnes
def LevinsonDurbin(Yn,p):
    return RecLevinsonDurbin(p,coeffsAuto(Yn));