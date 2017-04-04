# -*- coding: utf-8 -*-

from math import *
from numpy import *

# Autoc calcule les coefficients d'autocorellation des données extraites du fichier audio
# L est la liste des données
def autoc(L):
    n=len(L)
    # r stocke les coefficients auto regressifs
    r=[]
    r.append(1)
    for i in range(n-1):
        som=0
        for j in range(n-i):
            som+=L[j]*L[i+j]
        r.append(som)
    # M est la matrice d'autocorellation
    M=zeros((n-1,n-1))
    #initialisation de la diagonale
    for e in range(n-1):
        M[e,e]=1
    
    for d in range(len(r)-1):
        for k in range(n-d-2):
            M[d+k+1,k]=r[d+1]
            M[k,d+k+1]=r[d+1]
    M_inverse=linalg.inv(M)
    r2=array(r[1:])
    coeff=M.dot(-r2)
    
    return coeff
           
    
    
