<<<<<<< HEAD
# -*- coding: utf-8 -*-

from math import *
=======

>>>>>>> origin/master
from numpy import *
# Autoc calcule les coefficients d'autocorellation des données extraites du fichier audio
# L est la liste des données
# p est l'ordre du modèle auto regressif (ne prend en compte que
# les n premieres valeurs de l'echantillon
def autoc(L,p):
    L=L[:p+1]
    n=len(L)
    # r stocke les coefficients auto regressifs
    r=[]
    for i in range(n):
        # Calcul des coefficients d'auto-correlation r(i)
        som=0
        for j in range(n-i):
            som+=L[j]*L[i+j]
        r.append(som)
    # M est la matrice d'auto-correlation liée aux données de L
    M=zeros((n-1,n-1))
    for d in range(len(r)):
        for k in range(n-d-1):
            M[d+k,k]=r[d]
            M[k,d+k]=r[d]
    # M_inverse est l'inverse de M
    M_inverse=linalg.inv(M)
    # r2 est r privé de sa premiere valeur
    r2=array(r[1:])
    # coeff est un vecteur contenant les coefficients d'auto-régression
    coeff=M_inverse.dot(-r2)
    
    return coeff
           
    
    
