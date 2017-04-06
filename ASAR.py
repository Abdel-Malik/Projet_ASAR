from numpy import *
# Autoc calcule les coefficients d'autocorellation des données extraites du fichier audio
# L est la liste des données
# p est l'ordre du modèle auto regressif
def autoc(L,p):
    n=len(L)
    # r stocke les coefficients auto regressifs
    r=[]
    for i in range(p+1):
        # Calcul des coefficients d'auto-correlation r(i)
        som=0
        for j in range(n-i):
            som+=L[j]*L[i+j]
        r.append(som)
    # M est la matrice d'auto-correlation liée aux données de L
    
    m=len(r)
    M=zeros((m-1,m-1))
    for d in range(m):
        for k in range(m-d-1):
            M[d+k,k]=r[d]
            M[k,d+k]=r[d]
    # M_inverse est l'inverse de M
    M_inverse=linalg.inv(M)
    # r2 est r privé de sa premiere valeur
    r2=array(r[1:])
    # coeff est un vecteur contenant les coefficients d'auto-régression
    coeff=M_inverse.dot(-r2)
    
    return coeff
