from numpy import *
def Echange_Ligne(M,lig1,lig2):
    stockeLig=zeros(M[lig1:lig1+1].shape[1])
    for i in range(M[lig1:lig1+1].shape[1]):
        stockeLig[i]=M[lig1,i]
    M[lig1:lig1+1]=M[lig2:lig2+1]
    M[lig2:lig2+1]=stockeLig
    return M

##def Copie_Matrice(M):
##    n=M.shape[0]
##    N=zeros((n,n))
##    for i in range(n):
##        for j in range(n):
##            N[i,j]=M[i,j]
##    return N

# Gauss_Inverse retourne l'inverse de la matrice M en utilisant
# la méthode de Gauss-Jordan sur M augmentee de la matrice identité
def Gauss_Inverse(M):
    n=M.shape[0]
    # I est la matrice identité
    I=eye(n,n)
    #Lpiv est l'indice du dernier pivot trouvé
    Lpiv=-1
    for j in range(n):
        # On recherche le pivot maximal, k est sa ligne
        maxi=M[Lpiv,j]
        for i in range(Lpiv+1,n):
            if M[i,j]>=maxi:
                k=i
        if abs(M[k,j])>10**(-10): # pour éviter les erreurs d'approximation
        #if M[k,j]!=0:
            Lpiv+=1
            coeff=M[k,j]
            M[k:k+1]=M[k:k+1]/coeff
            I[k:k+1]=I[k:k+1]/coeff
            Echange_Ligne(M,k,Lpiv)
            Echange_Ligne(I,k,Lpiv)
            for i in range(n):
                if i!=Lpiv:
                    coeff2=M[i,j]
                    M[i:i+1]=M[i:i+1]-coeff2*M[Lpiv:Lpiv+1]
                    I[i:i+1]=I[i:i+1]-coeff2*I[Lpiv:Lpiv+1]
        
    return I

	