function res = Echange_Ligne(M,lig1,lig2)
    for i = (1:length(M(lig1)))
        stockeLig[i]=M[lig1,i]
	end
    M(lig1:lig1+1)=M(lig2:lig2+1)
    M(lig2:lig2+1)=stockeLig
    res = M
end

function res = Gauss_Inverse(M)
    n=length(M,1)
    % I est la matrice identité
    I=kDiag(n,n)
    %Lpiv est l'indice du dernier pivot trouvé
    Lpiv=-1
    for j =(1:n)
        %On recherche le pivot maximal, k est sa ligne
        maxi=M(Lpiv,j)
        for i =(Lpiv+1:n):
            if(M(i,j)>=maxi)
                k=i
			end
        if abs(M(k,j))>10**(-10) % pour éviter les erreurs d'approximation
            Lpiv+=1
            coeff=M(k,j)
            M(k:k+1)=M(k:k+1)/coeff
            I(k:k+1)=I(k:k+1)/coeff
            Echange_Ligne(M,k,Lpiv)
            Echange_Ligne(I,k,Lpiv)
            for i = (1:n)
                if i!=Lpiv
                    coeff2=M(i,j)
                    M(i:i+1)=M(i:i+1)-coeff2*M(Lpiv:Lpiv+1)
                    I(i:i+1)=I(i:i+1)-coeff2*I(Lpiv:Lpiv+1)
				end
			end
		end
	end
    res = I
end

% def autoc(L,p):
    % n=len(L)
    % # r stocke les coefficients auto regressifs
    % r=[]
    % for i in range(p+1):
        % # Calcul des coefficients d'auto-correlation r(i)
        % som=0
        % for j in range(n-i):
            % som+=L[j]*L[i+j]
        % r.append(som)
    % # M est la matrice d'auto-correlation liée aux données de L
    
    % m=len(r)
    % M=zeros((m-1,m-1))
    % for d in range(m):
        % for k in range(m-d-1):
            % M[d+k,k]=r[d]
            % M[k,d+k]=r[d]
    % # M_inverse est l'inverse de M
    % M_inverse=Gauss_Inverse(M)
    % # r2 est r privé de sa premiere valeur
    % r2=array(r[1:])
    % # coeff est un vecteur contenant les coefficients d'auto-régression
    % coeff=M_inverse.dot(-r2)
    
    % return coeff


    