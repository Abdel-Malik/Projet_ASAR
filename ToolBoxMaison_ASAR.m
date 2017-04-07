function res = ToolboxMP()
    res = autor([1,2,5,6,7],4);
end

function res = Echange_Ligne(M,lig1,lig2)
   t = M(lig1,:);
   M(lig1,:) = M(lig2,:);
   M(lig2,:) = t;
   res = M;
end

function res = Gauss_Inverse(M)
    I = eye(length(M));
    for e = 1:length(M)
        I(e,:) = I(e,:)/M(e,e);
        M(e,:) = M(e,:)/M(e,e);
        for i = (1:length(M))
            if(i~=e)
                I(i,:) = I(i,:)-I(e,:)*M(i,e);
                M(i,:) = M(i,:) - M(e,:)*M(i,e);
            end
        end
    end
    res = I;
end

function res = toeplitzAutoc(c,r)
    m=length(r);
    M = eye(m-1)*r(1);
    for d = 1:m
        for k = 1:m-d-1
            M(d+k,k) = r(d+1);
            M(k,d+k) = r(d+1);
        end
    end
	res = M;
end

function res = autoc(p,n)
	r=[];
    for i = 0:p
        % Calcul des coefficients d'auto-correlation r(i)
        som=0;
        for j = 1:n-i
            disp(i+j);
            som = som + L(j)*L(i+j);
        end
        r = [r,som];
    end
	res = r;
end

function res = autor(L,p)
    n=length(L);
    % r stocke les coefficients auto regressifs
    r=autoc(p,n);
    % M est la matrice d'auto-correlation liée aux données de L
    M=toeplitzAutoc(c,r);
    % M_inverse est l'inverse de M
    M_inverse=Gauss_Inverse(M);
    % r2 est r privé de sa premiere valeur
    r2= r(2:length(r));
    r2=r2';
    % coeff est un vecteur contenant les coefficients d'auto-régression
    coeff = M_inverse*(-r2);
    res = coeff;
end