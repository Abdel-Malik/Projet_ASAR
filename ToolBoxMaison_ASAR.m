function res = ToolBoxMaison_ASAR()
    %extractions du son ainsi que sa fréquence d'échantillonage
    [y2,Fs] = audioread('3Bonjours2.wav');
    %récupération d'un seul canal
    yi = y2(:,1);
    %récupération d'un ensemble de points
    y = yi(10001:11000);
    %Calcul des coefficients d'autoRégressions
    p=13;
    res = autor(y,p);
end

%Méthode de Gauss
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


%Création d'une matrice de Toeplitz symétrique en se servant du vecteur r
function res = toeplitzAutoc(r)
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



function res = autoc(p,L)
    n=length(L);
    r=[];
    for i = 0:p
        % Calcul des coefficients d'auto-correlation r(i)
        som=0;
        for j = 1:n-i
            som = som + L(j)*L(i+j);
        end
        r = [r,som];
    end
    res = r;
end


function res = autor(L,p)
    % r stocke les coefficients auto regressifs
    r=autoc(p,L);
    % M est la matrice d'auto-correlation liée aux données de L
    M=toeplitzAutoc(r);
    % M_inverse est l'inverse de M
    M_inverse=Gauss_Inverse(M);
    % r2 est r privé de sa premiere valeur
    r2= r(2:length(r));
    r2=r2';
    % coeff est un vecteur contenant les coefficients d'auto-régression
    coeff = M_inverse*(-r2);
    res = coeff;
end
