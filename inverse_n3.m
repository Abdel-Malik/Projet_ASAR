%----------------Inverse une matrice de toeplitz quelconque
%Calcule en O(n^3)
% T = [1 2 3;2 1 2;3 2 1];
%T = [3 4 5 8 9 27;4 3 4 5 8 9;5 4 3 4 5 8;8 5 4 3 4 5;9 8 5 4 3 5;27 9 8 5 4 3]
%T = [3 2;2 3]
%T = [3 2 1;2 3 2;1 2 3]
%T = [3 2 1 5;2 3 2 1;1 2 3 2;5 1 2 3]
%T = [37 2 5;3 37 2;1 3 37]
%T = [1 3 4 5;2 1 3 4;6 2 1 3;7 6 2 1]
%T = [3 2 1 5;2 3 2 1;1 2 3 2;5 1 2 3]
%T = [3 2 4 7 8 23;4 3 2 4 7 8;21 4 3 2 4 7;7 21 4 3 2 4;3 7 21 4 3 2;1 3 7 21 4 3]
%T = [4 8 23 7;8 4 8 23;23 8 4 8;7 23 8 4]
%T = [3 8 7 2 37;8 3 8 7 2;7 8 3 8 7;2 7 8 3 8;37 2 7 8 3]
%T = [6 8 7 4 1 38;8 6 8 7 4 1;7 8 6 8 7 4;4 7 8 6 8 7;1 4 7 8 6 8;38 1 4 7 8 6]
%!!!!!!
%T = [4 2 1 37;4 4 2 1;5 4 4 2;26 5 4 4]
%T = [3 2 9 17 2;40 3 2 9 17;22 40 3 2 9;1 22 40 3 2;8 1 22 40 3]
%T = [6 8 7 4 1 38;8 6 8 7 4 1;7 8 6 8 7 4;4 7 8 6 8 7;1 4 7 8 6 8;38 1 4 7 8 6]
%T = [3 4 5 8 9 27;4 3 4 5 8 9;5 4 3 4 5 8;8 5 4 3 4 5;9 8 5 4 3  4;27 9 8 5 4 3]
function res = inverse_n3()
   T = [3 8 7 2 37;8 3 8 7 2;7 8 3 8 7;2 7 8 3 8;37 2 7 8 3];
    inv(T)
    inv_n3(T);
    
end
%Parametre T : une matrice carr�
%Parametre n : indice de 1 � length(T)
%Sp�cification : renvoie un vecteur ligne de la premiere ligne de la matrice priv� du premier �l�ment T(1,1)
function y = s(T,n)
    S = zeros(1,n);
    for i=1:n
        S(1,i)=T(1,i+1);
    end
    y=S;
end
%Parametre T : une matrice carr�
%Parametre n : indice de 1 � length(T)
%Sp�cification : renvoie un vecteur colonne de la premiere colonne de la matrice priv� du premier �l�ment T(1,1)
function y=r(T,n)
    S = zeros(n,1);
    for i=1:n
        S(i,1)=T(i+1,1);
    end
    y=S;
end
%Parametre v : un vecteur 
%Sp�cification : renvoie un vecteur ligne correspondant au vecteur tel que v(i) se place a l'indice n-i du nouveau vecteur.
function y=retourne(v)
    n=length(v);
    S = zeros(1,n);
    for i=0:(n-1)
        S(1,i+1) = v(n-i);
    end
    y=S;
end
%Parametre u : un vecteur colonne de taille n-1
%Parametre v : un vecteur colonne de taille n-1
%Parametre T : une  sous-matrice carr� de taille n-1
%Sp�cification : Renvoie une matrice carr� de taille n.
function y=construireMatrice(u,v,T,n)
    M = zeros(n,n);
    M(1,1)= 1;
    M(2:n,1)=u;
    M(1,2:n)=v';
    M(2:n,2:n)=T;
    y=M;
    
end




function y=inv_n3(T)
r0 = T(1,1);
n=length(T(1,:))-1;
D = zeros(n+1);
    Ti = 1/r0;
    
    for i=1:n
        vi = -(Ti)'*(s(T,i))';
        ui = -Ti*r(T,i);
        alpha = r0 + vi'*r(T,i);
       sous_matrice = alpha*Ti+ ui*vi';
       Ti = construireMatrice(ui,vi,sous_matrice,i+1);
        Ti = 1/alpha * Ti;
    end
    y=Ti
end