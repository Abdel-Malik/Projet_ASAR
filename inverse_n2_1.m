%----------------Inverse une matrice de toeplitz quelconque
%Calcule en O(n^3)
% T = [1 2 3;2 1 2;3 2 1];
%T = [3 4 5 8 9 27;4 3 4 5 8 9;5 4 3 4 5 8;8 5 4 3 4 5;9 8 5 4 3 5;27 9 8 5 4 3]
%T = [6 8 7 4 1 38;8 6 8 7 4 1;7 8 6 8 7 4;4 7 8 6 8 7;1 4 7 8 6 8;38 1 4 7 8 6]
%T = [1 2;3 5]
function res = inverse_n2_1()
    T = [1 2 3;2 1 2;3 2 1];
    inv(T)
    inv_n2(T)
    
end

%Parametre T : une matrice carré
%Parametre n : indice de 1 à length(T)
%Spécification : renvoie un vecteur colonne de la premiere colonne de la matrice privé du premier élément T(1,1)
function y=r(T,n)
    S = zeros(n,1);
    for i=1:n
        S(i,1)=T(i+1,1);
    end
    y=S;
end
%Parametre v : un vecteur 
%Spécification : renvoie un vecteur ligne correspondant au vecteur tel que v(i) se place a l'indice n-i du nouveau vecteur.
function y=retourne(v)
    n=length(v);
    S = zeros(1,n);
    for i=0:(n-1)
        S(1,i+1) = v(n-i);
    end
    y=S;
end
%Parametre T : une matrice de Toeplit symétrique.
%Spécification : renvoie un réel : alpha, le vecteur colonne : u, le vecteur colonne : v  itérativement.
function[alpha,u]=Calcule(T)

    %Initialisation
    n=length(T(1,:));
    r0=T(1,1);
    Ti = 1/r0;
    Ui=cell(n);
    Ai =cell(n);
    i=1;
    ui = -Ti*r(T,i);
    alphai = r0 + ui'*r(T,i);
    Ui{i}=ui;
    Ai{i}=alphai;

    %itération
   for i=2:(length(T(1,:))-1)
       %Calcule de u(i)
       r1 = r(T,i);       
       ri = r1(i);
       ui = [Ui{i-1};0] - (ri + ui'*(retourne(r(T,i-1)))')/alphai * [(retourne(Ui{i-1}))';1];
       %Ajoute à la liste les nouveaux v(i+1) et u(i+1)
       Ui{i}=ui;
       %Calcule de alpha(i)
       alphai = r0 + Ui{i}'*r(T,i);        
       %Ajoute à la liste du nouveau alpha(i+1) 
       Ai{i}=alphai;
       
   end
    alpha=Ai{i};
    u=Ui{i};
end

function y=inve_n2(alpha,u,v)
    n=length(u)+1;
    S=zeros(n,n);
    S(1,1)=1/alpha;
    S(2:n,1)=u/alpha;
    S(1,2:n)=v'/alpha;
    U=1/alpha * (u*v' - retourne(v)'*retourne(u));
    for i=2:n
        for j=2:n
            S(i,j)=S(i-1,j-1)+U(i-1,j-1);
        end
    end
    
    y=S;
end

function y=inv_n2(T)
    [alpha,u]=Calcule(T);
    y=inve_n2(alpha,u,u);
end