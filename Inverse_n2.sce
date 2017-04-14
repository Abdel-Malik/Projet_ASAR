//----------------Inverse une matrice de toeplitz quelconque
//Calcule en O(n^3)

//Parametre T : une matrice carré
//Parametre n : indice de 1 à length(T)
//Spécification : renvoie un vecteur ligne de la premiere ligne de la matrice privé du premier élément T(1,1)
function y=s(T,n)
    S = zeros(1,n)
    for i=1:n
        S(1,i)=T(1,i+1)
    end
    y=S
endfunction
//Parametre T : une matrice carré
//Parametre n : indice de 1 à length(T)
//Spécification : renvoie un vecteur colonne de la premiere colonne de la matrice privé du premier élément T(1,1)
function y=r(T,n)
    S = zeros(n,1)
    //disp(S)
    for i=1:n
        S(i,1)=T(i+1,1)
        //disp(S)
    end
    y=S
endfunction
//Parametre v : un vecteur 
//Spécification : renvoie un vecteur ligne correspondant au vecteur tel que v(i) se place a l'indice n-i du nouveau vecteur.
function y=retourne(v)
    n=length(v)
    S = zeros(1,n)
    for i=0:(n-1)
        S(1,i+1) = v(n-i)
    end
    y=S
endfunction
//Parametre T : une matrice de Toeplit symétrique.
//Spécification : renvoie un réel : alpha, le vecteur colonne : u, le vecteur colonne : v  itérativement.
function[alpha,u,v]=Calcule(T)

//Initialisation
r0=T(1,1)
Ti = 1/r0
Vi=list()
Ui=list()
Ai =list()
i=1
vi = -(Ti)'*(s(T,i))'
ui = -Ti*r(T,i)
alphai = r0 + vi'*r(T,i)
Vi($+1)=vi
Ui($+1)=ui
Ai($+1)=alphai

//itération
for i=2:(length(T(1,:))-1) 
//Calcule de v(i)
r1 = s(T,i)
ri = r1(i)
vi = [Vi(1);0] - (ri + vi'*(retourne(s(T,i-1)))')/alphai * [(retourne(Ui(1)))';1]
//Calcule de u(i)
r1 = r(T,i)         
ri = r1(i)
ui = [Ui(1);0] - (ri + ui'*(retourne(r(T,i-1)))')/alphai * [(retourne(Vi(1)))';1]
//Ajoute à la liste les nouveaux v(i+1) et u(i+1)
Vi(0)=vi
Ui(0)=ui
//Calcule de alpha(i)
alphai = r0 + Vi(1)'*r(T,i)        
//Ajoute à la liste du nouveau alpha(i+1) 
Ai(0)=alphai
end


alpha=Ai(1)
u=Ui(1)
v=Vi(1)

disp(Vi(1))
disp(Ui(1))
disp(Ai(1))
endfunction
//Parametre alpha : réel
//Parametre u : vecteur colonne
//Parametre v : vecteur colonne
//Spécification : Renvoie la matrice auto-généré par alpha , u et v.
function y=generation_Matrice(alpha,u,v)
    n=length(u)+1
    S=zeros(n,n)
    S(1,1)=1/alpha
    S(2:n,1)=u/alpha
    S(1,2:n)=v'/alpha
    U=1/alpha * (u*v' - retourne(v)'*retourne(u))
    for i=2:n
        for j=2:n
            S(i,j)=S(i-1,j-1)+U(i-1,j-1)
        end
    end
    y=S
endfunction

T = [6 8 7 4 1 38;8 6 8 7 4 1;7 8 6 8 7 4;4 7 8 6 8 7;1 4 7 8 6 8;38 1 4 7 8 6]
//T = [3 8 7 2 37;8 3 8 7 2;7 8 3 8 7;2 7 8 3 8;37 2 7 8 3]
//T = [3 2 1 5;2 3 2 1;1 2 3 2;5 1 2 3]
//T = [4 8 23 7;8 4 8 23;23 8 4 8;7 23 8 4]
//T = [3 2 1;2 3 2;1 2 3]
//T = [4 2 1 37;4 4 2 1;5 4 4 2;26 5 4 4]
[a,u,v]=Calcule(T)
M=generation_Matrice(a,u,v)
disp(M)
