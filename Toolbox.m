%extractions du son ainsi que sa fréquence d'échantillonage
[y2,Fs] = audioread('Z:\GitHub\Projet_ASAR\3Bonjours2.wav');
%récupération d'un seul canal
yi = y2(:,1);

%récupération d'un ensemble de point
y = yi(10001:11000);
%Calcul des coefficients d'autoRégressions
p=3;
n1=size(y,1);
autoc=(xcorr(y,y));
autoc=autoc(n1:n1+p);
n=size(autoc,1);
% autoc est le vecteur des coefficients d'autocorrelation
T=toeplitz(autoc(1:n-1));
T_inv=inv(T);
coeff=T_inv*(-autoc(2:n));
% Ou coeff=T\(-autoc(2:n))';
