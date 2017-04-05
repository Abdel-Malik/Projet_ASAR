%extractions du son ainsi que sa fr�quence d'�chantillonage
[y2,Fs] = audioread('pianoLa.wav');
%r�cup�ration d'un seul canal
yi = y2(:,1);

%r�cup�ration d'un ensemble de point
y = yi(10001:11000);
%axe des abscisses
x=linspace(0,1000/Fs,1000);
%Calcul des coefficients d'autoR�gressions
a = lpc(y,2);
a2 = lpc(y);

%Calcul de la courbe approximant la fonction avec les coefficients
%d'autocorr�lation
est_x = filter([0 -a(2:end)],1,y);

est_x2 = filter([0 -a2(2:end)],1,y);

%Affichage
abs = linspace(0,10,10*Fs);
%plot(abs,yi(:,1))
plot(x,y,x,est_x,'--',x,est_x2,'.')