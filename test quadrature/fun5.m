function [v]=fun5(x)
% 1/x
global radius;
trsh=0.05; % SEUIL
polphi=[1/(trsh^6) 0 -3/trsh^4 0 3/(trsh^2) 0];
xtilde=x./radius;
axtilde=abs(xtilde);
axt1=max(axtilde,trsh);
v=(1/radius).*(polyval(polphi,xtilde).*(axtilde<=trsh)+ sign(xtilde).*(1./axt1).*(axtilde>trsh)); 