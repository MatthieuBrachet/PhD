function [v]=fun7(x)
% 1/abs(x)
global radius;
trsh=0.05;% SEUIL
polpsi=[3/(8*trsh^5) 0 -5/(4*trsh^3) 0 15/(8*trsh)];
xtilde=x./radius;
axtilde=abs(xtilde);
axt1=max(axtilde,trsh);
v=1./radius*(polyval(polpsi,xtilde).*(axtilde<=trsh)+ (1./axt1).*(axtilde>trsh)); 