function [v]=fun5(x)
polphi=[64 0 -48 0 12 0];
ax=abs(x);
ax1=max(ax,0.5);
v=polyval(polphi,x).*(ax<=0.5)+ sign(x).*(1./ax1).*(ax>0.5); 