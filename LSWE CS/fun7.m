function [v]=fun7(x)
polpsi=[12 0 -10 0 15/4];
ax=abs(x);
ax1=max(ax,0.5);
v=polyval(polpsi,x).*(ax<=0.5)+ (1./ax1).*(ax>0.5); 

