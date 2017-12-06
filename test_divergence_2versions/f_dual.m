function [h,u] = f_dual( x,y,z )
global radius
xx=x/radius; yy=y/radius; zz=z/radius;
h=sin(pi*xx)+sin(2*pi*yy)+sin(3*pi*zz);

nx=x./radius;
ny=y./radius;
nz=z./radius;
ux=1;
uy=2;
uz=3;
func=exp(xx)+exp(yy)+exp(zz);

u(:,:,1)=func.*(ny.*uz-nz.*uy);
u(:,:,2)=func.*(nz.*ux-nx.*uz);
u(:,:,3)=func.*(nx.*uy-ny.*ux);
end