function [ ff, divf ] = func( x,y,z )
global radius
% -------------------------------------------------------------------------
rr=radius;
p=1;
q=1;
r=1;
func=(x.^p).*(y.^q).*(z.^r);
proj=(p+q+r)./rr.*func;

gr(:,:,1)=p.*x.^(p-1).*y.^q.*z.^r-proj.*x./rr;
gr(:,:,2)=q.*x.^p.*y.^(q-1).*z.^r-proj.*y./rr;
gr(:,:,3)=r.*x.^p.*y.^q.*z.^(r-1)-proj.*z./rr;

nx=x./rr;
ny=y./rr;
nz=z./rr;
ux=1;
uy=1;
uz=1;

divx=gr(:,:,1).*(ny.*uz-nz.*uy);
divy=gr(:,:,2).*(nz.*ux-nx.*uz);
divz=gr(:,:,3).*(nx.*uy-ny.*ux);
divf=divx+divy+divz;

ff(:,:,1)=func.*(ny.*uz-nz.*uy);
ff(:,:,2)=func.*(nz.*ux-nx.*uz);
ff(:,:,3)=func.*(nx.*uy-ny.*ux);
end