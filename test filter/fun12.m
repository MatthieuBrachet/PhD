function [y] = fun12(x)
global radius omega u0
xe=0.3;
tetab=-pi/6;
tetae=pi/2;
[lambda, teta,r]=cart2sph(x,y,z);
[ lambdap, tetap ] = rotated_coord( lambda, teta );
x=xe.*(tetap-tetab).*(tetae-tetab).^(-1);
bx=exp(-1./x).*(x<0);
up=u0.*bx.*exp(-1./(xe-x)).*((xe-x)<0);
y=(2.*omega.*sin(x)+up.*tan(x)./radius).*up;