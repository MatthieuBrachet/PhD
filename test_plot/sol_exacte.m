function [ht,vt] = sol_exacte(x,y,z)
% exact solution for Williamson test.
global radius

gp=10;
h0=6500;
omega=7*10^-5;

% Rossby-Haurwitz wave (test case 6 of Williamson & al.).
a=radius;
R=4;
w=7.848*10^-6;
K=w;

[n1,n2]=size(x);
vt=zeros(n1,n2,3);
[lambda, teta,r]=cart2sph(x,y,z);
uu=a.*w.*cos(teta)+a.*K.*cos(teta).^(R-1).*(R.*sin(teta).^2-cos(teta).^2).*cos(R.*lambda);
vv=-a.*K.*R.*cos(teta).^(R-1).*sin(teta).*sin(R.*lambda);

elambda_x = -sin(lambda);
elambda_y =  cos(lambda);
elambda_z=zeros(size(x));
eteta_x = -sin(teta).*cos(lambda);
eteta_y = -sin(teta).*sin(lambda);
eteta_z =  cos(teta);

vt(:,:,1)=uu.*elambda_x + vv.*eteta_x;
vt(:,:,2)=uu.*elambda_y + vv.*eteta_y;
vt(:,:,3)=uu.*elambda_z + vv.*eteta_z;

A=(w/2).*(2*omega+w).*cos(teta).^2+.25*K^2.*cos(teta).^(2*R).*((R+1).*cos(teta).^2+(2*R^2-R-2)-2*R^2.*cos(teta).^(-2));
B=((2*(omega+w)*K)./((R+1).*(R+2))).*cos(teta).^R.*((R^2+2*R+2)-((R+1).^2).*cos(teta).^2);
C=.25*K^2*cos(teta).^(2*R).*((R+1).*cos(teta).^2-(R+2));

ght=gp*h0+a.^2.*A+a.^2.*B.*cos(R.*lambda)+a.^2.*C.*cos(2*R*lambda);
ht=ght./gp;

end

