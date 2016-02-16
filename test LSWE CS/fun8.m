function [v]=fun8(x,y,z,t)
% n=max(size(x,1),size(x,2));
% v=zeros(n,1);
% v1=zeros(n,1);v2=zeros(n,1);v3=zeros(n,1);
%
global nn
global radius
%
v=zeros(nn,nn,4);
[lambda,teta,radius1]=cart2sph(x,y,z);
%
%
elambda_x=zeros(nn,nn);elambda_y=zeros(nn,nn);elambda_z=zeros(nn,nn);
elambda_x = -sin(lambda);
elambda_y =  cos(lambda);
elambda_z=0;
%
eteta_x=zeros(nn,nn);eteta_y=zeros(nn,nn);eteta_z=zeros(nn,nn);
eteta_x = -sin(teta).*cos(lambda);
eteta_y = -sin(teta).*sin(lambda);
eteta_z =  cos(teta);
%
% case_lte='case_LTE1';
% --1--
% global kw Cw ku keta;
% kw=1;
% Cw=1;
% ku=1;
% keta=1;
% %
% eta0=zeros(nn,nn);
% u0=zeros(nn,nn); v0=zeros(nn,nn);
% eta0= keta*cos(teta).^2;
% u0=ku*(sin(0.5*lambda).^2).*sin(2*teta);
% v0=0.5*ku*sin(lambda).*cos(teta);
% %
% elambda_x=zeros(nn,nn);elambda_y=zeros(nn,nn);elambda_z=zeros(nn,nn);
% elambda_x = -sin(lambda);
% elambda_y =  cos(lambda);
% elambda_z=zeros(nn,nn);
% %
% eteta_x=zeros(nn,nn);eteta_y=zeros(nn,nn);eteta_z=zeros(nn,nn);
% eteta_x = -sin(teta).*cos(lambda);
% eteta_y = -sin(teta).*sin(lambda);
% eteta_z =  cos(teta);
% %
% for i=1:nn,
%     for j=1:nn,
%         v(i,j,1)=(u0(i,j)*elambda_x(i,j)+v0(i,j)*eteta_x(i,j))*cos(kw*(lambda(i,j)-Cw*t));
%         v(i,j,2)=(u0(i,j)*elambda_y(i,j)+v0(i,j)*eteta_y(i,j))*cos(kw*(lambda(i,j)-Cw*t));
%         v(i,j,3)=(u0(i,j)*elambda_z(i,j)+v0(i,j)*eteta_z(i,j))*cos(kw*(lambda(i,j)-Cw*t));
%         v(i,j,4)=eta0(i,j)*cos(kw*(lambda(i,j)-Cw*t));
%     end
% end
% --2--
 case_lte='case_LTE1';
n1=nn;
n2=nn;
global kv keta;
kv=1;keta=1;
kw1=1;kw2=1;kw3=1;
cv=kv*(sin(kw1*pi*x/radius)+sin(kw2*pi*y/radius)+sin(kw3*pi*z/radius));
gcv(1:n1,1:n2,1)=kv*kw1*(pi/radius)*cos(kw1*pi*x/radius);
gcv(1:n1,1:n2,2)=kv*kw2*(pi/radius)*cos(kw2*pi*y/radius);
gcv(1:n1,1:n2,3)=kv*kw3*(pi/radius)*cos(kw3*pi*z/radius);
 %
xnorm=zeros(n1,n2,3);
nnorm=sqrt(x.*x+y.*y+z.*z);
xnorm(1:n1,1:n2,1)=x./nnorm;
xnorm(1:n1,1:n2,2)=y./nnorm;
xnorm(1:n1,1:n2,3)=z./nnorm;
% SPHERICAL GRADIENT 
v(1:n1,1:n2,1)=gcv(1:n1,1:n2,1)...
                   -dot(xnorm(1:n1,1:n2,1:3),gcv(1:n1,1:n2,1:3),3).*xnorm(1:n1,1:n2,1);
v(1:n1,1:n2,2)=gcv(1:n1,1:n2,2)...
                   -dot(xnorm(1:n1,1:n2,1:3),gcv(1:n1,1:n2,1:3),3).*xnorm(1:n1,1:n2,2);
v(1:n1,1:n2,3)=gcv(1:n1,1:n2,3)...
                   -dot(xnorm(1:n1,1:n2,1:3),gcv(1:n1,1:n2,1:3),3).*xnorm(1:n1,1:n2,3);
%
v(1:n1,1:n2,4)=keta*(sin(kw1*pi*x/radius)+sin(kw2*pi*y/radius)+sin(kw3*pi*z/radius));


    