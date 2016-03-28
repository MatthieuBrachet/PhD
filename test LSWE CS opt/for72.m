function [mfor_I,mfor_II,mfor_III,mfor_IV,mfor_V,mfor_VI]=for72(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI,time,n,nn)
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;
%
global gprim hprim
global omega 
global radius
global kvit keta
% SOURCE TERM FOR EQUATIONS LTE%
%
% CASE SWITCH_LTE=1, cf aussi fun8.
% ------ face I
lambda=zeros(nn,nn);teta=zeros(nn,nn);radius1=zeros(nn,nn);
[lambda,teta,radius1]=cart2sph(x_fI,y_fI,z_fI);
%
elambda_x=zeros(nn,nn);elambda_y=zeros(nn,nn);elambda_z=zeros(nn,nn);
eteta_x=zeros(nn,nn);eteta_y=zeros(nn,nn);eteta_z=zeros(nn,nn);
%
elambda_x = -sin(lambda);
elambda_y =  cos(lambda);
elambda_z=0;
eteta_x = -sin(teta).*cos(lambda);
eteta_y = -sin(teta).*sin(lambda);
eteta_z =  cos(teta);
%
wex=fun8(x_fI,y_fI,z_fI,time);
v0=zeros(nn,nn,3);
v0=wex(:,:,1:3);
mfor_I(1:nn,1:nn,4)=0;
geta0=zeros(nn,nn,3);
geta0=-keta*sin(2*teta).*(eteta_x;eteta_y;eteta_z);
fclis=2*omega*sin(teta);
%
xnorm=zeros(nn,nn,3);
xnorm(:,:,1) = x_fI/radius;
xnorm(:,:,2) = y_fI/radius;
xnorm(:,:,3) = z_fI/radius;
%
mfor_I=zeros(nn,nn,4);
mfor_I(:,:,1:3)=gprim*geta0(:,:,1:3)+fclis.*cross(xnorm(:,:,1:3),v0(:,:,1:3));
mfor_I(1:nn,1:nn,4)=0;
% ------ face II
lambda=zeros(nn,nn);teta=zeros(nn,nn);radius1=zeros(nn,nn);
[lambda,teta,radius1]=cart2sph(x_fII,y_fII,z_fII);
%
elambda_x=zeros(nn,nn);elambda_y=zeros(nn,nn);elambda_z=zeros(nn,nn);
eteta_x=zeros(nn,nn);eteta_y=zeros(nn,nn);eteta_z=zeros(nn,nn);
%
elambda_x = -sin(lambda);
elambda_y =  cos(lambda);
elambda_z=0;
eteta_x = -sin(teta).*cos(lambda);
eteta_y = -sin(teta).*sin(lambda);
eteta_z =  cos(teta);
%
wex=fun8(x_fII,y_fII,z_fII,time);
v0=zeros(nn,nn,3);
v0=wex(:,:,1:3);
mfor_I(1:nn,1:nn,4)=0;
geta0=zeros(nn,nn,3);
geta0=-keta*sin(2*teta).*(eteta_x;eteta_y;eteta_z);
fclis=2*omega*sin(teta);
%
xnorm=zeros(nn,nn,3);
xnorm(:,:,1) = x_fII/radius;
xnorm(:,:,2) = y_fII/radius;
xnorm(:,:,3) = z_fII/radius;
%
mfor_II=zeros(nn,nn,4);
mfor_II(:,:,1:3)=gprim*geta0(:,:,1:3)+fclis.*cross(xnorm(:,:,1:3),v0(:,:,1:3));
mfor_II(1:nn,1:nn,4)=0;
% ------ face III
lambda=zeros(nn,nn);teta=zeros(nn,nn);radius1=zeros(nn,nn);
[lambda,teta,radius1]=cart2sph(x_fIII,y_fIII,z_fIII);
%
elambda_x=zeros(nn,nn);elambda_y=zeros(nn,nn);elambda_z=zeros(nn,nn);
eteta_x=zeros(nn,nn);eteta_y=zeros(nn,nn);eteta_z=zeros(nn,nn);
%
elambda_x = -sin(lambda);
elambda_y =  cos(lambda);
elambda_z=0;
eteta_x = -sin(teta).*cos(lambda);
eteta_y = -sin(teta).*sin(lambda);
eteta_z =  cos(teta);
%
wex=fun8(x_fIII,y_fIII,z_fIII,time);
v0=zeros(nn,nn,3);
v0=wex(:,:,1:3);
mfor_III(1:nn,1:nn,4)=0;
geta0=zeros(nn,nn,3);
geta0=-keta*sin(2*teta).*(eteta_x;eteta_y;eteta_z);
fclis=2*omega*sin(teta);
%
xnorm=zeros(nn,nn,3);
xnorm(:,:,1) = x_fIII/radius;
xnorm(:,:,2) = y_fIII/radius;
xnorm(:,:,3) = z_fIII/radius;
%
mfor_III=zeros(nn,nn,4);
mfor_III(:,:,1:3)=gprim*geta0(:,:,1:3)+fclis.*cross(xnorm(:,:,1:3),v0(:,:,1:3));
mfor_III(1:nn,1:nn,4)=0;
% ------ face IV
lambda=zeros(nn,nn);teta=zeros(nn,nn);radius1=zeros(nn,nn);
[lambda,teta,radius1]=cart2sph(x_fIV,y_fIV,z_fIV);
%
elambda_x=zeros(nn,nn);elambda_y=zeros(nn,nn);elambda_z=zeros(nn,nn);
eteta_x=zeros(nn,nn);eteta_y=zeros(nn,nn);eteta_z=zeros(nn,nn);
%
elambda_x = -sin(lambda);
elambda_y =  cos(lambda);
elambda_z=0;
eteta_x = -sin(teta).*cos(lambda);
eteta_y = -sin(teta).*sin(lambda);
eteta_z =  cos(teta);
%
wex=fun8(x_fIV,y_fIV,z_fIV,time);
v0=zeros(nn,nn,3);
v0=wex(:,:,1:3);
mfor_IV(1:nn,1:nn,4)=0;
geta0=zeros(nn,nn,3);
geta0=-keta*sin(2*teta).*(eteta_x;eteta_y;eteta_z);
fclis=2*omega*sin(teta);
%
xnorm=zeros(nn,nn,3);
xnorm(:,:,1) = x_fIV/radius;
xnorm(:,:,2) = y_fIV/radius;
xnorm(:,:,3) = z_fIV/radius;
%
mfor_IV=zeros(nn,nn,4);
mfor_IV(:,:,1:3)=gprim*geta0(:,:,1:3)+fclis.*cross(xnorm(:,:,1:3),v0(:,:,1:3));
mfor_IV(1:nn,1:nn,4)=0;
% ------ face V
lambda=zeros(nn,nn);teta=zeros(nn,nn);radius1=zeros(nn,nn);
[lambda,teta,radius1]=cart2sph(x_fV,y_fV,z_fV);
%
elambda_x=zeros(nn,nn);elambda_y=zeros(nn,nn);elambda_z=zeros(nn,nn);
eteta_x=zeros(nn,nn);eteta_y=zeros(nn,nn);eteta_z=zeros(nn,nn);
%
elambda_x = -sin(lambda);
elambda_y =  cos(lambda);
elambda_z=0;
eteta_x = -sin(teta).*cos(lambda);
eteta_y = -sin(teta).*sin(lambda);
eteta_z =  cos(teta);
%
wex=fun8(x_fV,y_fV,z_fV,time);
v0=zeros(nn,nn,3);
v0=wex(:,:,1:3);
mfor_V(1:nn,1:nn,4)=0;
geta0=zeros(nn,nn,3);
geta0=-keta*sin(2*teta).*(eteta_x;eteta_y;eteta_z);
fclis=2*omega*sin(teta);
%
xnorm=zeros(nn,nn,3);
xnorm(:,:,1) = x_fV/radius;
xnorm(:,:,2) = y_fV/radius;
xnorm(:,:,3) = z_fV/radius;
%
mfor_V=zeros(nn,nn,4);
mfor_V(:,:,1:3)=gprim*geta0(:,:,1:3)+fclis.*cross(xnorm(:,:,1:3),v0(:,:,1:3));
mfor_V(1:nn,1:nn,4)=0;
% ------ face VI
lambda=zeros(nn,nn);teta=zeros(nn,nn);radius1=zeros(nn,nn);
[lambda,teta,radius1]=cart2sph(x_fVI,y_fVI,z_fVI);
%
elambda_x=zeros(nn,nn);elambda_y=zeros(nn,nn);elambda_z=zeros(nn,nn);
eteta_x=zeros(nn,nn);eteta_y=zeros(nn,nn);eteta_z=zeros(nn,nn);
%
elambda_x = -sin(lambda);
elambda_y =  cos(lambda);
elambda_z=0;
eteta_x = -sin(teta).*cos(lambda);
eteta_y = -sin(teta).*sin(lambda);
eteta_z =  cos(teta);
%
wex=fun8(x_fVI,y_fVI,z_fVI,time);
v0=zeros(nn,nn,3);
v0=wex(:,:,1:3);
mfor_VI(1:nn,1:nn,4)=0;
geta0=zeros(nn,nn,3);
geta0=-keta*sin(2*teta).*(eteta_x;eteta_y;eteta_z);
fclis=2*omega*sin(teta);
%
xnorm=zeros(nn,nn,3);
xnorm(:,:,1) = x_fVI/radius;
xnorm(:,:,2) = y_fVI/radius;
xnorm(:,:,3) = z_fVI/radius;
%
mfor_VI=zeros(nn,nn,4);
mfor_VI(:,:,1:3)=gprim*geta0(:,:,1:3)+fclis.*cross(xnorm(:,:,1:3),v0(:,:,1:3));
mfor_VI(1:nn,1:nn,4)=0;
