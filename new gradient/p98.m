% MODULE PROBLEM FOR THE CUBED SPHERE
% ----------------------------------
clc; close all; clear all;

global n nn;
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;
global scheme

scheme='compact4';


n=255;
mod98; 

[funfI,funx,funy,funz]=fun3(x_fI,y_fI,z_fI);
[funfII,funx,funy,funz]=fun3(x_fII,y_fII,z_fII);
[funfIII,funx,funy,funz]=fun3(x_fIII,y_fIII,z_fIII);
[funfIV,funx,funy,funz]=fun3(x_fIV,y_fIV,z_fIV);
[funfV,funx,funy,funz]=fun3(x_fV,y_fV,z_fV);
[funfVI,funx,funy,funz]=fun3(x_fVI,y_fVI,z_fVI);

tic
for ite=1:50,
    ite
    [grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
      gr98(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn);
end
toc

figure(1);
plot(grad_I(:,11,2),'r'); grid;

%% FACE I
% gradient exact
[func,func_x,func_y,func_z]=fun3(x_fI,y_fI,z_fI);
% projection du gradient sur la sphere
norm_I=zeros(3,1);
gradse_I=zeros(nn,nn,3);gradset_I=zeros(nn,nn,3);
gradse_I(1:nn,1:nn,1)=func_x ; 
gradse_I(1:nn,1:nn,2)=func_y ; 
gradse_I(1:nn,1:nn,3)=func_z;
norm_xI=zeros(nn,nn);norm_yI=zeros(nn,nn);norm_zI=zeros(nn,nn);
for i=1:nn,
    for j=1:nn,
        norm_xI(i,j)=x_fI(i,j) ; norm_yI(i,j)=y_fI(i,j) ; norm_zI(i,j)=z_fI(i,j);
        xwk(1)=gradse_I(i,j,1) ; xwk(2)=gradse_I(i,j,2) ; xwk(3)=gradse_I(i,j,3);
        gradset_I(i,j,1)=gradse_I(i,j,1)-(xwk(1)*norm_xI(i,j)+xwk(2)*norm_yI(i,j)+xwk(3)*norm_zI(i,j))*norm_xI(i,j);
        gradset_I(i,j,2)=gradse_I(i,j,2)-(xwk(1)*norm_xI(i,j)+xwk(2)*norm_yI(i,j)+xwk(3)*norm_zI(i,j))*norm_yI(i,j);
        gradset_I(i,j,3)=gradse_I(i,j,3)-(xwk(1)*norm_xI(i,j)+xwk(2)*norm_yI(i,j)+xwk(3)*norm_zI(i,j))*norm_zI(i,j);
    end
end
% calcul erreur sur le gradient en X:
errg_fI=zeros(nn,nn);
for i=1:nn,
    for j=1:nn,
        errg_fI(i,j)=max(abs(grad_I(i,j,:)-gradset_I(i,j,:)));
    end
end
errI_15=max(max(abs(errg_fI)));

%% FACE II
% gradient exact
[func,func_x,func_y,func_z]=fun3(x_fII,y_fII,z_fII);
% projection du gradient sur la sphere
norm_II=zeros(3,1);
gradse_II=zeros(nn,nn,3);gradset_II=zeros(nn,nn,3); 
gradse_II(1:nn,1:nn,1)=func_x ; gradse_II(1:nn,1:nn,2)=func_y ; gradse_II(1:nn,1:nn,3)=func_z;
norm_xII=zeros(nn,nn);norm_yII=zeros(nn,nn);norm_zII=zeros(nn,nn);
for i=1:nn,
    for j=1:nn,
        norm_xII(i,j)=x_fII(i,j) ; norm_yII(i,j)=y_fII(i,j) ; norm_zII(i,j)=z_fII(i,j);
        xwk(1)=gradse_II(i,j,1) ; xwk(2)=gradse_II(i,j,2) ; xwk(3)=gradse_II(i,j,3);
        gradset_II(i,j,1)=gradse_II(i,j,1)-(xwk(1)*norm_xII(i,j)+xwk(2)*norm_yII(i,j)+xwk(3)*norm_zII(i,j))*norm_xII(i,j);
        gradset_II(i,j,2)=gradse_II(i,j,2)-(xwk(1)*norm_xII(i,j)+xwk(2)*norm_yII(i,j)+xwk(3)*norm_zII(i,j))*norm_yII(i,j);
        gradset_II(i,j,3)=gradse_II(i,j,3)-(xwk(1)*norm_xII(i,j)+xwk(2)*norm_yII(i,j)+xwk(3)*norm_zII(i,j))*norm_zII(i,j);
    end
end
% calcul erreur sur le gradient en X:
errg_fII=zeros(nn,nn);
for i=1:nn,
    for j=1:nn,
        errg_fII(i,j)=max(abs(grad_II(i,j,:)-gradset_II(i,j,:)));
    end
end
errII_16=max(max(abs(errg_fII)));

%% FACE III
% gradient exact
[func,func_x,func_y,func_z]=fun3(x_fIII,y_fIII,z_fIII);
% projection du gradient sur la sphere
norm_III=zeros(3,1);
gradse_III=zeros(nn,nn,3);gradset_III=zeros(nn,nn,3);
gradse_III(1:nn,1:nn,1)=func_x ; gradse_III(1:nn,1:nn,2)=func_y ; gradse_III(1:nn,1:nn,3)=func_z;
norm_xIII=zeros(nn,nn);norm_yIII=zeros(nn,nn);norm_zIII=zeros(nn,nn);
for i=1:nn,
    for j=1:nn,
        norm_xIII(i,j)=x_fIII(i,j) ; norm_yIII(i,j)=y_fIII(i,j) ; norm_zIII(i,j)=z_fIII(i,j);
        xwk(1)=gradse_III(i,j,1) ; xwk(2)=gradse_III(i,j,2) ; xwk(3)=gradse_III(i,j,3);
        gradset_III(i,j,1)=gradse_III(i,j,1)-(xwk(1)*norm_xIII(i,j)+xwk(2)*norm_yIII(i,j)+xwk(3)*norm_zIII(i,j))*norm_xIII(i,j);
        gradset_III(i,j,2)=gradse_III(i,j,2)-(xwk(1)*norm_xIII(i,j)+xwk(2)*norm_yIII(i,j)+xwk(3)*norm_zIII(i,j))*norm_yIII(i,j);
        gradset_III(i,j,3)=gradse_III(i,j,3)-(xwk(1)*norm_xIII(i,j)+xwk(2)*norm_yIII(i,j)+xwk(3)*norm_zIII(i,j))*norm_zIII(i,j);
    end
end
% calcul erreur sur le gradient
errg_fIII=zeros(nn,nn);
for i=1:nn,
    for j=1:nn,
        errg_fIII(i,j)=max(abs(grad_III(i,j,:)-gradset_III(i,j,:)));
    end
end
errIII_17=max(max(abs(errg_fIII)));

%% FACE IV
% gradient exact
[func,func_x,func_y,func_z]=fun3(x_fIV,y_fIV,z_fIV);
% projection du gradient sur la sphere
norm_IV=zeros(3,1);
gradse_IV=zeros(nn,nn,3);gradset_IV=zeros(nn,nn,3); 
gradse_IV(1:nn,1:nn,1)=func_x ; gradse_IV(1:nn,1:nn,2)=func_y ; gradse_IV(1:nn,1:nn,3)=func_z;
norm_xIV=zeros(nn,nn);norm_yIV=zeros(nn,nn);norm_zIV=zeros(nn,nn);
for i=1:nn,
    for j=1:nn,
        norm_xIV(i,j)=x_fIV(i,j) ; norm_yIV(i,j)=y_fIV(i,j) ; norm_zIV(i,j)=z_fIV(i,j);
        xwk(1)=gradse_IV(i,j,1) ; xwk(2)=gradse_IV(i,j,2) ; xwk(3)=gradse_IV(i,j,3);
        gradset_IV(i,j,1)=gradse_IV(i,j,1)-(xwk(1)*norm_xIV(i,j)+xwk(2)*norm_yIV(i,j)+xwk(3)*norm_zIV(i,j))*norm_xIV(i,j);
        gradset_IV(i,j,2)=gradse_IV(i,j,2)-(xwk(1)*norm_xIV(i,j)+xwk(2)*norm_yIV(i,j)+xwk(3)*norm_zIV(i,j))*norm_yIV(i,j);
        gradset_IV(i,j,3)=gradse_IV(i,j,3)-(xwk(1)*norm_xIV(i,j)+xwk(2)*norm_yIV(i,j)+xwk(3)*norm_zIV(i,j))*norm_zIV(i,j);
    end
end
% calcul erreur sur le gradient
errg_fIV=zeros(nn,nn);
for i=1:nn,
    for j=1:nn,
        errg_fIV(i,j)=max(abs(grad_IV(i,j,:)-gradset_IV(i,j,:)));
    end
end
errIV_18=max(max(abs(errg_fIV)));

%% FACE V
% gradient exact
[func,func_x,func_y,func_z]=fun3(x_fV,y_fV,z_fV);
% projection du gradient sur la sphere
norm_V=zeros(3,1);
gradse_V=zeros(nn,nn,3);gradset_V=zeros(nn,nn,3); 
gradse_V(1:nn,1:nn,1)=func_x ; gradse_V(1:nn,1:nn,2)=func_y ; gradse_V(1:nn,1:nn,3)=func_z;
norm_xV=zeros(nn,nn);norm_yV=zeros(nn,nn);norm_zV=zeros(nn,nn);
for i=1:nn,
    for j=1:nn,
        norm_xV(i,j)=x_fV(i,j) ; norm_yV(i,j)=y_fV(i,j) ; norm_zV(i,j)=z_fV(i,j);
        xwk(1)=gradse_V(i,j,1) ; xwk(2)=gradse_V(i,j,2) ; xwk(3)=gradse_V(i,j,3);
        gradset_V(i,j,1)=gradse_V(i,j,1)-(xwk(1)*norm_xV(i,j)+xwk(2)*norm_yV(i,j)+xwk(3)*norm_zV(i,j))*norm_xV(i,j);
        gradset_V(i,j,2)=gradse_V(i,j,2)-(xwk(1)*norm_xV(i,j)+xwk(2)*norm_yV(i,j)+xwk(3)*norm_zV(i,j))*norm_yV(i,j);
        gradset_V(i,j,3)=gradse_V(i,j,3)-(xwk(1)*norm_xV(i,j)+xwk(2)*norm_yV(i,j)+xwk(3)*norm_zV(i,j))*norm_zV(i,j);
    end
end
% calcul erreur sur le gradient en X:
errg_fV=zeros(nn,nn);
for i=1:nn,
    for j=1:nn,
        errg_fV(i,j)=max(abs(grad_V(i,j,:)-gradset_V(i,j,:)));
    end
end
errV_19=max(max(abs(errg_fV)));

%% FACE VI 
% gradient exact
[func,func_x,func_y,func_z]=fun3(x_fVI,y_fVI,z_fVI);
% projection du gradient sur la sphere
norm_VI=zeros(3,1);
gradse_VI=zeros(nn,nn,3);gradset_VI=zeros(nn,nn,3); 
gradse_VI(1:nn,1:nn,1)=func_x ; gradse_VI(1:nn,1:nn,2)=func_y ; gradse_VI(1:nn,1:nn,3)=func_z;
norm_xVI=zeros(nn,nn);norm_yVI=zeros(nn,nn);norm_zVI=zeros(nn,nn);
for i=1:nn,
    for j=1:nn,
        norm_xVI(i,j)=x_fVI(i,j) ; norm_yVI(i,j)=y_fVI(i,j) ; norm_zVI(i,j)=z_fVI(i,j);
        xwk(1)=gradse_VI(i,j,1) ; xwk(2)=gradse_VI(i,j,2) ; xwk(3)=gradse_VI(i,j,3);
        gradset_VI(i,j,1)=gradse_VI(i,j,1)-(xwk(1)*norm_xVI(i,j)+xwk(2)*norm_yVI(i,j)+xwk(3)*norm_zVI(i,j))*norm_xVI(i,j);
        gradset_VI(i,j,2)=gradse_VI(i,j,2)-(xwk(1)*norm_xVI(i,j)+xwk(2)*norm_yVI(i,j)+xwk(3)*norm_zVI(i,j))*norm_yVI(i,j);
        gradset_VI(i,j,3)=gradse_VI(i,j,3)-(xwk(1)*norm_xVI(i,j)+xwk(2)*norm_yVI(i,j)+xwk(3)*norm_zVI(i,j))*norm_zVI(i,j);
    end
end
% calcul erreur sur le gradient
errg_fVI=zeros(nn,nn);
for i=1:nn,
    for j=1:nn,
        errg_fVI(i,j)=max(abs(grad_VI(i,j,:)-gradset_VI(i,j,:)));
    end
end
errVI_20=max(max(abs(errg_fVI)));








err_grad_21=max([errVI_20,errV_19,errIV_18,errIII_17,errII_16,errI_15]);
for i=1:nn,
    for j=1:nn,
        errg_fI(i,j)=max(abs(grad_I(i,j,:)-gradset_I(i,j,:)));
    end
end
for i=1:nn,
    for j=1:nn,
        errg_fII(i,j)=max(abs(grad_II(i,j,:)-gradset_II(i,j,:)));
    end
end
for i=1:nn,
    for j=1:nn,
        errg_fIII(i,j)=max(abs(grad_III(i,j,:)-gradset_III(i,j,:)));
    end
end
for i=1:nn,
    for j=1:nn,
        errg_fIV(i,j)=max(abs(grad_IV(i,j,:)-gradset_IV(i,j,:)));
    end
end
for i=1:nn,
    for j=1:nn,
        errg_fV(i,j)=max(abs(grad_V(i,j,:)-gradset_V(i,j,:)));
    end
end
for i=1:nn,
    for j=1:nn,
        errg_fVI(i,j)=max(abs(grad_VI(i,j,:)-gradset_VI(i,j,:)));
    end
end
errI_25=max(max(abs(errg_fI)));
errII_26=max(max(abs(errg_fII)));
errIII_27=max(max(abs(errg_fIII)));
errIV_28=max(max(abs(errg_fIV)));
errV_29=max(max(abs(errg_fV)));
errVI_30=max(max(abs(errg_fVI)));
err_grad_31=max([errVI_30,errV_29,errIV_28,errIII_27,errII_26,errI_25])

  
