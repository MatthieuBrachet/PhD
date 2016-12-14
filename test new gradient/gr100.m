function [grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI,vb_fI]=...
    gr100(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn)
global xi eta dxi deta
global beta betacr;
global pg kg;
global gxi_I gxi_II gxi_III gxi_IV gxi_V gxi_VI;
global geta_I geta_II geta_III geta_IV geta_V geta_VI;
global kmat
global eta_c jeta_c xi_c ixi_c
global alfasp betasp gamasp
global m1 m2
global umat1 lmat1

grad_I=zeros(nn,nn,3);grad_II=zeros(nn,nn,3);grad_III=zeros(nn,nn,3);
grad_IV=zeros(nn,nn,3);grad_V=zeros(nn,nn,3);grad_VI=zeros(nn,nn,3);
funspline=zeros(nn,1);
xc0=zeros(nn,1);xc1=zeros(nn,1);xc2=zeros(nn,1);xc3=zeros(nn,1);

%% RESEAU 1 : 
va_fI=zeros(4*(nn-1),nn); 

% PANEL I :
va_fI(1:nn-1,1:nn)=funfI(1:nn-1,1:nn);

% PANEL II: 
fun=funfII(1:nn-1,1:nn)';
[xc0,xc1,xc2,xc3] = spl100(fun);
for i=1:nn-1
    jj=jeta_c(i,[1:nn]);
    xeta1=eta_c(i,[1:nn]);
    xeta2=xeta1-eta(jj);
    xeta2=xeta2';
    funbII1=((xc3(jj,i).*xeta2+xc2(jj,i)).*xeta2+xc1(jj,i)).*xeta2+xc0(jj,i);
    va_fI(nn-1+i,1:nn)=funbII1(1:nn);
end

% PANEL III
va_fI(2*nn-1:3*nn-3,1:nn)=funfIII(1:nn-1,nn:-1:1);

%  PANEL IV:
fun=funfIV(1:nn-1,1:nn)';
[xc0,xc1,xc2,xc3] = spl100(fun);
for i=1:nn-1, 
    jj=jeta_c(i,nn+1-[1:nn]);
    xeta1=eta_c(i,nn+1-[1:nn]);
    xeta2=xeta1-eta(jj);
    xeta2=xeta2';
    funbIV1=((xc3(jj,i).*xeta2+xc2(jj,i)).*xeta2+xc1(jj,i)).*xeta2+xc0(jj,i);

    va_fI(3*nn-3+i,1:nn)=funbIV1(1:nn);
end

funaI=(3/dxi)*va_fI(:,1:nn);
xsmIa=kg*funaI;
vad_fI(:,1:nn)=pg\xsmIa;

%% SET 2:

vb_fI=zeros(nn,4*(nn-1));
% PANEL I :
vb_fI(1:nn,1:nn-1)=funfI(1:nn,1:nn-1);

% FACE V:
% for j=1:nn-1, 
%     funspline(1:nn)=funfV(1:nn,j);
%     sm1=(3/deta)*kmat*funspline(2:nn-1);
%     sm1(1)=sm1(1)+(m1/(2*deta))*funspline(1);
%     sm1(n)=sm1(n)+(m2/(2*deta))*funspline(nn);
%     funspline_x=zeros(nn,1);
%     funspline_x(2:nn-1)=umat1\(lmat1\sm1);
%     funspline_x(1)=(1/alfasp)*((betasp/deta)*(funspline(2)-funspline(1))... 
%         +(gamasp/(2*deta)*(funspline(3)-funspline(1)))...
%         -betasp*funspline_x(2)-gamasp*funspline_x(3));
%     funspline_x(nn)=(1/alfasp)*(-(betasp/alfasp)*(funspline(nn-1)-funspline(nn))...
%         -(gamasp/(2*deta))*(funspline(nn-2)-funspline(nn))...
%         +betasp*funspline_x(nn-1)+gamasp*funspline_x(nn-2));
% 
%      % COEFFTS C_0,C_1, C_2 AND C_3 OF CUBIC SPLINE J+1/2
%     xc0(1:nn)=funspline(1:nn);
%     xc1(1:nn)=funspline_x(1:nn);
%     xc2(1:nn-1)=3*(funspline(2:nn)-funspline(1:nn-1))/deta^2 ...
%                 -(2*funspline_x(1:nn-1)+funspline_x(2:nn))/deta;
%     xc2(nn)=0;        
%     xc3(1:nn-1)=2*(funspline(1:nn-1)-funspline(2:nn))/deta^3 ...
%                 + (funspline_x(1:nn-1)+funspline_x(2:nn))/deta^2;
%     xc3(nn)=0;

fun=funfV(1:nn,1:nn-1);
[xc0,xc1,xc2,xc3] = spl100(fun);
for j=1:nn-1
    ii=ixi_c([1:nn],j);
    xxi1=xi_c([1:nn],j);
    xxi2=xxi1'-xi(ii);
    xxi2=xxi2';
    funaV1=((xc3(ii).*xxi2+xc2(ii)).*xxi2+xc1(ii)).*xxi2+xc0(ii);

    vb_fI(1:nn,nn-1+j)=funaV1(1:nn);
end


% FACE III:
vb_fI(1:nn,2*nn-1:3*nn-3)=funfIII(1:nn,nn:-1:2); 
 
% FACE VI: 
for j=1:nn-1,
    funspline(1:nn)=funfVI(1:nn,j); 
     %
    sm1=(3/dxi)*kmat*funspline(2:nn-1);
    sm1(1)=sm1(1)+(m1/(2*dxi))*funspline(1);
    sm1(n)=sm1(n)+(m2/(2*dxi))*funspline(nn); 
    %
    funspline_x=zeros(nn,1);

    funspline_x(2:nn-1)=umat1\(lmat1\sm1);

    funspline_x(1)=(1/alfasp)*((betasp/dxi)*(funspline(2)-funspline(1))... % U_x,0
        +(gamasp/(2*dxi)*(funspline(3)-funspline(1)))...
        - betasp*funspline_x(2)-gamasp*funspline_x(3));
    funspline_x(nn)=(1/alfasp)*(-(betasp/alfasp)*(funspline(nn-1)-funspline(nn))... % U_x,N
        -(gamasp/(2*dxi))*(funspline(nn-2)-funspline(nn))...
        + betasp*funspline_x(nn-1)+gamasp*funspline_x(nn-2));

     % COEFFTS C_0,C_1, C_2 AND C_3 OF CUBIC SPLINE J+1/2
    xc0(1:nn)=funspline(1:nn);
    xc1(1:nn)=funspline_x(1:nn);
    xc2(1:nn-1)=3*(funspline(2:nn)-funspline(1:nn-1))/dxi^2 ...
                -(2*funspline_x(1:nn-1)+funspline_x(2:nn))/dxi;
    xc2(nn)=0;        
    xc3(1:nn-1)=2*(funspline(1:nn-1)-funspline(2:nn))/dxi^3 ...
                + (funspline_x(1:nn-1)+funspline_x(2:nn))/dxi^2;
    xc3(nn)=0;
    
    ii=ixi_c(nn+1-[1:nn],j);
    xxi1=xi_c(nn+1-[1:nn],j);
    xxi2=xxi1'-xi(ii);
    xxi2=xxi2';
    funaVI1=((xc3(ii).*xxi2+xc2(ii)).*xxi2+xc1(ii)).*xxi2+xc0(ii);

    vb_fI(1:nn,3*nn-3+j)=funaVI1(1:nn);
end

% CALCUL DES DERIVEES BETA SUR LE RESEAU DE GRANDS CERCLES I-BETA
funbI=(3/deta)*vb_fI(1:nn,:)';
xsmIb=kg*funbI;
vbd_fI(1:nn,:)=(pg\xsmIb)';
