function [grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI,va1,va2]=...
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

%% SET 1 : 
va_fI=zeros(4*(nn-1),nn); 

% PANEL I :
va_fI(1:nn-1,1:nn)=funfI(1:nn-1,1:nn);

% Panel II:
funspline(1:nn,1:nn-1)=funfII(1:nn-1,1:nn)';
sm1=(3/deta)*kmat*funspline(2:nn-1,1:nn-1);
sm1(1,1:nn-1)=sm1(1,1:nn-1)+(m1/(2*deta))*funspline(1,1:nn-1);
sm1(n,1:nn-1)=sm1(n,1:nn-1)+(m2/(2*deta))*funspline(nn,1:nn-1);
funspline_x(2:nn-1,1:nn-1)=umat1\(lmat1\sm1);
funspline_x(1,1:nn-1)=(1/alfasp)*((betasp/deta)*(funspline(2,1:nn-1)-funspline(1,1:nn-1))... % U_x,0
    +(gamasp/(2*deta)*(funspline(3,1:nn-1)-funspline(1,1:nn-1)))...
    -betasp*funspline_x(2,1:nn-1)-gamasp*funspline_x(3,1:nn-1));
funspline_x(nn,1:nn-1)=(1/alfasp)*(-(betasp/alfasp)*(funspline(nn-1,1:nn-1)-funspline(nn,1:nn-1))... % U_x,N
    -(gamasp/(2*deta))*(funspline(nn-2,1:nn-1)-funspline(nn,1:nn-1))...
    +betasp*funspline_x(nn-1,1:nn-1)+gamasp*funspline_x(nn-2,1:nn-1));
% COEFFTS C_0,C_1, C_2 AND C_3 OF CUBIC SPLINE J+1/2
xc0(1:nn,1:nn-1)=funspline(1:nn,1:nn-1);
xc1(1:nn,1:nn-1)=funspline_x(1:nn,1:nn-1);
xc2(1:nn-1,1:nn-1)=3*(funspline(2:nn,1:nn-1)-funspline(1:nn-1,1:nn-1))/deta^2 ...
            -(2*funspline_x(1:nn-1,1:nn-1)+funspline_x(2:nn,1:nn-1))/deta;
xc2(nn,1:nn-1)=0;        
xc3(1:nn-1,1:nn-1)=2*(funspline(1:nn-1,1:nn-1)-funspline(2:nn,1:nn-1))/deta^3 ...
            + (funspline_x(1:nn-1,1:nn-1)+funspline_x(2:nn,1:nn-1))/deta^2;
xc3(nn,1:nn-1)=0;
jj=jeta_c(1:nn-1,[1:nn]);
xeta1=eta_c(1:nn-1,[1:nn]);
xeta2=xeta1-eta(jj);
xeta2=xeta2';
funbII1=((xc3(jj,1:nn-1).*xeta2+xc2(jj,1:nn-1)).*xeta2+xc1(jj,1:nn-1)).*xeta2+xc0(jj,1:nn-1);
va_fI(nn-1+[1:nn-1],1:nn)=funbII1(1:nn-1,1:nn)';
va1=va_fI;



for i=1:nn-1, 
    funspline(1:nn)=funfII(i,1:nn);
    sm1=(3/deta)*kmat*funspline(2:nn-1);
    sm1(1)=sm1(1)+(m1/(2*deta))*funspline(1);
    sm1(n)=sm1(n)+(m2/(2*deta))*funspline(nn);
    funspline_x=zeros(nn,1);
    funspline_x(2:nn-1)=umat1\(lmat1\sm1);
    funspline_x(1)=(1/alfasp)*((betasp/deta)*(funspline(2)-funspline(1))... % U_x,0
        +(gamasp/(2*deta)*(funspline(3)-funspline(1)))...
        -betasp*funspline_x(2)-gamasp*funspline_x(3));
    funspline_x(nn)=(1/alfasp)*(-(betasp/alfasp)*(funspline(nn-1)-funspline(nn))... % U_x,N
        -(gamasp/(2*deta))*(funspline(nn-2)-funspline(nn))...
        +betasp*funspline_x(nn-1)+gamasp*funspline_x(nn-2));
    
     % COEFFTS C_0,C_1, C_2 AND C_3 OF CUBIC SPLINE J+1/2
    xc0(1:nn)=funspline(1:nn);
    xc1(1:nn)=funspline_x(1:nn);
    xc2(1:nn-1)=3*(funspline(2:nn)-funspline(1:nn-1))/deta^2 ...
                -(2*funspline_x(1:nn-1)+funspline_x(2:nn))/deta;
    xc2(nn)=0;        
    xc3(1:nn-1)=2*(funspline(1:nn-1)-funspline(2:nn))/deta^3 ...
                + (funspline_x(1:nn-1)+funspline_x(2:nn))/deta^2;
    xc3(nn)=0;
    jj=jeta_c(i,[1:nn]);
    xeta1=eta_c(i,[1:nn]);
    xeta2=xeta1-eta(jj);
    xeta2=xeta2';
    funbII1=((xc3(jj).*xeta2+xc2(jj)).*xeta2+xc1(jj)).*xeta2+xc0(jj);
    va_fI(nn-1+i,1:nn)=funbII1(1:nn);
end
va2=va_fI;



