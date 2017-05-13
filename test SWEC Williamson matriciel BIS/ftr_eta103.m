function [funftI,funftII,funftIII,funftIV,funftV,funftVI]=...
    ftr_eta103(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn)

global xi eta dxi deta
global ftr
global kmat
global eta_c jeta_c xi_c ixi_c
global alfasp betasp gamasp
global m1 m2
global umat1 lmat1

betaspline=zeros(nn,1);
alfaspline=zeros(nn,1);
funspline=zeros(nn,1);
funftI(1:nn,1:nn)=funfI(1:nn,1:nn);
funftII(1:nn,1:nn)=funfII(1:nn,1:nn);
funftIII(1:nn,1:nn)=funfIII(1:nn,1:nn);
funftIV(1:nn,1:nn)=funfIV(1:nn,1:nn);
funftV(1:nn,1:nn)=funfV(1:nn,1:nn);
funftVI(1:nn,1:nn)=funfVI(1:nn,1:nn);

%% RESEAU 1 : 

va_fI=zeros(4*(nn-1),nn); 

% PANEL I :
va_fI(1:nn-1,1:nn)=funfI(1:nn-1,1:nn);

% PANEL II :
fun=funfII(1:nn-1,1:nn)';
funspline(1:nn,1:nn-1)=fun;

sm1=(3/deta)*kmat*funspline(2:nn-1,1:nn-1);
sm1(1,1:nn-1)=sm1(1,1:nn-1)+(m1/(2*deta))*funspline(1,1:nn-1);
sm1(n,1:nn-1)=sm1(n,1:nn-1)+(m2/(2*deta))*funspline(nn,1:nn-1);

funspline_x(2:nn-1,1:nn-1)=umat1\(lmat1\sm1);
funspline_x(1,1:nn-1)=(1/alfasp)*((betasp/deta)*(funspline(2,1:nn-1)-funspline(1,1:nn-1))... % U_x,0
    +(gamasp/(2*deta)*(funspline(3,1:nn-1)-funspline(1,1:nn-1)))...
    -betasp*funspline_x(2,1:nn-1)-gamasp*funspline_x(3,1:nn-1));
funspline_x(nn,1:nn-1)=(1/alfasp)*(-(betasp/deta)*(funspline(nn-1,1:nn-1)-funspline(nn,1:nn-1))... % U_x,N
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

for i=1:nn-1
    jj=jeta_c(i,1:nn);
    xeta1=eta_c(i,1:nn);
    xeta2=xeta1-eta(jj);
    xeta2=xeta2';
    funbII1=((xc3(jj,i).*xeta2+xc2(jj,i)).*xeta2+xc1(jj,i)).*xeta2+xc0(jj,i);
    va_fI(nn-1+i,1:nn)=funbII1(1:nn);
end

% PANEL III
va_fI(2*nn-1:3*nn-3,1:nn)=funfIII(1:nn-1,nn:-1:1);

%  PANEL IV:
fun=funfIV(1:nn-1,1:nn)';
funspline(1:nn,1:nn-1)=fun;

sm1=(3/deta)*kmat*funspline(2:nn-1,1:nn-1);
sm1(1,1:nn-1)=sm1(1,1:nn-1)+(m1/(2*deta))*funspline(1,1:nn-1);
sm1(n,1:nn-1)=sm1(n,1:nn-1)+(m2/(2*deta))*funspline(nn,1:nn-1);

funspline_x(2:nn-1,1:nn-1)=umat1\(lmat1\sm1);
funspline_x(1,1:nn-1)=(1/alfasp)*((betasp/deta)*(funspline(2,1:nn-1)-funspline(1,1:nn-1))... % U_x,0
    +(gamasp/(2*deta)*(funspline(3,1:nn-1)-funspline(1,1:nn-1)))...
    -betasp*funspline_x(2,1:nn-1)-gamasp*funspline_x(3,1:nn-1));
funspline_x(nn,1:nn-1)=(1/alfasp)*(-(betasp/deta)*(funspline(nn-1,1:nn-1)-funspline(nn,1:nn-1))... % U_x,N
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

for i=1:nn-1, 
    jj=jeta_c(i,nn+1-[1:nn]);
    xeta1=eta_c(i,nn+1-[1:nn]);
    xeta2=xeta1-eta(jj);
    xeta2=xeta2';
    funbIV1=((xc3(jj,i).*xeta2+xc2(jj,i)).*xeta2+xc1(jj,i)).*xeta2+xc0(jj,i);

    va_fI(3*nn-3+i,1:nn)=funbIV1(1:nn);
end

va_fI(:,1:nn)=ftr*va_fI(:,1:nn);

% AFFECTATIONS I+III
funftI(1:nn,1:nn)= va_fI(1:nn,1:nn);
funftIII(1:nn,nn-[1:nn]+1)= va_fI(2*nn-1:3*nn-2,1:nn);



%% RESEAU 3: 
va_fII=zeros(4*(nn-1),nn); 

% PANEL II :
va_fII(1:nn-1,1:nn)=funfII(1:nn-1,1:nn);

% PANEL III

fun=funfIII(1:nn-1,1:nn)';
funspline(1:nn,1:nn-1)=fun;
    %
    sm1=(3/dxi)*kmat*funspline(2:nn-1,1:nn-1);
    sm1(1,1:nn-1)=sm1(1,1:nn-1)+(m1/(2*dxi))*funspline(1,1:nn-1);
    sm1(n,1:nn-1)=sm1(n,1:nn-1)+(m2/(2*dxi))*funspline(nn,1:nn-1);


    funspline_x(2:nn-1,1:nn-1)=umat1\(lmat1\sm1);

    funspline_x(1,1:nn-1)=(1/alfasp)*((betasp/dxi)*(funspline(2,1:nn-1)-funspline(1,1:nn-1))... % U_x,0
        +(gamasp/(2*dxi)*(funspline(3,1:nn-1)-funspline(1,1:nn-1)))...
        -betasp*funspline_x(2,1:nn-1)-gamasp*funspline_x(3,1:nn-1));
    funspline_x(nn,1:nn-1)=(1/alfasp)*(-(betasp/dxi)*(funspline(nn-1,1:nn-1)-funspline(nn,1:nn-1))... % U_x,N
        -(gamasp/(2*dxi))*(funspline(nn-2,1:nn-1)-funspline(nn,1:nn-1))...
        +betasp*funspline_x(nn-1,1:nn-1)+gamasp*funspline_x(nn-2,1:nn-1));

     % COEFFTS C_0,C_1, C_2 AND C_3 OF CUBIC SPLINE J+1/2
    xc0(1:nn,1:nn-1)=funspline(1:nn,1:nn-1);
    xc1(1:nn,1:nn-1)=funspline_x(1:nn,1:nn-1);
    xc2(1:nn-1,1:nn-1)=3*(funspline(2:nn,1:nn-1)-funspline(1:nn-1,1:nn-1))/dxi^2 ...
                -(2*funspline_x(1:nn-1,1:nn-1)+funspline_x(2:nn,1:nn-1))/dxi;
    xc2(nn,1:nn-1)=0;        
    xc3(1:nn-1,1:nn-1)=2*(funspline(1:nn-1,1:nn-1)-funspline(2:nn,1:nn-1))/dxi^3 ...
                + (funspline_x(1:nn-1,1:nn-1)+funspline_x(2:nn,1:nn-1))/dxi^2;
    xc3(nn,1:nn-1)=0;

for i=1:nn-1
    jj=jeta_c(i,1:nn);
    xeta1=eta_c(i,1:nn);
    xeta2=xeta1-eta(jj);
    xeta2=xeta2';
    funbIII1=((xc3(jj,i).*xeta2+xc2(jj,i)).*xeta2+xc1(jj,i)).*xeta2+xc0(jj,i);

    va_fII(nn-1+i,1:nn)=funbIII1(1:nn);
end

% FACE IV: 
va_fII(2*nn-1:3*nn-3,1:nn)=funfIV(1:nn-1,nn:-1:1);   

% FACE I:

fun=funfI(1:nn-1,1:nn)';
funspline(1:nn,1:nn-1)=fun;

sm1=(3/deta)*kmat*funspline(2:nn-1,1:nn-1);
sm1(1,1:nn-1)=sm1(1,1:nn-1)+(m1/(2*deta))*funspline(1,1:nn-1);
sm1(n,1:nn-1)=sm1(n,1:nn-1)+(m2/(2*deta))*funspline(nn,1:nn-1);


funspline_x(2:nn-1,1:nn-1)=umat1\(lmat1\sm1);

funspline_x(1,1:nn-1)=(1/alfasp)*((betasp/deta)*(funspline(2,1:nn-1)-funspline(1,1:nn-1))... % U_x,0
    +(gamasp/(2*deta)*(funspline(3,1:nn-1)-funspline(1,1:nn-1)))...
    -betasp*funspline_x(2,1:nn-1)-gamasp*funspline_x(3,1:nn-1));
funspline_x(nn,1:nn-1)=(1/alfasp)*(-(betasp/deta)*(funspline(nn-1,1:nn-1)-funspline(nn,1:nn-1))... % U_x,N
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

for i=1:nn-1
    jj=jeta_c(i,nn+1-[1:nn]);
    xeta1=eta_c(i,nn+1-[1:nn]);
    xeta2=xeta1-eta(jj);
    xeta2=xeta2';
    funbI1=((xc3(jj,i).*xeta2+xc2(jj,i)).*xeta2+xc1(jj,i)).*xeta2+xc0(jj,i);

    va_fII(3*nn-3+i,1:nn)=funbI1(1:nn);
end

va_fII(:,1:nn)=ftr*va_fII;
% TRANSFERT FACE II+IV
funftII(1:nn,1:nn) = va_fII(1:nn,1:nn);
funftIV(1:nn,nn-[1:nn]+1) = va_fII(2*nn-1:3*nn-2,1:nn);

%% RESEAU 5: 
va_fV=zeros(4*(nn-1),nn);
funaII1=zeros(nn,1);
funaIV1=zeros(nn,1);

% PANEL V
va_fV(1:nn-1,1:nn)=funfV(1:nn-1,1:nn);

% PANEL II: 
fun=funfII(1:nn,nn:-1:2);
funspline(1:nn,1:nn-1)=fun;
    
sm1=(3/dxi)*kmat*funspline(2:nn-1,1:nn-1);
sm1(1,1:nn-1)=sm1(1,1:nn-1)+(m1/(2*dxi))*funspline(1,1:nn-1);
sm1(n,1:nn-1)=sm1(n,1:nn-1)+(m2/(2*dxi))*funspline(nn,1:nn-1);

funspline_x(2:nn-1,1:nn-1)=umat1\(lmat1\sm1);

funspline_x(1,1:nn-1)=(1/alfasp)*((betasp/dxi)*(funspline(2,1:nn-1)-funspline(1,1:nn-1))... % U_x,0
    +(gamasp/(2*dxi)*(funspline(3,1:nn-1)-funspline(1,1:nn-1)))...
    -betasp*funspline_x(2,1:nn-1)-gamasp*funspline_x(3,1:nn-1));
funspline_x(nn,1:nn-1)=(1/alfasp)*(-(betasp/dxi)*(funspline(nn-1,1:nn-1)-funspline(nn,1:nn-1))... % U_x,N
    -(gamasp/(2*dxi))*(funspline(nn-2,1:nn-1)-funspline(nn,1:nn-1))...
    +betasp*funspline_x(nn-1,1:nn-1)+gamasp*funspline_x(nn-2,1:nn-1));

 % COEFFTS C_0,C_1, C_2 AND C_3 OF CUBIC SPLINE J+1/2
xc0(1:nn,1:nn-1)=funspline(1:nn,1:nn-1);
xc1(1:nn,1:nn-1)=funspline_x(1:nn,1:nn-1);
xc2(1:nn-1,1:nn-1)=3*(funspline(2:nn,1:nn-1)-funspline(1:nn-1,1:nn-1))/dxi^2 ...
            -(2*funspline_x(1:nn-1,1:nn-1)+funspline_x(2:nn,1:nn-1))/dxi;
xc2(nn,1:nn-1)=0;        
xc3(1:nn-1,1:nn-1)=2*(funspline(1:nn-1,1:nn-1)-funspline(2:nn,1:nn-1))/dxi^3 ...
            + (funspline_x(1:nn-1,1:nn-1)+funspline_x(2:nn,1:nn-1))/dxi^2;
xc3(nn,1:nn-1)=0; 
  
for j=1:nn-1
    ii=ixi_c(1:nn,j);
    xxi1=xi_c(1:nn,j);
    xxi2=xxi1'-xi(ii);
    xxi2=xxi2';
    funaII1(1:nn)=((xc3(ii,j).*xxi2+xc2(ii,j)).*xxi2+xc1(ii,j)).*xxi2+xc0(ii,j);

    va_fV(nn-1+j,1:nn)=funaII1(1:nn);
end

% PANEL VI: 
va_fV(2*nn-1:3*nn-3,1:nn)=funfVI(nn:-1:2,1:nn);

% PANEL IV:
fun=funfIV(1:nn,1:nn-1);
funspline(1:nn,1:nn-1)=fun;
    
sm1=(3/dxi)*kmat*funspline(2:nn-1,1:nn-1);
sm1(1,1:nn-1)=sm1(1,1:nn-1)+(m1/(2*dxi))*funspline(1,1:nn-1);
sm1(n,1:nn-1)=sm1(n,1:nn-1)+(m2/(2*dxi))*funspline(nn,1:nn-1);

funspline_x(2:nn-1,1:nn-1)=umat1\(lmat1\sm1);

funspline_x(1,1:nn-1)=(1/alfasp)*((betasp/dxi)*(funspline(2,1:nn-1)-funspline(1,1:nn-1))... % U_x,0
    +(gamasp/(2*dxi)*(funspline(3,1:nn-1)-funspline(1,1:nn-1)))...
    -betasp*funspline_x(2,1:nn-1)-gamasp*funspline_x(3,1:nn-1));
funspline_x(nn,1:nn-1)=(1/alfasp)*(-(betasp/dxi)*(funspline(nn-1,1:nn-1)-funspline(nn,1:nn-1))... % U_x,N
    -(gamasp/(2*dxi))*(funspline(nn-2,1:nn-1)-funspline(nn,1:nn-1))...
    +betasp*funspline_x(nn-1,1:nn-1)+gamasp*funspline_x(nn-2,1:nn-1));

 % COEFFTS C_0,C_1, C_2 AND C_3 OF CUBIC SPLINE J+1/2
xc0(1:nn,1:nn-1)=funspline(1:nn,1:nn-1);
xc1(1:nn,1:nn-1)=funspline_x(1:nn,1:nn-1);
xc2(1:nn-1,1:nn-1)=3*(funspline(2:nn,1:nn-1)-funspline(1:nn-1,1:nn-1))/dxi^2 ...
            -(2*funspline_x(1:nn-1,1:nn-1)+funspline_x(2:nn,1:nn-1))/dxi;
xc2(nn,1:nn-1)=0;        
xc3(1:nn-1,1:nn-1)=2*(funspline(1:nn-1,1:nn-1)-funspline(2:nn,1:nn-1))/dxi^3 ...
            + (funspline_x(1:nn-1,1:nn-1)+funspline_x(2:nn,1:nn-1))/dxi^2;
xc3(nn,1:nn-1)=0; 

for j=1:nn-1
    ii=ixi_c(1:nn,j);
    xxi1=xi_c(1:nn,j);
    xxi2=xxi1'-xi(ii);
    xxi2=xxi2';
    funaIV1(1:nn)=((xc3(ii,j).*xxi2+xc2(ii,j)).*xxi2+xc1(ii,j)).*xxi2+xc0(ii,j);
    
    va_fV(3*nn-3+j,1:nn)=funaIV1(1:nn);
end

va_fV=ftr*va_fV;

% TRANSFERT FACE V+VI
funftV(1:nn,1:nn) = va_fV(1:nn,1:nn);
funftVI(nn:-1:1,1:nn) = va_fV(2*nn-1:3*nn-2,1:nn);


%% RESEAU 6: ASSEMBLAGE DES DONNEES SUR LE RESEAU V-BETA
vb_fV=zeros(nn,4*(nn-1));
funaIII1=zeros(nn,1);
funaI1=zeros(nn,1);

% PANEL V :
vb_fV(1:nn,1:nn-1)=funftV(1:nn,1:nn-1);

% PANEL III
fun=funfIII(1:nn,nn+1-[1:nn-1]);
funspline(1:nn,1:nn-1)=fun;

sm1=(3/dxi)*kmat*funspline(2:nn-1,1:nn-1);
sm1(1,1:nn-1)=sm1(1,1:nn-1)+(m1/(2*dxi))*funspline(1,1:nn-1);
sm1(n,1:nn-1)=sm1(n,1:nn-1)+(m2/(2*dxi))*funspline(nn,1:nn-1);  

funspline_x(2:nn-1,1:nn-1)=umat1\(lmat1\sm1);

funspline_x(1,1:nn-1)=(1/alfasp)*((betasp/dxi)*(funspline(2,1:nn-1)-funspline(1,1:nn-1))... % U_x,0
    +(gamasp/(2*dxi)*(funspline(3,1:nn-1)-funspline(1,1:nn-1)))...
    -betasp*funspline_x(2,1:nn-1)-gamasp*funspline_x(3,1:nn-1));
funspline_x(nn,1:nn-1)=(1/alfasp)*(-(betasp/dxi)*(funspline(nn-1,1:nn-1)-funspline(nn,1:nn-1))... % U_x,N
    -(gamasp/(2*dxi))*(funspline(nn-2,1:nn-1)-funspline(nn,1:nn-1))...
    +betasp*funspline_x(nn-1,1:nn-1)+gamasp*funspline_x(nn-2,1:nn-1));

% COEFFTS C_0,C_1, C_2 AND C_3 OF CUBIC SPLINE J+1/2
xc0(1:nn,1:nn-1)=funspline(1:nn,1:nn-1);
xc1(1:nn,1:nn-1)=funspline_x(1:nn,1:nn-1);
xc2(1:nn-1,1:nn-1)=3*(funspline(2:nn,1:nn-1)-funspline(1:nn-1,1:nn-1))/dxi^2 ...
            -(2*funspline_x(1:nn-1,1:nn-1)+funspline_x(2:nn,1:nn-1))/dxi;
xc2(nn,1:nn-1)=0;        
xc3(1:nn-1,1:nn-1)=2*(funspline(1:nn-1,1:nn-1)-funspline(2:nn,1:nn-1))/dxi^3 ...
            + (funspline_x(1:nn-1,1:nn-1)+funspline_x(2:nn,1:nn-1))/dxi^2;
xc3(nn,1:nn-1)=0; 

for j=1:nn-1
    ii=ixi_c(nn+1-[1:nn],j);
    xxi1=xi_c(nn+1-[1:nn],j);
    xxi2=xxi1'-xi(ii);
    xxi2=xxi2';
    funaIII1(1:nn)=((xc3(ii,j).*xxi2+xc2(ii,j)).*xxi2+xc1(ii,j)).*xxi2+xc0(ii,j);
    
    vb_fV(1:nn,nn-1+j)=funaIII1(1:nn);
end

% PANEL VI :
vb_fV(1:nn,2*nn-1:3*nn-3)=funftVI(nn:-1:1,1:nn-1); 

% PANEL I
fun=funfI(1:nn,1:nn-1);
funspline(1:nn,1:nn-1)=fun;

sm1=(3/dxi)*kmat*funspline(2:nn-1,1:nn-1);
sm1(1,1:nn-1)=sm1(1,1:nn-1)+(m1/(2*dxi))*funspline(1,1:nn-1);
sm1(n,1:nn-1)=sm1(n,1:nn-1)+(m2/(2*dxi))*funspline(nn,1:nn-1);  

funspline_x(2:nn-1,1:nn-1)=umat1\(lmat1\sm1);

funspline_x(1,1:nn-1)=(1/alfasp)*((betasp/dxi)*(funspline(2,1:nn-1)-funspline(1,1:nn-1))... % U_x,0
    +(gamasp/(2*dxi)*(funspline(3,1:nn-1)-funspline(1,1:nn-1)))...
    -betasp*funspline_x(2,1:nn-1)-gamasp*funspline_x(3,1:nn-1));
funspline_x(nn,1:nn-1)=(1/alfasp)*(-(betasp/dxi)*(funspline(nn-1,1:nn-1)-funspline(nn,1:nn-1))... % U_x,N
    -(gamasp/(2*dxi))*(funspline(nn-2,1:nn-1)-funspline(nn,1:nn-1))...
    +betasp*funspline_x(nn-1,1:nn-1)+gamasp*funspline_x(nn-2,1:nn-1));

% COEFFTS C_0,C_1, C_2 AND C_3 OF CUBIC SPLINE J+1/2
xc0(1:nn,1:nn-1)=funspline(1:nn,1:nn-1);
xc1(1:nn,1:nn-1)=funspline_x(1:nn,1:nn-1);
xc2(1:nn-1,1:nn-1)=3*(funspline(2:nn,1:nn-1)-funspline(1:nn-1,1:nn-1))/dxi^2 ...
            -(2*funspline_x(1:nn-1,1:nn-1)+funspline_x(2:nn,1:nn-1))/dxi;
xc2(nn,1:nn-1)=0;        
xc3(1:nn-1,1:nn-1)=2*(funspline(1:nn-1,1:nn-1)-funspline(2:nn,1:nn-1))/dxi^3 ...
            + (funspline_x(1:nn-1,1:nn-1)+funspline_x(2:nn,1:nn-1))/dxi^2;
xc3(nn,1:nn-1)=0; 

for j=1:nn-1
    ii=ixi_c(nn+1-[1:nn],j);
    xxi1=xi_c(nn+1-[1:nn],j);
    xxi2=xxi1'-xi(ii);
    xxi2=xxi2';
    funaI1(1:nn)=((xc3(ii,j).*xxi2+xc2(ii,j)).*xxi2+xc1(ii,j)).*xxi2+xc0(ii,j);
    
    vb_fV(1:nn,3*nn-3+j)=funaI1(1:nn);
end

vb_fV=(ftr*vb_fV')';

% TRANSFERT FACE V+VI

funftV(1:nn,1:nn) = vb_fV(1:nn,1:nn);
funftVI(nn-[1:nn]+1,1:nn) = vb_fV(1:nn,2*nn-1:3*nn-2);
