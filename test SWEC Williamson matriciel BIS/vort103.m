function [vort_I,vort_II,vort_III,vort_IV,vort_V,vort_VI]=...
    vort101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn)
global xi eta dxi deta
global gxi_I gxi_II gxi_III gxi_IV gxi_V gxi_VI;
global geta_I geta_II geta_III geta_IV geta_V geta_VI;
global gr_I gr_II gr_III gr_IV gr_V gr_VI
global pg kg;
global kmat umat1 lmat1
global eta_c jeta_c xi_c ixi_c
global alfasp betasp gamasp
global m1 m2

%% --- COORDINATE 1 -------------------------------------------------------

%% RESEAU 1 : 
% PANEL I :
va_fI(1:nn-1,1:nn)=funfI(1:nn-1,1:nn,1);

% PANEL II :
fun=funfII(1:nn-1,1:nn,1)';
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
va_fI(2*nn-1:3*nn-3,1:nn)=funfIII(1:nn-1,nn:-1:1,1);

%  PANEL IV:
fun=funfIV(1:nn-1,1:nn,1)';
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

funaI=va_fI(:,1:nn);
xsmIa=kg*funaI;
vad_fI(:,1:nn)=pg\xsmIa;

%% RESEAU 2

% PANEL I :
vb_fI(1:nn,1:nn-1)=funfI(1:nn,1:nn-1,1);

% FACE V:
fun=funfV(1:nn,1:nn-1,1);
funspline(1:nn,1:nn-1)=fun;
sm1=(3/dxi)*kmat*funspline(2:nn-1,1:nn-1);
sm1(1,1:nn-1)=sm1(1,1:nn-1)+(m1/(2*dxi))*funspline(1,1:nn-1);
sm1(n,1:nn-1)=sm1(n,1:nn-1)+(m2/(2*dxi))*funspline(nn,1:nn-1);    

funspline_x(2:nn-1,1:nn-1)=umat1\(lmat1\sm1);

funspline_x(1,1:nn-1)=(1/alfasp)*((betasp/dxi)*(funspline(2,1:nn-1)-funspline(1,1:nn-1))... % U_x,0
    +(gamasp/(2*dxi)*(funspline(3,1:nn-1)-funspline(1,1:nn-1)))...
    - betasp*funspline_x(2,1:nn-1)-gamasp*funspline_x(3,1:nn-1));
funspline_x(nn,1:nn-1)=(1/alfasp)*(-(betasp/dxi)*(funspline(nn-1,1:nn-1)-funspline(nn,1:nn-1))... % U_x,N
    -(gamasp/(2*dxi))*(funspline(nn-2,1:nn-1)-funspline(nn,1:nn-1))...
    + betasp*funspline_x(nn-1,1:nn-1)+gamasp*funspline_x(nn-2,1:nn-1));

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
    funaV1=((xc3(ii,j).*xxi2+xc2(ii,j)).*xxi2+xc1(ii,j)).*xxi2+xc0(ii,j);

    vb_fI(1:nn,nn-1+j)=funaV1(1:nn);
end

% FACE III:
vb_fI(1:nn,2*nn-1:3*nn-3)=funfIII(1:nn,nn:-1:2,1); 
 
% FACE VI: 
fun=funfVI(1:nn,1:nn-1,1); 
funspline(1:nn,1:nn-1)=fun;
sm1=(3/dxi)*kmat*funspline(2:nn-1,1:nn-1);
sm1(1,1:nn-1)=sm1(1,1:nn-1)+(m1/(2*dxi))*funspline(1,1:nn-1);
sm1(n,1:nn-1)=sm1(n,1:nn-1)+(m2/(2*dxi))*funspline(nn,1:nn-1);    

funspline_x(2:nn-1,1:nn-1)=umat1\(lmat1\sm1);

funspline_x(1,1:nn-1)=(1/alfasp)*((betasp/dxi)*(funspline(2,1:nn-1)-funspline(1,1:nn-1))... % U_x,0
    +(gamasp/(2*dxi)*(funspline(3,1:nn-1)-funspline(1,1:nn-1)))...
    - betasp*funspline_x(2,1:nn-1)-gamasp*funspline_x(3,1:nn-1));
funspline_x(nn,1:nn-1)=(1/alfasp)*(-(betasp/dxi)*(funspline(nn-1,1:nn-1)-funspline(nn,1:nn-1))... % U_x,N
    -(gamasp/(2*dxi))*(funspline(nn-2,1:nn-1)-funspline(nn,1:nn-1))...
    + betasp*funspline_x(nn-1,1:nn-1)+gamasp*funspline_x(nn-2,1:nn-1));

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
    funaVI1=((xc3(ii,j).*xxi2+xc2(ii,j)).*xxi2+xc1(ii,j)).*xxi2+xc0(ii,j);

    vb_fI(1:nn,3*nn-3+j)=funaVI1(1:nn);
end

% CALCUL DES DERIVEES BETA SUR LE RESEAU DE GRANDS CERCLES I-BETA
funbI=vb_fI(1:nn,:)';
xsmIb=kg*funbI;
vbd_fI(1:nn,:)=(pg\xsmIb)';

%% RESEAU 3: 

% PANEL II :
va_fII(1:nn-1,1:nn)=funfII(1:nn-1,1:nn,1);

% PANEL III

fun=funfIII(1:nn-1,1:nn,1)';
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
va_fII(2*nn-1:3*nn-3,1:nn)=funfIV(1:nn-1,nn:-1:1,1);   

% FACE I:

fun=funfI(1:nn-1,1:nn,1)';
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

funaII=va_fII(:,1:nn);
xsmIIa=kg*funaII;
vad_fII(:,1:nn)=pg\xsmIIa;

%% SET 4
vb_fII=zeros(nn,4*(nn-1)); 

% PANEL II : 
vb_fII(1:nn,1:nn-1)=funfII(1:nn,1:nn-1,1);

% PANEL V
fun=funfV(nn:-1:2,1:nn,1)';
funspline(1:nn,1:nn-1)=fun;

sm1=(3/dxi)*kmat*funspline(2:nn-1,1:nn-1);
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
    funbV1=((xc3(jj,i).*xeta2+xc2(jj,i)).*xeta2+xc1(jj,i)).*xeta2+xc0(jj,i);
     vb_fII(1:nn,nn-1+i)=funbV1(1:nn);
end

% PANEL III:
vb_fII(1:nn,2*nn-1:3*nn-3)=funfIV(1:nn,nn:-1:2,1); 

% PANEL VI: 
fun=funfVI(1:nn-1,1:nn,1)';
funspline(1:nn,1:nn-1)=fun;

sm1=(3/dxi)*kmat*funspline(2:nn-1,1:nn-1);
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
    funbVI1=((xc3(jj,i).*xeta2+xc2(jj,i)).*xeta2+xc1(jj,i)).*xeta2+xc0(jj,i);
    %
    vb_fII(1:nn,3*nn-3+i)=funbVI1(1:nn);
end

funbII=vb_fII(1:nn,:)';
xsmIIb=kg*funbII;
vbd_fII(1:nn,:)=(pg\xsmIIb)';

%% SET 5
% PANEL V
va_fV(1:nn-1,1:nn)=funfV(1:nn-1,1:nn,1);

% PANEL II: 
fun=funfII(1:nn,nn:-1:2,1);
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
va_fV(2*nn-1:3*nn-3,1:nn)=funfVI(nn:-1:2,1:nn,1);

% PANEL IV:
fun=funfIV(1:nn,1:nn-1,1);
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

funaV=va_fV(:,1:nn);
xsmVa=kg*funaV;
vad_fV(:,1:nn)=pg\xsmVa;

%% SET 6 :

% PANEL V :
vb_fV(1:nn,1:nn-1)=funfV(1:nn,1:nn-1,1);

% PANEL III
fun=funfIII(1:nn,nn+1-[1:nn-1],1);
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
vb_fV(1:nn,2*nn-1:3*nn-3)=funfVI(nn:-1:1,1:nn-1,1); 

% PANEL I
fun=funfI(1:nn,1:nn-1,1);
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


funbV=vb_fV(1:nn,:)';
xsmVb=kg*funbV;
vbd_fV(1:nn,:)=(pg\xsmVb)';

%% AFFECTATIONS
dg_xi_fI(1:nn,1:nn,1)=vad_fI([1:nn],[1:nn]);
dg_eta_fI(1:nn,1:nn,1)=vbd_fI([1:nn],[1:nn]);

dg_xi_fII(1:nn,1:nn,1)=vad_fII([1:nn],[1:nn]);
dg_eta_fII(1:nn,1:nn,1)=vbd_fII([1:nn],[1:nn]);

dg_xi_fIII(1:nn,1:nn,1)=vad_fI(2*(nn-1)+[1:nn],nn+1-[1:nn]);
dg_eta_fIII(1:nn,1:nn,1)=-vbd_fI([1:nn],2*(nn-1)+nn-[1:nn]+1); 

dg_xi_fIV(1:nn,1:nn,1)=vad_fII(2*(nn-1)+[1:nn],nn-[1:nn]+1);
dg_eta_fIV(1:nn,1:nn,1)=-vbd_fII([1:nn],2*(nn-1)+nn-[1:nn]+1);

dg_xi_fV(1:nn,1:nn,1)=vad_fV([1:nn],[1:nn]);
dg_eta_fV(1:nn,1:nn,1)=vbd_fV([1:nn],[1:nn]);

dg_xi_fVI(1:nn,1:nn,1)=-vad_fV(2*(nn-1)+nn+1-[1:nn],[1:nn]);
dg_eta_fVI(1:nn,1:nn,1)=vbd_fV(nn-[1:nn]+1,2*(nn-1)+[1:nn]);

%% --- COORDINATE 1 -------------------------------------------------------

%% RESEAU 1 : 
% PANEL I :
va_fI(1:nn-1,1:nn)=funfI(1:nn-1,1:nn,2);

% PANEL II :
fun=funfII(1:nn-1,1:nn,2)';
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
va_fI(2*nn-1:3*nn-3,1:nn)=funfIII(1:nn-1,nn:-1:1,2);

%  PANEL IV:
fun=funfIV(1:nn-1,1:nn,2)';
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

funaI=va_fI(:,1:nn);
xsmIa=kg*funaI;
vad_fI(:,1:nn)=pg\xsmIa;

%% RESEAU 2

% PANEL I :
vb_fI(1:nn,1:nn-1)=funfI(1:nn,1:nn-1,2);

% FACE V:
fun=funfV(1:nn,1:nn-1,2);
funspline(1:nn,1:nn-1)=fun;
sm1=(3/dxi)*kmat*funspline(2:nn-1,1:nn-1);
sm1(1,1:nn-1)=sm1(1,1:nn-1)+(m1/(2*dxi))*funspline(1,1:nn-1);
sm1(n,1:nn-1)=sm1(n,1:nn-1)+(m2/(2*dxi))*funspline(nn,1:nn-1);    

funspline_x(2:nn-1,1:nn-1)=umat1\(lmat1\sm1);

funspline_x(1,1:nn-1)=(1/alfasp)*((betasp/dxi)*(funspline(2,1:nn-1)-funspline(1,1:nn-1))... % U_x,0
    +(gamasp/(2*dxi)*(funspline(3,1:nn-1)-funspline(1,1:nn-1)))...
    - betasp*funspline_x(2,1:nn-1)-gamasp*funspline_x(3,1:nn-1));
funspline_x(nn,1:nn-1)=(1/alfasp)*(-(betasp/dxi)*(funspline(nn-1,1:nn-1)-funspline(nn,1:nn-1))... % U_x,N
    -(gamasp/(2*dxi))*(funspline(nn-2,1:nn-1)-funspline(nn,1:nn-1))...
    + betasp*funspline_x(nn-1,1:nn-1)+gamasp*funspline_x(nn-2,1:nn-1));

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
    funaV1=((xc3(ii,j).*xxi2+xc2(ii,j)).*xxi2+xc1(ii,j)).*xxi2+xc0(ii,j);

    vb_fI(1:nn,nn-1+j)=funaV1(1:nn);
end

% FACE III:
vb_fI(1:nn,2*nn-1:3*nn-3)=funfIII(1:nn,nn:-1:2,2); 
 
% FACE VI: 
fun=funfVI(1:nn,1:nn-1,2); 
funspline(1:nn,1:nn-1)=fun;
sm1=(3/dxi)*kmat*funspline(2:nn-1,1:nn-1);
sm1(1,1:nn-1)=sm1(1,1:nn-1)+(m1/(2*dxi))*funspline(1,1:nn-1);
sm1(n,1:nn-1)=sm1(n,1:nn-1)+(m2/(2*dxi))*funspline(nn,1:nn-1);    

funspline_x(2:nn-1,1:nn-1)=umat1\(lmat1\sm1);

funspline_x(1,1:nn-1)=(1/alfasp)*((betasp/dxi)*(funspline(2,1:nn-1)-funspline(1,1:nn-1))... % U_x,0
    +(gamasp/(2*dxi)*(funspline(3,1:nn-1)-funspline(1,1:nn-1)))...
    - betasp*funspline_x(2,1:nn-1)-gamasp*funspline_x(3,1:nn-1));
funspline_x(nn,1:nn-1)=(1/alfasp)*(-(betasp/dxi)*(funspline(nn-1,1:nn-1)-funspline(nn,1:nn-1))... % U_x,N
    -(gamasp/(2*dxi))*(funspline(nn-2,1:nn-1)-funspline(nn,1:nn-1))...
    + betasp*funspline_x(nn-1,1:nn-1)+gamasp*funspline_x(nn-2,1:nn-1));

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
    funaVI1=((xc3(ii,j).*xxi2+xc2(ii,j)).*xxi2+xc1(ii,j)).*xxi2+xc0(ii,j);

    vb_fI(1:nn,3*nn-3+j)=funaVI1(1:nn);
end

% CALCUL DES DERIVEES BETA SUR LE RESEAU DE GRANDS CERCLES I-BETA
funbI=vb_fI(1:nn,:)';
xsmIb=kg*funbI;
vbd_fI(1:nn,:)=(pg\xsmIb)';

%% RESEAU 3: 

% PANEL II :
va_fII(1:nn-1,1:nn)=funfII(1:nn-1,1:nn,2);

% PANEL III

fun=funfIII(1:nn-1,1:nn,2)';
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
va_fII(2*nn-1:3*nn-3,1:nn)=funfIV(1:nn-1,nn:-1:1,2);   

% FACE I:

fun=funfI(1:nn-1,1:nn,2)';
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

funaII=va_fII(:,1:nn);
xsmIIa=kg*funaII;
vad_fII(:,1:nn)=pg\xsmIIa;

%% SET 4
vb_fII=zeros(nn,4*(nn-1)); 

% PANEL II : 
vb_fII(1:nn,1:nn-1)=funfII(1:nn,1:nn-1,2);

% PANEL V
fun=funfV(nn:-1:2,1:nn,2)';
funspline(1:nn,1:nn-1)=fun;

sm1=(3/dxi)*kmat*funspline(2:nn-1,1:nn-1);
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
    funbV1=((xc3(jj,i).*xeta2+xc2(jj,i)).*xeta2+xc1(jj,i)).*xeta2+xc0(jj,i);
     vb_fII(1:nn,nn-1+i)=funbV1(1:nn);
end

% PANEL III:
vb_fII(1:nn,2*nn-1:3*nn-3)=funfIV(1:nn,nn:-1:2,2); 

% PANEL VI: 
fun=funfVI(1:nn-1,1:nn,2)';
funspline(1:nn,1:nn-1)=fun;

sm1=(3/dxi)*kmat*funspline(2:nn-1,1:nn-1);
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
    funbVI1=((xc3(jj,i).*xeta2+xc2(jj,i)).*xeta2+xc1(jj,i)).*xeta2+xc0(jj,i);
    %
    vb_fII(1:nn,3*nn-3+i)=funbVI1(1:nn);
end

funbII=vb_fII(1:nn,:)';
xsmIIb=kg*funbII;
vbd_fII(1:nn,:)=(pg\xsmIIb)';

%% SET 5
% PANEL V
va_fV(1:nn-1,1:nn)=funfV(1:nn-1,1:nn,2);

% PANEL II: 
fun=funfII(1:nn,nn:-1:2,2);
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
va_fV(2*nn-1:3*nn-3,1:nn)=funfVI(nn:-1:2,1:nn,2);

% PANEL IV:
fun=funfIV(1:nn,1:nn-1,2);
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

funaV=va_fV(:,1:nn);
xsmVa=kg*funaV;
vad_fV(:,1:nn)=pg\xsmVa;

%% SET 6 :

% PANEL V :
vb_fV(1:nn,1:nn-1)=funfV(1:nn,1:nn-1,2);

% PANEL III
fun=funfIII(1:nn,nn+1-[1:nn-1],2);
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
vb_fV(1:nn,2*nn-1:3*nn-3)=funfVI(nn:-1:1,1:nn-1,2); 

% PANEL I
fun=funfI(1:nn,1:nn-1,2);
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


funbV=vb_fV(1:nn,:)';
xsmVb=kg*funbV;
vbd_fV(1:nn,:)=(pg\xsmVb)';

%% AFFECTATIONS
dg_xi_fI(1:nn,1:nn,2)=vad_fI([1:nn],[1:nn]);
dg_eta_fI(1:nn,1:nn,2)=vbd_fI([1:nn],[1:nn]);

dg_xi_fII(1:nn,1:nn,2)=vad_fII([1:nn],[1:nn]);
dg_eta_fII(1:nn,1:nn,2)=vbd_fII([1:nn],[1:nn]);

dg_xi_fIII(1:nn,1:nn,2)=vad_fI(2*(nn-1)+[1:nn],nn+1-[1:nn]);
dg_eta_fIII(1:nn,1:nn,2)=-vbd_fI([1:nn],2*(nn-1)+nn-[1:nn]+1); 

dg_xi_fIV(1:nn,1:nn,2)=vad_fII(2*(nn-1)+[1:nn],nn-[1:nn]+1);
dg_eta_fIV(1:nn,1:nn,2)=-vbd_fII([1:nn],2*(nn-1)+nn-[1:nn]+1);

dg_xi_fV(1:nn,1:nn,2)=vad_fV([1:nn],[1:nn]);
dg_eta_fV(1:nn,1:nn,2)=vbd_fV([1:nn],[1:nn]);

dg_xi_fVI(1:nn,1:nn,2)=-vad_fV(2*(nn-1)+nn+1-[1:nn],[1:nn]);
dg_eta_fVI(1:nn,1:nn,2)=vbd_fV(nn-[1:nn]+1,2*(nn-1)+[1:nn]);

%% --- COORDINATE 3 -------------------------------------------------------

%% RESEAU 1 : 
% PANEL I :
va_fI(1:nn-1,1:nn)=funfI(1:nn-1,1:nn,3);

% PANEL II :
fun=funfII(1:nn-1,1:nn,3)';
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
va_fI(2*nn-1:3*nn-3,1:nn)=funfIII(1:nn-1,nn:-1:1,3);

%  PANEL IV:
fun=funfIV(1:nn-1,1:nn,3)';
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

funaI=va_fI(:,1:nn);
xsmIa=kg*funaI;
vad_fI(:,1:nn)=pg\xsmIa;

%% RESEAU 2

% PANEL I :
vb_fI(1:nn,1:nn-1)=funfI(1:nn,1:nn-1,3);

% FACE V:
fun=funfV(1:nn,1:nn-1,3);
funspline(1:nn,1:nn-1)=fun;
sm1=(3/dxi)*kmat*funspline(2:nn-1,1:nn-1);
sm1(1,1:nn-1)=sm1(1,1:nn-1)+(m1/(2*dxi))*funspline(1,1:nn-1);
sm1(n,1:nn-1)=sm1(n,1:nn-1)+(m2/(2*dxi))*funspline(nn,1:nn-1);    

funspline_x(2:nn-1,1:nn-1)=umat1\(lmat1\sm1);

funspline_x(1,1:nn-1)=(1/alfasp)*((betasp/dxi)*(funspline(2,1:nn-1)-funspline(1,1:nn-1))... % U_x,0
    +(gamasp/(2*dxi)*(funspline(3,1:nn-1)-funspline(1,1:nn-1)))...
    - betasp*funspline_x(2,1:nn-1)-gamasp*funspline_x(3,1:nn-1));
funspline_x(nn,1:nn-1)=(1/alfasp)*(-(betasp/dxi)*(funspline(nn-1,1:nn-1)-funspline(nn,1:nn-1))... % U_x,N
    -(gamasp/(2*dxi))*(funspline(nn-2,1:nn-1)-funspline(nn,1:nn-1))...
    + betasp*funspline_x(nn-1,1:nn-1)+gamasp*funspline_x(nn-2,1:nn-1));

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
    funaV1=((xc3(ii,j).*xxi2+xc2(ii,j)).*xxi2+xc1(ii,j)).*xxi2+xc0(ii,j);

    vb_fI(1:nn,nn-1+j)=funaV1(1:nn);
end

% FACE III:
vb_fI(1:nn,2*nn-1:3*nn-3)=funfIII(1:nn,nn:-1:2,3); 
 
% FACE VI: 
fun=funfVI(1:nn,1:nn-1,3); 
funspline(1:nn,1:nn-1)=fun;
sm1=(3/dxi)*kmat*funspline(2:nn-1,1:nn-1);
sm1(1,1:nn-1)=sm1(1,1:nn-1)+(m1/(2*dxi))*funspline(1,1:nn-1);
sm1(n,1:nn-1)=sm1(n,1:nn-1)+(m2/(2*dxi))*funspline(nn,1:nn-1);    

funspline_x(2:nn-1,1:nn-1)=umat1\(lmat1\sm1);

funspline_x(1,1:nn-1)=(1/alfasp)*((betasp/dxi)*(funspline(2,1:nn-1)-funspline(1,1:nn-1))... % U_x,0
    +(gamasp/(2*dxi)*(funspline(3,1:nn-1)-funspline(1,1:nn-1)))...
    - betasp*funspline_x(2,1:nn-1)-gamasp*funspline_x(3,1:nn-1));
funspline_x(nn,1:nn-1)=(1/alfasp)*(-(betasp/dxi)*(funspline(nn-1,1:nn-1)-funspline(nn,1:nn-1))... % U_x,N
    -(gamasp/(2*dxi))*(funspline(nn-2,1:nn-1)-funspline(nn,1:nn-1))...
    + betasp*funspline_x(nn-1,1:nn-1)+gamasp*funspline_x(nn-2,1:nn-1));

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
    funaVI1=((xc3(ii,j).*xxi2+xc2(ii,j)).*xxi2+xc1(ii,j)).*xxi2+xc0(ii,j);

    vb_fI(1:nn,3*nn-3+j)=funaVI1(1:nn);
end

% CALCUL DES DERIVEES BETA SUR LE RESEAU DE GRANDS CERCLES I-BETA
funbI=vb_fI(1:nn,:)';
xsmIb=kg*funbI;
vbd_fI(1:nn,:)=(pg\xsmIb)';

%% RESEAU 3: 

% PANEL II :
va_fII(1:nn-1,1:nn)=funfII(1:nn-1,1:nn,3);

% PANEL III

fun=funfIII(1:nn-1,1:nn,3)';
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
va_fII(2*nn-1:3*nn-3,1:nn)=funfIV(1:nn-1,nn:-1:1,3);   

% FACE I:

fun=funfI(1:nn-1,1:nn,3)';
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

funaII=va_fII(:,1:nn);
xsmIIa=kg*funaII;
vad_fII(:,1:nn)=pg\xsmIIa;

%% SET 4
vb_fII=zeros(nn,4*(nn-1)); 

% PANEL II : 
vb_fII(1:nn,1:nn-1)=funfII(1:nn,1:nn-1,3);

% PANEL V
fun=funfV(nn:-1:2,1:nn,3)';
funspline(1:nn,1:nn-1)=fun;

sm1=(3/dxi)*kmat*funspline(2:nn-1,1:nn-1);
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
    funbV1=((xc3(jj,i).*xeta2+xc2(jj,i)).*xeta2+xc1(jj,i)).*xeta2+xc0(jj,i);
     vb_fII(1:nn,nn-1+i)=funbV1(1:nn);
end

% PANEL III:
vb_fII(1:nn,2*nn-1:3*nn-3)=funfIV(1:nn,nn:-1:2,3); 

% PANEL VI: 
fun=funfVI(1:nn-1,1:nn,3)';
funspline(1:nn,1:nn-1)=fun;

sm1=(3/dxi)*kmat*funspline(2:nn-1,1:nn-1);
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
    funbVI1=((xc3(jj,i).*xeta2+xc2(jj,i)).*xeta2+xc1(jj,i)).*xeta2+xc0(jj,i);
    %
    vb_fII(1:nn,3*nn-3+i)=funbVI1(1:nn);
end

funbII=vb_fII(1:nn,:)';
xsmIIb=kg*funbII;
vbd_fII(1:nn,:)=(pg\xsmIIb)';

%% SET 5
% PANEL V
va_fV(1:nn-1,1:nn)=funfV(1:nn-1,1:nn,3);

% PANEL II: 
fun=funfII(1:nn,nn:-1:2,3);
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
va_fV(2*nn-1:3*nn-3,1:nn)=funfVI(nn:-1:2,1:nn,3);

% PANEL IV:
fun=funfIV(1:nn,1:nn-1,3);
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

funaV=va_fV(:,1:nn);
xsmVa=kg*funaV;
vad_fV(:,1:nn)=pg\xsmVa;

%% SET 6 :

% PANEL V :
vb_fV(1:nn,1:nn-1)=funfV(1:nn,1:nn-1,3);

% PANEL III
fun=funfIII(1:nn,nn+1-[1:nn-1],3);
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
vb_fV(1:nn,2*nn-1:3*nn-3)=funfVI(nn:-1:1,1:nn-1,3); 

% PANEL I
fun=funfI(1:nn,1:nn-1,3);
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


funbV=vb_fV(1:nn,:)';
xsmVb=kg*funbV;
vbd_fV(1:nn,:)=(pg\xsmVb)';

%% AFFECTATIONS
dg_xi_fI(1:nn,1:nn,3)=vad_fI([1:nn],[1:nn]);
dg_eta_fI(1:nn,1:nn,3)=vbd_fI([1:nn],[1:nn]);

dg_xi_fII(1:nn,1:nn,3)=vad_fII([1:nn],[1:nn]);
dg_eta_fII(1:nn,1:nn,3)=vbd_fII([1:nn],[1:nn]);

dg_xi_fIII(1:nn,1:nn,3)=vad_fI(2*(nn-1)+[1:nn],nn+1-[1:nn]);
dg_eta_fIII(1:nn,1:nn,3)=-vbd_fI([1:nn],2*(nn-1)+nn-[1:nn]+1); 

dg_xi_fIV(1:nn,1:nn,3)=vad_fII(2*(nn-1)+[1:nn],nn-[1:nn]+1);
dg_eta_fIV(1:nn,1:nn,3)=-vbd_fII([1:nn],2*(nn-1)+nn-[1:nn]+1);

dg_xi_fV(1:nn,1:nn,3)=vad_fV([1:nn],[1:nn]);
dg_eta_fV(1:nn,1:nn,3)=vbd_fV([1:nn],[1:nn]);

dg_xi_fVI(1:nn,1:nn,3)=-vad_fV(2*(nn-1)+nn+1-[1:nn],[1:nn]);
dg_eta_fVI(1:nn,1:nn,3)=vbd_fV(nn-[1:nn]+1,2*(nn-1)+[1:nn]);

%% ETAPE 3 - ASSEMBLAGE vortATIONNEL EN FONCTION DES DERIVEES XI/ETA
%% -------------------------------------------------------------------------------------------------------------------------------------------------------PARTIE A ECRIRE EN MATRICIEL
% 8 - FACE I: 
dx=(gxi_I(1:nn,1:nn,2).*dg_xi_fI(1:nn,1:nn,3)-gxi_I(1:nn,1:nn,3).*dg_xi_fI(1:nn,1:nn,2)).*gr_I(1:nn,1:nn,1);
dy=(gxi_I(1:nn,1:nn,3).*dg_xi_fI(1:nn,1:nn,1)-gxi_I(1:nn,1:nn,1).*dg_xi_fI(1:nn,1:nn,3)).*gr_I(1:nn,1:nn,2);
dz=(gxi_I(1:nn,1:nn,1).*dg_xi_fI(1:nn,1:nn,2)-gxi_I(1:nn,1:nn,2).*dg_xi_fI(1:nn,1:nn,1)).*gr_I(1:nn,1:nn,3);
comp_xi=dx+dy+dz;
dx=(geta_I(1:nn,1:nn,2).*dg_eta_fI(1:nn,1:nn,3)-geta_I(1:nn,1:nn,3).*dg_eta_fI(1:nn,1:nn,2)).*gr_I(1:nn,1:nn,1);
dy=(geta_I(1:nn,1:nn,3).*dg_eta_fI(1:nn,1:nn,1)-geta_I(1:nn,1:nn,1).*dg_eta_fI(1:nn,1:nn,3)).*gr_I(1:nn,1:nn,2);
dz=(geta_I(1:nn,1:nn,1).*dg_eta_fI(1:nn,1:nn,2)-geta_I(1:nn,1:nn,2).*dg_eta_fI(1:nn,1:nn,1)).*gr_I(1:nn,1:nn,3);
comp_eta=dx+dy+dz;
vort_I(1:nn,1:nn)=comp_xi+comp_eta;

% 8 - FACE II: 
dx=(gxi_II(1:nn,1:nn,2).*dg_xi_fII(1:nn,1:nn,3)-gxi_II(1:nn,1:nn,3).*dg_xi_fII(1:nn,1:nn,2)).*gr_II(1:nn,1:nn,1);
dy=(gxi_II(1:nn,1:nn,3).*dg_xi_fII(1:nn,1:nn,1)-gxi_II(1:nn,1:nn,1).*dg_xi_fII(1:nn,1:nn,3)).*gr_II(1:nn,1:nn,2);
dz=(gxi_II(1:nn,1:nn,1).*dg_xi_fII(1:nn,1:nn,2)-gxi_II(1:nn,1:nn,2).*dg_xi_fII(1:nn,1:nn,1)).*gr_II(1:nn,1:nn,3);
comp_xi=dx+dy+dz;
dx=(geta_II(1:nn,1:nn,2).*dg_eta_fII(1:nn,1:nn,3)-geta_II(1:nn,1:nn,3).*dg_eta_fII(1:nn,1:nn,2)).*gr_II(1:nn,1:nn,1);
dy=(geta_II(1:nn,1:nn,3).*dg_eta_fII(1:nn,1:nn,1)-geta_II(1:nn,1:nn,1).*dg_eta_fII(1:nn,1:nn,3)).*gr_II(1:nn,1:nn,2);
dz=(geta_II(1:nn,1:nn,1).*dg_eta_fII(1:nn,1:nn,2)-geta_II(1:nn,1:nn,2).*dg_eta_fII(1:nn,1:nn,1)).*gr_II(1:nn,1:nn,3);
comp_eta=dx+dy+dz;
vort_II(1:nn,1:nn)=comp_xi+comp_eta;

% 8 - FACE III: 
dx=(gxi_III(1:nn,1:nn,2).*dg_xi_fIII(1:nn,1:nn,3)-gxi_III(1:nn,1:nn,3).*dg_xi_fIII(1:nn,1:nn,2)).*gr_III(1:nn,1:nn,1);
dy=(gxi_III(1:nn,1:nn,3).*dg_xi_fIII(1:nn,1:nn,1)-gxi_III(1:nn,1:nn,1).*dg_xi_fIII(1:nn,1:nn,3)).*gr_III(1:nn,1:nn,2);
dz=(gxi_III(1:nn,1:nn,1).*dg_xi_fIII(1:nn,1:nn,2)-gxi_III(1:nn,1:nn,2).*dg_xi_fIII(1:nn,1:nn,1)).*gr_III(1:nn,1:nn,3);
comp_xi=dx+dy+dz;
dx=(geta_III(1:nn,1:nn,2).*dg_eta_fIII(1:nn,1:nn,3)-geta_III(1:nn,1:nn,3).*dg_eta_fIII(1:nn,1:nn,2)).*gr_III(1:nn,1:nn,1);
dy=(geta_III(1:nn,1:nn,3).*dg_eta_fIII(1:nn,1:nn,1)-geta_III(1:nn,1:nn,1).*dg_eta_fIII(1:nn,1:nn,3)).*gr_III(1:nn,1:nn,2);
dz=(geta_III(1:nn,1:nn,1).*dg_eta_fIII(1:nn,1:nn,2)-geta_III(1:nn,1:nn,2).*dg_eta_fIII(1:nn,1:nn,1)).*gr_III(1:nn,1:nn,3);
comp_eta=dx+dy+dz;
vort_III(1:nn,1:nn)=comp_xi+comp_eta;

% 8 - FACE IV:
dx=(gxi_IV(1:nn,1:nn,2).*dg_xi_fIV(1:nn,1:nn,3)-gxi_IV(1:nn,1:nn,3).*dg_xi_fIV(1:nn,1:nn,2)).*gr_IV(1:nn,1:nn,1);
dy=(gxi_IV(1:nn,1:nn,3).*dg_xi_fIV(1:nn,1:nn,1)-gxi_IV(1:nn,1:nn,1).*dg_xi_fIV(1:nn,1:nn,3)).*gr_IV(1:nn,1:nn,2);
dz=(gxi_IV(1:nn,1:nn,1).*dg_xi_fIV(1:nn,1:nn,2)-gxi_IV(1:nn,1:nn,2).*dg_xi_fIV(1:nn,1:nn,1)).*gr_IV(1:nn,1:nn,3);
comp_xi=dx+dy+dz;
dx=(geta_IV(1:nn,1:nn,2).*dg_eta_fIV(1:nn,1:nn,3)-geta_IV(1:nn,1:nn,3).*dg_eta_fIV(1:nn,1:nn,2)).*gr_IV(1:nn,1:nn,1);
dy=(geta_IV(1:nn,1:nn,3).*dg_eta_fIV(1:nn,1:nn,1)-geta_IV(1:nn,1:nn,1).*dg_eta_fIV(1:nn,1:nn,3)).*gr_IV(1:nn,1:nn,2);
dz=(geta_IV(1:nn,1:nn,1).*dg_eta_fIV(1:nn,1:nn,2)-geta_IV(1:nn,1:nn,2).*dg_eta_fIV(1:nn,1:nn,1)).*gr_IV(1:nn,1:nn,3);
comp_eta=dx+dy+dz;
vort_IV(1:nn,1:nn)=comp_xi+comp_eta;

% 8 - FACE V: 
dx=(gxi_V(1:nn,1:nn,2).*dg_xi_fV(1:nn,1:nn,3)-gxi_V(1:nn,1:nn,3).*dg_xi_fV(1:nn,1:nn,2)).*gr_V(1:nn,1:nn,1);
dy=(gxi_V(1:nn,1:nn,3).*dg_xi_fV(1:nn,1:nn,1)-gxi_V(1:nn,1:nn,1).*dg_xi_fV(1:nn,1:nn,3)).*gr_V(1:nn,1:nn,2);
dz=(gxi_V(1:nn,1:nn,1).*dg_xi_fV(1:nn,1:nn,2)-gxi_V(1:nn,1:nn,2).*dg_xi_fV(1:nn,1:nn,1)).*gr_V(1:nn,1:nn,3);
comp_xi=dx+dy+dz;
dx=(geta_V(1:nn,1:nn,2).*dg_eta_fV(1:nn,1:nn,3)-geta_V(1:nn,1:nn,3).*dg_eta_fV(1:nn,1:nn,2)).*gr_V(1:nn,1:nn,1);
dy=(geta_V(1:nn,1:nn,3).*dg_eta_fV(1:nn,1:nn,1)-geta_V(1:nn,1:nn,1).*dg_eta_fV(1:nn,1:nn,3)).*gr_V(1:nn,1:nn,2);
dz=(geta_V(1:nn,1:nn,1).*dg_eta_fV(1:nn,1:nn,2)-geta_V(1:nn,1:nn,2).*dg_eta_fV(1:nn,1:nn,1)).*gr_V(1:nn,1:nn,3);
comp_eta=dx+dy+dz;
vort_V(1:nn,1:nn)=comp_xi+comp_eta;

% 8 - FACE VI: 
dx=(gxi_VI(1:nn,1:nn,2).*dg_xi_fVI(1:nn,1:nn,3)-gxi_VI(1:nn,1:nn,3).*dg_xi_fVI(1:nn,1:nn,2)).*gr_VI(1:nn,1:nn,1);
dy=(gxi_VI(1:nn,1:nn,3).*dg_xi_fVI(1:nn,1:nn,1)-gxi_VI(1:nn,1:nn,1).*dg_xi_fVI(1:nn,1:nn,3)).*gr_VI(1:nn,1:nn,2);
dz=(gxi_VI(1:nn,1:nn,1).*dg_xi_fVI(1:nn,1:nn,2)-gxi_VI(1:nn,1:nn,2).*dg_xi_fVI(1:nn,1:nn,1)).*gr_VI(1:nn,1:nn,3);
comp_xi=dx+dy+dz;
dx=(geta_VI(1:nn,1:nn,2).*dg_eta_fVI(1:nn,1:nn,3)-geta_VI(1:nn,1:nn,3).*dg_eta_fVI(1:nn,1:nn,2)).*gr_VI(1:nn,1:nn,1);
dy=(geta_VI(1:nn,1:nn,3).*dg_eta_fVI(1:nn,1:nn,1)-geta_VI(1:nn,1:nn,1).*dg_eta_fVI(1:nn,1:nn,3)).*gr_VI(1:nn,1:nn,2);
dz=(geta_VI(1:nn,1:nn,1).*dg_eta_fVI(1:nn,1:nn,2)-geta_VI(1:nn,1:nn,2).*dg_eta_fVI(1:nn,1:nn,1)).*gr_VI(1:nn,1:nn,3);
comp_eta=dx+dy+dz;
vort_VI(1:nn,1:nn)=comp_xi+comp_eta;


%% ETAPE 4 - demi-somme + 1/3 somme de la divergence 
uwk_I=vort_I(1:nn,1:nn);uwk_II=vort_II(1:nn,1:nn);uwk_III=vort_III(1:nn,1:nn);
uwk_IV=vort_IV(1:nn,1:nn);uwk_V=vort_V(1:nn,1:nn);uwk_VI=vort_VI(1:nn,1:nn);
[vort_I(1:nn,1:nn),vort_II(1:nn,1:nn),vort_III(1:nn,1:nn),vort_IV(1:nn,1:nn),vort_V(1:nn,1:nn),vort_VI(1:nn,1:nn,1)]=...
    ds103(uwk_I,uwk_II,uwk_III,uwk_IV,uwk_V,uwk_VI,n,nn);
end

