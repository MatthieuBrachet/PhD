function [div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
    div103(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI,n,nn)
global radius;
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;
global gxi_I gxi_II gxi_III gxi_IV gxi_V gxi_VI;
global geta_I geta_II geta_III geta_IV geta_V geta_VI;
global xi eta dxi deta
global pg kg;
global kmat
global eta_c jeta_c xi_c ixi_c
global alfasp betasp gamasp
global m1 m2
global umat1 lmat1

%% ************************************************************************ FRONT and BOTTOM on XI
% FACE I
xxtIa_I=y_fI./x_fI;
yytIa_I=z_fI./x_fI;

% FACE II
xwk1=fun5(x_fII);
xwk2=fun7(x_fII);
xxtIa_II=y_fII.*xwk1;
yytIa_II=z_fII.*xwk2;

% FACE III
xxtIa_III=y_fIII./x_fIII;
yytIa_III=-z_fIII./x_fIII;

% FACE IV
xwk1=fun5(x_fII); 
xwk2=fun7(x_fII); 
xxtIa_IV=y_fIV.*xwk1;
yytIa_IV=z_fIV.*xwk2;

deltatIa_I=(1+xxtIa_I.^2+yytIa_I.^2); 
deltatIa_II=(1+xxtIa_II.^2+yytIa_II.^2); 
deltatIa_III=(1+xxtIa_III.^2+yytIa_III.^2); 
deltatIa_IV=(1+xxtIa_IV.^2+yytIa_IV.^2);
deltabtIa_I=sqrt(deltatIa_I);
deltabtIa_II=sqrt(deltatIa_II);
deltabtIa_III=sqrt(deltatIa_III);
deltabtIa_IV=sqrt(deltatIa_IV);
gtIa_I=(radius^2)*(1+xxtIa_I.^2).*(1+yytIa_I.^2)./(deltabtIa_I.^3);
gtIa_II=(radius^2)*(1+xxtIa_II.^2).*(1+yytIa_II.^2)./(deltabtIa_II.^3);
gtIa_III=(radius^2)*(1+xxtIa_III.^2).*(1+yytIa_III.^2)./(deltabtIa_III.^3);
gtIa_IV=(radius^2)*(1+xxtIa_IV.^2).*(1+yytIa_IV.^2)./(deltabtIa_IV.^3);

gxiIa_II=zeros(nn,nn,3); gxiIa_IV=zeros(nn,nn,3);

% - FACE I -
gxiIa_I=gxi_I;
% - FACE II -
xwk1=fun5(x_fII);
xwk2=fun7(x_fII);
gxiIa_II(1:nn,1:nn,1)= -y_fII.*xwk1.*xwk1;
gxiIa_II(1:nn,1:nn,2)= xwk1;
gxiIa_II(1:nn,1:nn,3)= 0;

gxiIa_II(1:nn,1:nn,1)=gxiIa_II(1:nn,1:nn,1)./(1+(y_fII.*xwk2).^2);
gxiIa_II(1:nn,1:nn,2)=gxiIa_II(1:nn,1:nn,2)./(1+(y_fII.*xwk2).^2);
gxiIa_II(1:nn,1:nn,3)=gxiIa_II(1:nn,1:nn,3)./(1+(y_fII.*xwk2).^2);

% FACE III -
gxiIa_III=gxi_III;
% - FACE IV -
xwk1=fun5(x_fIV);
xwk2=fun7(x_fIV);
gxiIa_IV(1:nn,1:nn,1)= -y_fIV.*xwk1.*xwk1;
gxiIa_IV(1:nn,1:nn,2)= xwk1;
gxiIa_IV(1:nn,1:nn,3)= 0;
%
gxiIa_IV(1:nn,1:nn,1)=gxiIa_IV(1:nn,1:nn,1)./(1+(y_fIV.*xwk2).^2);
gxiIa_IV(1:nn,1:nn,2)=gxiIa_IV(1:nn,1:nn,2)./(1+(y_fIV.*xwk2).^2);
gxiIa_IV(1:nn,1:nn,3)=gxiIa_IV(1:nn,1:nn,3)./(1+(y_fIV.*xwk2).^2);

funfI=mfunfI(1:nn,1:nn,1).*gxiIa_I(1:nn,1:nn,1)+mfunfI(1:nn,1:nn,2).*gxiIa_I(1:nn,1:nn,2)+mfunfI(1:nn,1:nn,3).*gxiIa_I(1:nn,1:nn,3);
funfI=funfI.*gtIa_I;

funfII=mfunfII(1:nn,1:nn,1).*gxiIa_II(1:nn,1:nn,1)+mfunfII(1:nn,1:nn,2).*gxiIa_II(1:nn,1:nn,2)+mfunfII(1:nn,1:nn,3).*gxiIa_II(1:nn,1:nn,3);
funfII=funfII.*gtIa_II(1:nn,1:nn);

funfIII=mfunfIII(1:nn,1:nn,1).*gxiIa_III(1:nn,1:nn,1)+mfunfIII(1:nn,1:nn,2).*gxiIa_III(1:nn,1:nn,2)+mfunfIII(1:nn,1:nn,3).*gxiIa_III(1:nn,1:nn,3);
funfIII=funfIII.*gtIa_III;

funfIV=mfunfIV(1:nn,1:nn,1).*gxiIa_IV(1:nn,1:nn,1)+mfunfIV(1:nn,1:nn,2).*gxiIa_IV(1:nn,1:nn,2)+mfunfIV(1:nn,1:nn,3).*gxiIa_IV(1:nn,1:nn,3);
funfIV=funfIV.*gtIa_IV;



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

%% ************************************************************************ FRONT and BOTTOM on ETA
% FACE I
xxtIb_I=y_fI./x_fI;
yytIb_I=z_fI./x_fI;

% FACE V
xwk1=fun5(x_fV);
xwk2=fun7(x_fV);
xxtIb_V=y_fV.*xwk1;
yytIb_V=z_fV.*xwk2;

% FACE III
xxtIb_III=y_fIII./x_fIII;
yytIb_III=-z_fIII./x_fIII;

% FACE VI
xwk1=fun5(x_fVI);
xwk2=fun7(x_fVI);      
xxtIb_VI=y_fVI.*xwk1;
yytIb_VI=z_fVI.*xwk2;
        

deltatIb_I=(1+xxtIb_I.^2+yytIb_I.^2); 
deltatIb_V=(1+xxtIb_V.^2+yytIb_V.^2); 
deltatIb_III=(1+xxtIb_III.^2+yytIb_III.^2); 
deltatIb_VI=(1+xxtIb_VI.^2+yytIb_VI.^2);

deltabtIb_I=sqrt(deltatIb_I);
deltabtIb_V=sqrt(deltatIb_V);
deltabtIb_III=sqrt(deltatIb_III);
deltabtIb_VI=sqrt(deltatIb_VI);

gtIb_I=(radius^2)*(1+xxtIb_I.^2).*(1+yytIb_I.^2)./(deltabtIb_I.^3);
gtIb_V=(radius^2)*(1+xxtIb_V.^2).*(1+yytIb_V.^2)./(deltabtIb_V.^3);
gtIb_III=(radius^2)*(1+xxtIb_III.^2).*(1+yytIb_III.^2)./(deltabtIb_III.^3);
gtIb_VI=(radius^2)*(1+xxtIb_VI.^2).*(1+yytIb_VI.^2)./(deltabtIb_VI.^3);

getaIb_V=zeros(nn,nn,3); getaIb_VI=zeros(nn,nn,3);

% - FACE I -
getaIb_I=geta_I;
% - FACE V -
xwk1=fun5(x_fV);
xwk2=fun7(x_fV);
getaIb_V(1:nn,1:nn,1)= -z_fV.*xwk1.*xwk2;
getaIb_V(1:nn,1:nn,2)= 0;
getaIb_V(1:nn,1:nn,3)= xwk2;

getaIb_V(1:nn,1:nn,1)=getaIb_V(1:nn,1:nn,1)./(1+(z_fV.*xwk2).^2);
getaIb_V(1:nn,1:nn,2)=getaIb_V(1:nn,1:nn,2)./(1+(z_fV.*xwk2).^2);
getaIb_V(1:nn,1:nn,3)=getaIb_V(1:nn,1:nn,3)./(1+(z_fV.*xwk2).^2);

% - FACE III -
getaIb_III=geta_III;
% - FACE VI -
xwk1=fun5(x_fVI);
xwk2=fun7(x_fVI);
getaIb_VI(1:nn,1:nn,1)= -z_fVI.*xwk1.*xwk2;
getaIb_VI(1:nn,1:nn,2)= 0;
getaIb_VI(1:nn,1:nn,3)= xwk2;

getaIb_VI(1:nn,1:nn,1)=getaIb_VI(1:nn,1:nn,1)./(1+(z_fVI.*xwk2).^2);
getaIb_VI(1:nn,1:nn,2)=getaIb_VI(1:nn,1:nn,2)./(1+(z_fVI.*xwk2).^2);
getaIb_VI(1:nn,1:nn,3)=getaIb_VI(1:nn,1:nn,3)./(1+(z_fVI.*xwk2).^2);


funfI=mfunfI(1:nn,1:nn,1).*getaIb_I(1:nn,1:nn,1)+mfunfI(1:nn,1:nn,2).*getaIb_I(1:nn,1:nn,2)+mfunfI(1:nn,1:nn,3).*getaIb_I(1:nn,1:nn,3);
funfI=funfI.*gtIb_I;

funfV=mfunfV(1:nn,1:nn,1).*getaIb_V(1:nn,1:nn,1)+mfunfV(1:nn,1:nn,2).*getaIb_V(1:nn,1:nn,2)+mfunfV(1:nn,1:nn,3).*getaIb_V(1:nn,1:nn,3);
funfV=funfV.*gtIb_V;

funfIII=mfunfIII(1:nn,1:nn,1).*getaIb_III(1:nn,1:nn,1)+mfunfIII(1:nn,1:nn,2).*getaIb_III(1:nn,1:nn,2)+mfunfIII(1:nn,1:nn,3).*getaIb_III(1:nn,1:nn,3);
funfIII=funfIII.*gtIb_III;

funfVI=mfunfVI(1:nn,1:nn,1).*getaIb_VI(1:nn,1:nn,1)+mfunfVI(1:nn,1:nn,2).*getaIb_VI(1:nn,1:nn,2)+mfunfVI(1:nn,1:nn,3).*getaIb_VI(1:nn,1:nn,3);
funfVI=funfVI.*gtIb_VI;

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

%% *** assemblage *********************************************************
% FACE I
dg_alfa(1:nn,1:nn,1) = vad_fI(1:nn,1:nn);
dg_beta(1:nn,1:nn,1) = vbd_fI(1:nn,1:nn);

% FACE III
dg_alfa(1:nn,1:nn,3)=vad_fI(2*(nn-1)+[1:nn],nn-[1:nn]+1);
dg_beta(1:nn,1:nn,3)=vbd_fI(1:nn,2*(nn-1)+nn-[1:nn]+1); 

gt_I=gtIa_I;
div_fI(1:nn,1:nn)=(dg_alfa(1:nn,1:nn,1) + dg_beta(1:nn,1:nn,1))./gt_I(1:nn,1:nn);

gt_III=gtIa_III; 
div_fIII(1:nn,1:nn)=(dg_alfa(1:nn,1:nn,3) - dg_beta(1:nn,1:nn,3))./gt_III(1:nn,1:nn);

%% ************************************************************************ EAST and WEST on XI
% FACE II
xxtIIa_II=-x_fII./y_fII;
yytIIa_II=z_fII./y_fII;

% FACE III
xwk1=fun5(y_fIII);
xwk2=fun7(y_fIII);
xxtIIa_III=-x_fIII.*xwk1;
yytIIa_III=z_fIII.*xwk2;

% FACE IV
xxtIIa_IV=-x_fIV./y_fIV;
yytIIa_IV=-z_fIV./y_fIV;

% FACE I
xwk1=fun5(y_fI);
xwk2=fun7(y_fI);
xxtIIa_I=-x_fI.*xwk1;
yytIIa_I=z_fI.*xwk2;


deltatIIa_II=(1+xxtIIa_II.^2+yytIIa_II.^2); 
deltatIIa_III=(1+xxtIIa_III.^2+yytIIa_III.^2); 
deltatIIa_IV=(1+xxtIIa_IV.^2+yytIIa_IV.^2); 
deltatIIa_I=(1+xxtIIa_I.^2+yytIIa_I.^2);

deltabtIIa_II=sqrt(deltatIIa_II);
deltabtIIa_III=sqrt(deltatIIa_III);
deltabtIIa_IV=sqrt(deltatIIa_IV);
deltabtIIa_I=sqrt(deltatIIa_I);

gtIIa_II=(radius^2)*(1+xxtIIa_II.^2).*(1+yytIIa_II.^2)./(deltabtIIa_II.^3);
gtIIa_III=(radius^2)*(1+xxtIIa_III.^2).*(1+yytIIa_III.^2)./(deltabtIIa_III.^3);
gtIIa_IV=(radius^2)*(1+xxtIIa_IV.^2).*(1+yytIIa_IV.^2)./(deltabtIIa_IV.^3);
gtIIa_I=(radius^2)*(1+xxtIIa_I.^2).*(1+yytIIa_I.^2)./(deltabtIIa_I.^3);

% - FACE II -
gxiIIa_II=gxi_II;
% - FACE III -
xwk1=fun5(y_fIII);
xwk2=fun7(y_fIII);
gxiIIa_III(1:nn,1:nn,1)= -xwk1;
gxiIIa_III(1:nn,1:nn,2)= x_fIII.*xwk1.*xwk1;
gxiIIa_III(1:nn,1:nn,3)= 0;

gxiIIa_III(1:nn,1:nn,1)=gxiIIa_III(1:nn,1:nn,1)./(1+(x_fIII.*xwk2).^2);
gxiIIa_III(1:nn,1:nn,2)=gxiIIa_III(1:nn,1:nn,2)./(1+(x_fIII.*xwk2).^2);
gxiIIa_III(1:nn,1:nn,3)=gxiIIa_III(1:nn,1:nn,3)./(1+(x_fIII.*xwk2).^2);

% - FACE IV -
gxiIIa_IV=gxi_IV;

% - FACE I -
xwk1=fun5(y_fI);
xwk2=fun7(y_fI);
gxiIIa_I(1:nn,1:nn,1)= -xwk1;
gxiIIa_I(1:nn,1:nn,2)= x_fI.*xwk1.*xwk1;
gxiIIa_I(1:nn,1:nn,3)= 0;

gxiIIa_I(1:nn,1:nn,1)=gxiIIa_I(1:nn,1:nn,1)./(1+(x_fI.*xwk2).^2);
gxiIIa_I(1:nn,1:nn,2)=gxiIIa_I(1:nn,1:nn,2)./(1+(x_fI.*xwk2).^2);
gxiIIa_I(1:nn,1:nn,3)=gxiIIa_I(1:nn,1:nn,3)./(1+(x_fI.*xwk2).^2);


funfII=mfunfII(1:nn,1:nn,1).*gxiIIa_II(1:nn,1:nn,1)+mfunfII(1:nn,1:nn,2).*gxiIIa_II(1:nn,1:nn,2)+mfunfII(1:nn,1:nn,3).*gxiIIa_II(1:nn,1:nn,3);
funfII=funfII.*gtIIa_II;

funfIII=mfunfIII(1:nn,1:nn,1).*gxiIIa_III(1:nn,1:nn,1)+mfunfIII(1:nn,1:nn,2).*gxiIIa_III(1:nn,1:nn,2)+mfunfIII(1:nn,1:nn,3).*gxiIIa_III(1:nn,1:nn,3);
funfIII=funfIII.*gtIIa_III;

funfIV=mfunfIV(1:nn,1:nn,1).*gxiIIa_IV(1:nn,1:nn,1)+mfunfIV(1:nn,1:nn,2).*gxiIIa_IV(1:nn,1:nn,2)+mfunfIV(1:nn,1:nn,3).*gxiIIa_IV(1:nn,1:nn,3);
funfIV=funfIV.*gtIIa_IV;

funfI=mfunfI(1:nn,1:nn,1).*gxiIIa_I(1:nn,1:nn,1)+mfunfI(1:nn,1:nn,2).*gxiIIa_I(1:nn,1:nn,2)+mfunfI(1:nn,1:nn,3).*gxiIIa_I(1:nn,1:nn,3);
funfI=funfI.*gtIIa_I;

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

%% ************************************************************************ EAST and WEST on ETA
% Face II
xxtIIb_II=-x_fII./y_fII;
yytIIb_II=z_fII./y_fII;

% FACE V
xwk=fun7(y_fV);
xxtIIb_V=-x_fV.*xwk;
yytIIb_V=z_fV.*xwk;

% FACE IV
xxtIIb_IV=-x_fIV./y_fIV;
yytIIb_IV=-z_fIV./y_fIV;

% FACE VI
xwk=fun7(y_fVI);
xxtIIb_VI=-x_fVI.*xwk;
yytIIb_VI=z_fVI.*xwk;


deltatIIb_II=(1+xxtIIb_II.^2+yytIIb_II.^2); 
deltatIIb_V=(1+xxtIIb_V.^2+yytIIb_V.^2); 
deltatIIb_IV=(1+xxtIIb_IV.^2+yytIIb_IV.^2); 
deltatIIb_VI=(1+xxtIIb_VI.^2+yytIIb_VI.^2);

deltabtIIb_II=sqrt(deltatIIb_II);
deltabtIIb_V=sqrt(deltatIIb_V);
deltabtIIb_IV=sqrt(deltatIIb_IV);
deltabtIIb_VI=sqrt(deltatIIb_VI);

gtIIb_II=(radius^2)*(1+xxtIIb_II.^2).*(1+yytIIb_II.^2)./(deltabtIIb_II.^3);
gtIIb_V=(radius^2)*(1+xxtIIb_V.^2).*(1+yytIIb_V.^2)./(deltabtIIb_V.^3);
gtIIb_IV=(radius^2)*(1+xxtIIb_IV.^2).*(1+yytIIb_IV.^2)./(deltabtIIb_IV.^3);
gtIIb_VI=(radius^2)*(1+xxtIIb_VI.^2).*(1+yytIIb_VI.^2)./(deltabtIIb_VI.^3);

getaIIb_V=zeros(nn,nn,3); getaIIb_VI=zeros(nn,nn,3);
% - FACE II -
getaIIb_II=geta_II;

% - FACE V -
xwk1=fun5(y_fV);
xwk2=fun7(y_fV);
getaIIb_V(1:nn,1:nn,1)=0;
getaIIb_V(1:nn,1:nn,2)=-z_fV.*xwk1;
getaIIb_V(1:nn,1:nn,3)= 1;

getaIIb_V(1:nn,1:nn,1)=xwk2.*getaIIb_V(1:nn,1:nn,1)./(1+(z_fV.*xwk2).^2);
getaIIb_V(1:nn,1:nn,2)=xwk2.*getaIIb_V(1:nn,1:nn,2)./(1+(z_fV.*xwk2).^2);
getaIIb_V(1:nn,1:nn,3)=xwk2.*getaIIb_V(1:nn,1:nn,3)./(1+(z_fV.*xwk2).^2);
      
% - FACE IV -
getaIIb_IV=geta_IV;

% - FACE VI -
xwk1=fun5(y_fVI);
xwk2=fun7(y_fVI);
getaIIb_VI(1:nn,1:nn,1)= 0;
getaIIb_VI(1:nn,1:nn,2)=-z_fVI.*xwk1.*xwk2;
getaIIb_VI(1:nn,1:nn,3)= xwk2;

getaIIb_VI(1:nn,1:nn,1)=getaIIb_VI(1:nn,1:nn,1)./(1+(z_fVI.*xwk2).^2);
getaIIb_VI(1:nn,1:nn,2)=getaIIb_VI(1:nn,1:nn,2)./(1+(z_fVI.*xwk2).^2);
getaIIb_VI(1:nn,1:nn,3)=getaIIb_VI(1:nn,1:nn,3)./(1+(z_fVI.*xwk2).^2);


funfII=mfunfII(:,:,1).*getaIIb_II(:,:,1)+mfunfII(:,:,2).*getaIIb_II(:,:,2)+mfunfII(:,:,3).*getaIIb_II(:,:,3);
funfII=funfII.*gtIIb_II;

funfV=mfunfV(:,:,1).*getaIIb_V(:,:,1)+mfunfV(:,:,2).*getaIIb_V(:,:,2)+mfunfV(:,:,3).*getaIIb_V(:,:,3);
funfV=funfV.*gtIIb_V;

funfIV=mfunfIV(:,:,1).*getaIIb_IV(:,:,1)+mfunfIV(:,:,2).*getaIIb_IV(:,:,2)+mfunfIV(:,:,3).*getaIIb_IV(:,:,3);
funfIV=funfIV.*gtIIb_IV;

funfVI=mfunfVI(:,:,1).*getaIIb_VI(:,:,1)+mfunfVI(:,:,2).*getaIIb_VI(:,:,2)+mfunfVI(:,:,3).*getaIIb_VI(:,:,3);
funfVI=funfVI.*gtIIb_VI;

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

%% *** Assemblage *********************************************************

% FACE II
dg_alfa(1:nn,1:nn,2)=vad_fII(1:nn,1:nn);
dg_beta(1:nn,1:nn,2)=vbd_fII(1:nn,1:nn);

% FACE IV
dg_alfa(1:nn,1:nn,4)=vad_fII(2*(nn-1)+[1:nn],nn-[1:nn]+1);
dg_beta(1:nn,1:nn,4)=vbd_fII(1:nn,2*(nn-1)+nn-[1:nn]+1);


gt_II=gtIIa_II; 
div_fII=(dg_alfa(1:nn,1:nn,2) + dg_beta(1:nn,1:nn,2))./gt_II;

gt_IV=gtIIa_IV; 
div_fIV=(dg_alfa(1:nn,1:nn,4) - dg_beta(1:nn,1:nn,4))./gt_IV;


%% ************************************************************************ NORTH AND SOUTH on XI
% FACE V
xxtVa_V=y_fV./z_fV;
yytVa_V=-x_fV./z_fV;

% FACE II
xwk1=fun5(z_fII);
xwk2=fun7(z_fII);
xxtVa_II=y_fII.*xwk2;
yytVa_II=-x_fII.*xwk1;

% FACE VI
xxtVa_VI=-y_fVI./z_fVI;
yytVa_VI=-x_fVI./z_fVI;

% FACE IV
xwk1=fun5(z_fIV);
xwk2=fun7(z_fIV);
xxtVa_IV=y_fIV.*xwk2;
yytVa_IV=-x_fIV.*xwk1;


deltatVa_V=(1+xxtVa_V.^2+yytVa_V.^2); 
deltatVa_II=(1+xxtVa_II.^2+yytVa_II.^2); 
deltatVa_VI=(1+xxtVa_VI.^2+yytVa_VI.^2); 
deltatVa_IV=(1+xxtVa_IV.^2+yytVa_IV.^2);

deltabtVa_V=sqrt(deltatVa_V);
deltabtVa_II=sqrt(deltatVa_II);
deltabtVa_VI=sqrt(deltatVa_VI);
deltabtVa_IV=sqrt(deltatVa_IV);

gtVa_V =(radius^2)*(1+xxtVa_V.^2).*(1+yytVa_V.^2)./(deltabtVa_V.^3);
gtVa_II=(radius^2)*(1+xxtVa_II.^2).*(1+yytVa_II.^2)./(deltabtVa_II.^3);
gtVa_VI=(radius^2)*(1+xxtVa_VI.^2).*(1+yytVa_VI.^2)./(deltabtVa_VI.^3);
gtVa_IV=(radius^2)*(1+xxtVa_IV.^2).*(1+yytVa_IV.^2)./(deltabtVa_IV.^3);

gxiVa_II=zeros(nn,nn,3); gxiVa_IV=zeros(nn,nn,3);

% - FACE V -
gxiVa_V=gxi_V;
% - FACE II -
xwk1=fun5(z_fII);
xwk2=fun7(z_fII);
gxiVa_II(1:nn,1:nn,1)= 0;
gxiVa_II(1:nn,1:nn,2)= xwk2;
gxiVa_II(1:nn,1:nn,3)= -y_fII.*xwk1.*xwk2;

gxiVa_II(1:nn,1:nn,1)=gxiVa_II(1:nn,1:nn,1)./(1+(y_fII.*xwk2).^2);
gxiVa_II(1:nn,1:nn,2)=gxiVa_II(1:nn,1:nn,2)./(1+(y_fII.*xwk2).^2);
gxiVa_II(1:nn,1:nn,3)=gxiVa_II(1:nn,1:nn,3)./(1+(y_fII.*xwk2).^2);

% - FACE VI -
gxiVa_VI=gxi_VI;

% - FACE IV -
xwk1=fun5(z_fIV);
xwk2=fun7(z_fIV);
gxiVa_IV(1:nn,1:nn,1)= 0;
gxiVa_IV(1:nn,1:nn,2)= xwk2;
gxiVa_IV(1:nn,1:nn,3)= -y_fIV.*xwk1.*xwk2;
gxiVa_IV(1:nn,1:nn,1)=gxiVa_IV(1:nn,1:nn,1)./(1+(y_fIV.*xwk2).^2);
gxiVa_IV(1:nn,1:nn,2)=gxiVa_IV(1:nn,1:nn,2)./(1+(y_fIV.*xwk2).^2);
gxiVa_IV(1:nn,1:nn,3)=gxiVa_IV(1:nn,1:nn,3)./(1+(y_fIV.*xwk2).^2);


funfV=mfunfV(1:nn,1:nn,1).*gxiVa_V(1:nn,1:nn,1)+mfunfV(1:nn,1:nn,2).*gxiVa_V(1:nn,1:nn,2)+mfunfV(1:nn,1:nn,3).*gxiVa_V(1:nn,1:nn,3);
funfV=funfV.*gtVa_V;

funfII=mfunfII(1:nn,1:nn,1).*gxiVa_II(1:nn,1:nn,1)+mfunfII(1:nn,1:nn,2).*gxiVa_II(1:nn,1:nn,2)+mfunfII(1:nn,1:nn,3).*gxiVa_II(1:nn,1:nn,3);
funfII=funfII.*gtVa_II;

funfVI=mfunfVI(1:nn,1:nn,1).*gxiVa_VI(1:nn,1:nn,1)+mfunfVI(1:nn,1:nn,2).*gxiVa_VI(1:nn,1:nn,2)+mfunfVI(1:nn,1:nn,3).*gxiVa_VI(1:nn,1:nn,3);
funfVI=funfVI.*gtVa_VI;

funfIV=mfunfIV(1:nn,1:nn,1).*gxiVa_IV(1:nn,1:nn,1)+mfunfIV(1:nn,1:nn,2).*gxiVa_IV(1:nn,1:nn,2)+mfunfIV(1:nn,1:nn,3).*gxiVa_IV(1:nn,1:nn,3);
funfIV=funfIV.*gtVa_IV;

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

%% ************************************************************************ NORTH and SOUTH on ETA
% FACE V
xxtVb_V=y_fV./z_fV;
yytVb_V=-x_fV./z_fV;

% FACE III
xwk1=fun5(z_fIII);
xwk2=fun7(z_fIII);
xxtVb_III=y_fIII.*xwk2;
yytVb_III=-x_fIII.*xwk1;

% FACE VI
xxtVb_VI=-y_fVI./z_fVI;
yytVb_VI=-x_fVI./z_fVI;

% FACE I
xwk1=fun5(z_fI);
xwk2=fun7(z_fI);
xxtVb_I=y_fI.*xwk2;
yytVb_I=-x_fI.*xwk1;


deltatVb_V=(1+xxtVb_V.^2+yytVb_V.^2); 
deltatVb_III=(1+xxtVb_III.^2+yytVb_III.^2); 
deltatVb_VI=(1+xxtVb_VI.^2+yytVb_VI.^2); 
deltatVb_I=(1+xxtVb_I.^2+yytVb_I.^2);

deltabtVb_V=sqrt(deltatVb_V);
deltabtVb_III=sqrt(deltatVb_III);
deltabtVb_VI=sqrt(deltatVb_VI);
deltabtVb_I=sqrt(deltatVb_I);

gtVb_V=(radius^2)*(1+xxtVb_V.^2).*(1+yytVb_V.^2)./(deltabtVb_V.^3);
gtVb_III=(radius^2)*(1+xxtVb_III.^2).*(1+yytVb_III.^2)./(deltabtVb_III.^3);
gtVb_VI=(radius^2)*(1+xxtVb_VI.^2).*(1+yytVb_VI.^2)./(deltabtVb_VI.^3);
gtVb_I=(radius^2)*(1+xxtVb_I.^2).*(1+yytVb_I.^2)./(deltabtVb_I.^3);

getaVb_III=zeros(nn,nn,3); getaVb_I=zeros(nn,nn,3);
% - FACE V -
getaVb_V=geta_V;
% - FACE III -
xwk1=fun5(z_fIII);
xwk2=fun7(z_fIII);
getaVb_III(1:nn,1:nn,1)= -xwk1; 
getaVb_III(1:nn,1:nn,2)= 0;
getaVb_III(1:nn,1:nn,3)= x_fIII.*xwk1.*xwk1;

getaVb_III(1:nn,1:nn,1)=getaVb_III(1:nn,1:nn,1)./(1+(x_fIII.*xwk2).^2);
getaVb_III(1:nn,1:nn,2)=getaVb_III(1:nn,1:nn,2)./(1+(x_fIII.*xwk2).^2);
getaVb_III(1:nn,1:nn,3)=getaVb_III(1:nn,1:nn,3)./(1+(x_fIII.*xwk2).^2);

% - FACE VI -
getaVb_VI=geta_VI;
% - FACE I -
xwk1=fun5(z_fI);
xwk2=fun7(z_fI);
getaVb_I(1:nn,1:nn,1)= -xwk1; 
getaVb_I(1:nn,1:nn,2)= 0;
getaVb_I(1:nn,1:nn,3)= x_fI.*xwk1.*xwk1;
%
getaVb_I(1:nn,1:nn,1)=getaVb_I(1:nn,1:nn,1)./(1+(x_fI.*xwk2).^2);
getaVb_I(1:nn,1:nn,2)=getaVb_I(1:nn,1:nn,2)./(1+(x_fI.*xwk2).^2);
getaVb_I(1:nn,1:nn,3)=getaVb_I(1:nn,1:nn,3)./(1+(x_fI.*xwk2).^2);


funfV=mfunfV(1:nn,1:nn,1).*getaVb_V(1:nn,1:nn,1)+mfunfV(1:nn,1:nn,2).*getaVb_V(1:nn,1:nn,2)+mfunfV(1:nn,1:nn,3).*getaVb_V(1:nn,1:nn,3);
funfV=funfV.*gtVb_V;

funfIII=mfunfIII(1:nn,1:nn,1).*getaVb_III(1:nn,1:nn,1)+mfunfIII(1:nn,1:nn,2).*getaVb_III(1:nn,1:nn,2)+mfunfIII(1:nn,1:nn,3).*getaVb_III(1:nn,1:nn,3);
funfIII=funfIII.*gtVb_III;

funfVI=mfunfVI(1:nn,1:nn,1).*getaVb_VI(1:nn,1:nn,1)+mfunfVI(1:nn,1:nn,2).*getaVb_VI(1:nn,1:nn,2)+mfunfVI(1:nn,1:nn,3).*getaVb_VI(1:nn,1:nn,3);
funfVI=funfVI.*gtVb_VI;

funfI=mfunfI(1:nn,1:nn,1).*getaVb_I(1:nn,1:nn,1)+mfunfI(1:nn,1:nn,2).*getaVb_I(1:nn,1:nn,2)+mfunfI(1:nn,1:nn,3).*getaVb_I(1:nn,1:nn,3);
funfI=funfI.*gtVb_I;


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


%% *** Assemblage *********************************************************

% FACE V
dg_alfa(1:nn,1:nn,5)=vad_fV(1:nn,1:nn);
dg_beta(1:nn,1:nn,5)=vbd_fV(1:nn,1:nn);

% FACE VI
dg_alfa(1:nn,1:nn,6)=vad_fV(2*(nn-1)+nn-[1:nn]+1,1:nn);
dg_beta(1:nn,1:nn,6)=vbd_fV(nn-[1:nn]+1,2*(nn-1)+[1:nn]);


gt_V=gtVb_V;
div_fV(1:nn,1:nn)=(dg_alfa(1:nn,1:nn,5) + dg_beta(1:nn,1:nn,5))./gt_V;

gt_VI=gtVa_VI; 
div_fVI(1:nn,1:nn)=(-dg_alfa(1:nn,1:nn,6) + dg_beta(1:nn,1:nn,6))./gt_VI;


%% *** demi-somme + 1/3 somme de la divergence ****************************
uwk_I=div_fI(1:nn,1:nn,1);uwk_II=div_fII(1:nn,1:nn,1);uwk_III=div_fIII(1:nn,1:nn,1);
uwk_IV=div_fIV(1:nn,1:nn,1);uwk_V=div_fV(1:nn,1:nn,1);uwk_VI=div_fVI(1:nn,1:nn,1);
[div_fI(1:nn,1:nn),div_fII(1:nn,1:nn),div_fIII(1:nn,1:nn),div_fIV(1:nn,1:nn),div_fV(1:nn,1:nn),div_fVI(1:nn,1:nn,1)]=...
    ds101(uwk_I,uwk_II,uwk_III,uwk_IV,uwk_V,uwk_VI,n,nn);
