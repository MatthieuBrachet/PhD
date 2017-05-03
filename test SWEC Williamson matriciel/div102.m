function [div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
    div102(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI,n,nn)
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
xxtIa_I=zeros(nn,nn);xxtIa_II=zeros(nn,nn);
xxtIa_III=zeros(nn,nn);xxtIa_IV=zeros(nn,nn);
yytIa_I=zeros(nn,nn);yytIa_II=zeros(nn,nn);
yytIa_III=zeros(nn,nn);yytIa_IV=zeros(nn,nn);

% FACE I
for i=1:nn,
    for j=1:nn,
        xxtIa_I(i,j)=y_fI(i,j)/x_fI(i,j);
        yytIa_I(i,j)=z_fI(i,j)/x_fI(i,j);
    end
end
% FACE II
for i=1:nn,
    for j=1:nn,
        xwk1=fun5(x_fII(i,j));
        xwk2=fun7(x_fII(i,j));
        xxtIa_II(i,j)=y_fII(i,j)*xwk1;
        yytIa_II(i,j)=z_fII(i,j)*xwk2;
    end
end
% FACE III
for i=1:nn,
    for j=1:nn,
        xxtIa_III(i,j)=y_fIII(i,j)/x_fIII(i,j);
        yytIa_III(i,j)=-z_fIII(i,j)/x_fIII(i,j);
    end
end
% FACE IV
for i=1:nn,
    for j=1:nn,
        xwk1=fun5(x_fII(i,j)); 
        xwk2=fun7(x_fII(i,j)); 
        xxtIa_IV(i,j)=y_fIV(i,j)*xwk1;
        yytIa_IV(i,j)=z_fIV(i,j)*xwk2;
    end
end

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
for i=1:nn,
    for j=1:nn,
      xwk1=fun5(x_fII(i,j));
      xwk2=fun7(x_fII(i,j));
      gxiIa_II(i,j,1)= -y_fII(i,j)*xwk1*xwk1;
      gxiIa_II(i,j,2)= xwk1;
      gxiIa_II(i,j,3)= 0;
      
      gxiIa_II(i,j,1:3)=gxiIa_II(i,j,1:3)/(1+(y_fII(i,j)*xwk2)^2);
    end
end
% - FACE III -
gxiIa_III=gxi_III;
% - FACE IV -
for i=1:nn,
    for j=1:nn,
      xwk1=fun5(x_fIV(i,j));
      xwk2=fun7(x_fIV(i,j));
      gxiIa_IV(i,j,1)= -y_fIV(i,j)*xwk1*xwk1;
      gxiIa_IV(i,j,2)= xwk1;
      gxiIa_IV(i,j,3)= 0;
      %
      gxiIa_IV(i,j,1:3)=gxiIa_IV(i,j,1:3)/(1+(y_fIV(i,j)*xwk2)^2);
    end
end

funfI=zeros(nn,nn);funfII=zeros(nn,nn);funfIII=zeros(nn,nn);
funfIV=zeros(nn,nn);funfV=zeros(nn,nn);funfVI=zeros(nn,nn);
%
for i=1:nn,
    for j=1:nn,
      funfI(i,j)=dot(mfunfI(i,j,:),gxiIa_I(i,j,:));
      funfI(i,j)=funfI(i,j)*gtIa_I(i,j);
     end
end
for i=1:nn,
    for j=1:nn,
      funfII(i,j)=dot(mfunfII(i,j,:),gxiIa_II(i,j,:));
      funfII(i,j)=funfII(i,j)*gtIa_II(i,j);
     end
end
for i=1:nn,
    for j=1:nn,
      funfIII(i,j)=dot(mfunfIII(i,j,:),gxiIa_III(i,j,:));
      funfIII(i,j)=funfIII(i,j)*gtIa_III(i,j);
     end
end
for i=1:nn,
    for j=1:nn,
      funfIV(i,j)=dot(mfunfIV(i,j,:),gxiIa_IV(i,j,:));
      funfIV(i,j)=funfIV(i,j)*gtIa_IV(i,j);
     end
end

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
xxtIb_I=zeros(nn,nn);xxtIb_V=zeros(nn,nn);
xxtIb_III=zeros(nn,nn);xxtIb_VI=zeros(nn,nn);
yytIb_I=zeros(nn,nn);yytIb_V=zeros(nn,nn);
yytIb_III=zeros(nn,nn);yytIb_VI=zeros(nn,nn);

% FACE I
for i=1:nn,
    for j=1:nn,
        xxtIb_I(i,j)=y_fI(i,j)/x_fI(i,j);
        yytIb_I(i,j)=z_fI(i,j)/x_fI(i,j);
    end
end
% FACE V
for i=1:nn,
    for j=1:nn,
        xwk1=fun5(x_fV(i,j));
        xwk2=fun7(x_fV(i,j));
        xxtIb_V(i,j)=y_fV(i,j)*xwk1;
        yytIb_V(i,j)=z_fV(i,j)*xwk2;
    end
end
% FACE III
for i=1:nn,
    for j=1:nn,
        xxtIb_III(i,j)=y_fIII(i,j)/x_fIII(i,j);
        yytIb_III(i,j)=-z_fIII(i,j)/x_fIII(i,j);
    end
end
% FACE VI
for i=1:nn,
    for j=1:nn,
        xwk1=fun5(x_fVI(i,j));
        xwk2=fun7(x_fVI(i,j));      
        xxtIb_VI(i,j)=y_fVI(i,j)*xwk1;
        yytIb_VI(i,j)=z_fVI(i,j)*xwk2;
    end
end

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
for i=1:nn,
    for j=1:nn,
      xwk1=fun5(x_fV(i,j));
      xwk2=fun7(x_fV(i,j));
      getaIb_V(i,j,1)= -z_fV(i,j)*xwk1*xwk2;
      getaIb_V(i,j,2)= 0;
      getaIb_V(i,j,3)= xwk2;

      getaIb_V(i,j,1:3)=getaIb_V(i,j,1:3)/(1+(z_fV(i,j)*xwk2)^2);
    end
end
% - FACE III -
getaIb_III=geta_III;
% - FACE VI -
for i=1:nn,
    for j=1:nn,
      xwk1=fun5(x_fVI(i,j));
      xwk2=fun7(x_fVI(i,j));
      getaIb_VI(i,j,1)= -z_fVI(i,j)*xwk1*xwk2;
      getaIb_VI(i,j,2)= 0;
      getaIb_VI(i,j,3)= xwk2;

      getaIb_VI(i,j,1:3)=getaIb_VI(i,j,1:3)/(1+(z_fVI(i,j)*xwk2)^2);
    end
end

for i=1:nn,
    for j=1:nn,
      funfI(i,j)=dot(mfunfI(i,j,:),getaIb_I(i,j,:));
      funfI(i,j)=funfI(i,j)*gtIb_I(i,j);
     end
end
for i=1:nn,
    for j=1:nn,
      funfV(i,j)=dot(mfunfV(i,j,:),getaIb_V(i,j,:));
      funfV(i,j)=funfV(i,j)*gtIb_V(i,j);
     end
end
for i=1:nn,
    for j=1:nn,
      funfIII(i,j)=dot(mfunfIII(i,j,:),getaIb_III(i,j,:));
      funfIII(i,j)=funfIII(i,j)*gtIb_III(i,j);
     end
end
for i=1:nn,
    for j=1:nn,
      funfVI(i,j)=dot(mfunfVI(i,j,:),getaIb_VI(i,j,:));
      funfVI(i,j)=funfVI(i,j)*gtIb_VI(i,j);
     end
end

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
dg_alfa=zeros(nn,nn,6);dg_beta=zeros(nn,nn,6);
% FACE I
for i=1:nn,
    for j=1:nn,
        dg_alfa(i,j,1) = vad_fI(i,j);
        dg_beta(i,j,1) = vbd_fI(i,j);
    end
end
% FACE III
for i=1:nn,
    for j=1:nn,
        dg_alfa(i,j,3)=vad_fI(2*(nn-1)+i,nn-j+1);
        dg_beta(i,j,3)=vbd_fI(i,2*(nn-1)+nn-j+1); 
    end
end

div_fI=zeros(nn,nn);

gt_I=gtIa_I;
for i=1:nn,
    for j=1:nn,
        div_fI(i,j)=(dg_alfa(i,j,1) + dg_beta(i,j,1))/gt_I(i,j);
    end
end

div_fIII=zeros(nn,nn);

gt_III=gtIa_III; 
for i=1:nn,
    for j=1:nn,
        div_fIII(i,j)=(dg_alfa(i,j,3) - dg_beta(i,j,3))/gt_III(i,j);
    end
end

%% ************************************************************************ EAST and WEST on XI
xxtIIa_II=zeros(nn,nn);xxtIIa_III=zeros(nn,nn);
xxtIIa_IV=zeros(nn,nn);xxtIIa_I=zeros(nn,nn);
yytIIa_II=zeros(nn,nn);yytIIa_III=zeros(nn,nn);
yytIIa_IV=zeros(nn,nn);yytIIa_I=zeros(nn,nn);

% FACE II
for i=1:nn,
    for j=1:nn,
        xxtIIa_II(i,j)=-x_fII(i,j)/y_fII(i,j);
        yytIIa_II(i,j)=z_fII(i,j)/y_fII(i,j);
    end
end
% FACE III
xwk1=fun5(y_fIII);
xwk2=fun7(y_fIII);
for i=1:nn,
    for j=1:nn,
        xxtIIa_III(i,j)=-x_fIII(i,j)*xwk1(i,j);
        yytIIa_III(i,j)=z_fIII(i,j)*xwk2(i,j);
    end
end
% FACE IV
for i=1:nn,
    for j=1:nn,
        xxtIIa_IV(i,j)=-x_fIV(i,j)/y_fIV(i,j);
        yytIIa_IV(i,j)=-z_fIV(i,j)/y_fIV(i,j);
    end
end
% FACE I
xwk1=fun5(y_fI);
xwk2=fun7(y_fI);
for i=1:nn,
    for j=1:nn,
        xxtIIa_I(i,j)=-x_fI(i,j)*xwk1(i,j);
        yytIIa_I(i,j)=z_fI(i,j)*xwk2(i,j);
    end
end

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

gxiIIa_III=zeros(nn,nn,3); gxiIIa_I=zeros(nn,nn,3);
% - FACE II -
gxiIIa_II=gxi_II;
% - FACE III -
for i=1:nn,
    for j=1:nn,
      xwk1=fun5(y_fIII(i,j));
      xwk2=fun7(y_fIII(i,j));
      gxiIIa_III(i,j,1)= -xwk1;
      gxiIIa_III(i,j,2)= x_fIII(i,j)*xwk1*xwk1;
      gxiIIa_III(i,j,3)= 0;

      gxiIIa_III(i,j,1:3)=gxiIIa_III(i,j,1:3)/(1+(x_fIII(i,j)*xwk2)^2);
    end
end
% - FACE IV -
gxiIIa_IV=gxi_IV;
% - FACE I -
for i=1:nn,
    for j=1:nn,
      xwk1=fun5(y_fI(i,j));
      xwk2=fun7(y_fI(i,j));
      gxiIIa_I(i,j,1)= -xwk1;
      gxiIIa_I(i,j,2)= x_fI(i,j)*xwk1*xwk1;
      gxiIIa_I(i,j,3)= 0;

      gxiIIa_I(i,j,1:3)=gxiIIa_I(i,j,1:3)/(1+(x_fI(i,j)*xwk2)^2);
    end
end

funfI=zeros(nn,nn);funfII=zeros(nn,nn);funfIII=zeros(nn,nn);
funfIV=zeros(nn,nn);funfV=zeros(nn,nn);funfVI=zeros(nn,nn);

for i=1:nn,
    for j=1:nn,
      funfII(i,j)=dot(mfunfII(i,j,:),gxiIIa_II(i,j,:));
      funfII(i,j)=funfII(i,j)*gtIIa_II(i,j);
     end
end
for i=1:nn,
    for j=1:nn,
      funfIII(i,j)=dot(mfunfIII(i,j,:),gxiIIa_III(i,j,:));
      funfIII(i,j)=funfIII(i,j)*gtIIa_III(i,j);
     end
end
for i=1:nn,
    for j=1:nn,
      funfIV(i,j)=dot(mfunfIV(i,j,:),gxiIIa_IV(i,j,:));
      funfIV(i,j)=funfIV(i,j)*gtIIa_IV(i,j);
     end
end
for i=1:nn,
    for j=1:nn,
      funfI(i,j)=dot(mfunfI(i,j,:),gxiIIa_I(i,j,:));
      funfI(i,j)=funfI(i,j)*gtIIa_I(i,j);
     end
end

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
xxtIIb_II=zeros(nn,nn);xxtIIb_V=zeros(nn,nn);
xxtIIb_IV=zeros(nn,nn);xxtIIb_VI=zeros(nn,nn);
yytIIb_II=zeros(nn,nn);yytIIb_V=zeros(nn,nn);
yytIIb_IV=zeros(nn,nn);yytIIb_VI=zeros(nn,nn);

% FACE II
for i=1:nn,
    for j=1:nn,
        xxtIIb_II(i,j)=-x_fII(i,j)/y_fII(i,j);
        yytIIb_II(i,j)=z_fII(i,j)/y_fII(i,j);
    end
end
% FACE V
xwk=fun7(y_fV);
for i=1:nn,
    for j=1:nn,
        xxtIIb_V(i,j)=-x_fV(i,j)*xwk(i,j);
        yytIIb_V(i,j)=z_fV(i,j)*xwk(i,j);
    end
end
% FACE IV
for i=1:nn,
    for j=1:nn,
        xxtIIb_IV(i,j)=-x_fIV(i,j)/y_fIV(i,j);
        yytIIb_IV(i,j)=-z_fIV(i,j)/y_fIV(i,j);
    end
end
% FACE VI
xwk=fun7(y_fVI);
for i=1:nn,
    for j=1:nn,
        xxtIIb_VI(i,j)=-x_fVI(i,j)*xwk(i,j);
        yytIIb_VI(i,j)=z_fVI(i,j)*xwk(i,j);
    end
end

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
for i=1:nn,
    for j=1:nn,
      xwk1=fun5(y_fV(i,j));
      xwk2=fun7(y_fV(i,j));
      getaIIb_V(i,j,1)=0;
      getaIIb_V(i,j,2)=-z_fV(i,j)*xwk1;
      getaIIb_V(i,j,3)= 1;
      
      getaIIb_V(i,j,1:3)=xwk2*getaIIb_V(i,j,1:3)/(1+(z_fV(i,j)*xwk2)^2);
    end
end
% - FACE IV -
getaIIb_IV=geta_IV;
% - FACE VI -
for i=1:nn,
    for j=1:nn,
      xwk1=fun5(y_fVI(i,j));
      xwk2=fun7(y_fVI(i,j));
      getaIIb_VI(i,j,1)= 0;
      getaIIb_VI(i,j,2)=-z_fVI(i,j)*xwk1*xwk2;
      getaIIb_VI(i,j,3)= xwk2;
      %
      getaIIb_VI(i,j,1:3)=getaIIb_VI(i,j,1:3)/(1+(z_fVI(i,j)*xwk2)^2);
    end
end

for i=1:nn,
    for j=1:nn,
      funfII(i,j)=dot(mfunfII(i,j,:),getaIIb_II(i,j,:));
      funfII(i,j)=funfII(i,j)*gtIIb_II(i,j);
     end
end
for i=1:nn,
    for j=1:nn,
      funfV(i,j)=dot(mfunfV(i,j,:),getaIIb_V(i,j,:));
      funfV(i,j)=funfV(i,j)*gtIIb_V(i,j);
     end
end
for i=1:nn,
    for j=1:nn,
      funfIV(i,j)=dot(mfunfIV(i,j,:),getaIIb_IV(i,j,:));
      funfIV(i,j)=funfIV(i,j)*gtIIb_IV(i,j);
     end
end
for i=1:nn,
    for j=1:nn,
      funfVI(i,j)=dot(mfunfVI(i,j,:),getaIIb_VI(i,j,:));
      funfVI(i,j)=funfVI(i,j)*gtIIb_VI(i,j);
     end
end

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
for i=1:nn,
    for j=1:nn,
        dg_alfa(i,j,2)=vad_fII(i,j);
        dg_beta(i,j,2)=vbd_fII(i,j);
    end
end
% FACE IV
for i=1:nn,
    for j=1:nn,
        dg_alfa(i,j,4)=vad_fII(2*(nn-1)+i,nn-j+1);
        dg_beta(i,j,4)=vbd_fII(i,2*(nn-1)+nn-j+1);
    end
end

div_fII=zeros(nn,nn);

gt_II=gtIIa_II; 
for i=1:nn,
    for j=1:nn,
        div_fII(i,j)=(dg_alfa(i,j,2) + dg_beta(i,j,2))/gt_II(i,j);
    end
end

div_fIV=zeros(nn,nn);
gt_IV=gtIIa_IV; 
for i=1:nn,
    for j=1:nn,
        div_fIV(i,j)=(dg_alfa(i,j,4) - dg_beta(i,j,4))/gt_IV(i,j);
    end
end

%% ************************************************************************ NORTH AND SOUTH on XI
xxtVa_V=zeros(nn,nn);xxtVa_II=zeros(nn,nn);
xxtVa_VI=zeros(nn,nn);xxtVa_IV=zeros(nn,nn);

yytVa_V=zeros(nn,nn);yytVa_II=zeros(nn,nn);
yytVa_VI=zeros(nn,nn);yytVa_IV=zeros(nn,nn);

% FACE V
for i=1:nn,
    for j=1:nn,
        xxtVa_V(i,j)=y_fV(i,j)/z_fV(i,j);
        yytVa_V(i,j)=-x_fV(i,j)/z_fV(i,j);
    end
end
% FACE II
for i=1:nn,
    for j=1:nn,
        xwk1=fun5(z_fII(i,j));
        xwk2=fun7(z_fII(i,j));
        xxtVa_II(i,j)=y_fII(i,j)*xwk2;
        yytVa_II(i,j)=-x_fII(i,j)*xwk1;
    end
end
% FACE VI
for i=1:nn,
    for j=1:nn,
        xxtVa_VI(i,j)=-y_fVI(i,j)/z_fVI(i,j);
        yytVa_VI(i,j)=-x_fVI(i,j)/z_fVI(i,j);
    end
end
% FACE IV
for i=1:nn,
    for j=1:nn,
        xwk1=fun5(z_fIV(i,j));
        xwk2=fun7(z_fIV(i,j));
        xxtVa_IV(i,j)=y_fIV(i,j)*xwk2;
        yytVa_IV(i,j)=-x_fIV(i,j)*xwk1;
    end
end

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
for i=1:nn,
    for j=1:nn,
      xwk1=fun5(z_fII(i,j));
      xwk2=fun7(z_fII(i,j));
      gxiVa_II(i,j,1)= 0;
      gxiVa_II(i,j,2)= xwk2;
      gxiVa_II(i,j,3)= -y_fII(i,j)*xwk1*xwk2;
      gxiVa_II(i,j,1:3)=gxiVa_II(i,j,1:3)/(1+(y_fII(i,j)*xwk2)^2);
    end
end
% - FACE VI -
gxiVa_VI=gxi_VI;
% - FACE IV -
for i=1:nn,
    for j=1:nn,
      xwk1=fun5(z_fIV(i,j));
      xwk2=fun7(z_fIV(i,j));
      gxiVa_IV(i,j,1)= 0;
      gxiVa_IV(i,j,2)= xwk2;
      gxiVa_IV(i,j,3)= -y_fIV(i,j)*xwk1*xwk2;
      gxiVa_IV(i,j,1:3)=gxiVa_IV(i,j,1:3)/(1+(y_fIV(i,j)*xwk2)^2);
    end
end

funfI=zeros(nn,nn);funfII=zeros(nn,nn);funfIII=zeros(nn,nn);
funfIV=zeros(nn,nn);funfV=zeros(nn,nn);funfVI=zeros(nn,nn);

for i=1:nn,
    for j=1:nn,
      funfV(i,j)=dot(mfunfV(i,j,:),gxiVa_V(i,j,:));
      funfV(i,j)=funfV(i,j)*gtVa_V(i,j);
     end
end
for i=1:nn,
    for j=1:nn,
      funfII(i,j)=dot(mfunfII(i,j,:),gxiVa_II(i,j,:));
      funfII(i,j)=funfII(i,j)*gtVa_II(i,j);
     end
end
for i=1:nn,
    for j=1:nn,
      funfVI(i,j)=dot(mfunfVI(i,j,:),gxiVa_VI(i,j,:));
      funfVI(i,j)=funfVI(i,j)*gtVa_VI(i,j);
     end
end
for i=1:nn,
    for j=1:nn,
      funfIV(i,j)=dot(mfunfIV(i,j,:),gxiVa_IV(i,j,:));
      funfIV(i,j)=funfIV(i,j)*gtVa_IV(i,j);
     end
end

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
xxtVb_V=zeros(nn,nn);xxtVb_III=zeros(nn,nn);
xxtVb_VI=zeros(nn,nn);xxtVb_I=zeros(nn,nn);

yytVb_V=zeros(nn,nn);yytVb_III=zeros(nn,nn);
yytVb_VI=zeros(nn,nn);yytVb_I=zeros(nn,nn);

% FACE V
for i=1:nn,
    for j=1:nn,
        xxtVb_V(i,j)=y_fV(i,j)/z_fV(i,j);
        yytVb_V(i,j)=-x_fV(i,j)/z_fV(i,j);
    end
end
% FACE III
for i=1:nn,
    for j=1:nn,
        xwk1=fun5(z_fIII(i,j));
        xwk2=fun7(z_fIII(i,j));
        xxtVb_III(i,j)=y_fIII(i,j)*xwk2;
        yytVb_III(i,j)=-x_fIII(i,j)*xwk1;
    end
end
% FACE VI
for i=1:nn,
    for j=1:nn,
        xxtVb_VI(i,j)=-y_fVI(i,j)/z_fVI(i,j);
        yytVb_VI(i,j)=-x_fVI(i,j)/z_fVI(i,j);
    end
end
% FACE I
for i=1:nn,
    for j=1:nn,
        xwk1=fun5(z_fI(i,j));
        xwk2=fun7(z_fI(i,j));
        xxtVb_I(i,j)=y_fI(i,j)*xwk2;
        yytVb_I(i,j)=-x_fI(i,j)*xwk1;
    end
end

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
for i=1:nn,
    for j=1:nn,
      xwk1=fun5(z_fIII(i,j));
      xwk2=fun7(z_fIII(i,j));
      getaVb_III(i,j,1)= -xwk1; 
      getaVb_III(i,j,2)= 0;
      getaVb_III(i,j,3)= x_fIII(i,j)*xwk1*xwk1;

      getaVb_III(i,j,1:3)=getaVb_III(i,j,1:3)/(1+(x_fIII(i,j)*xwk2)^2);
    end
end
% - FACE VI -
getaVb_VI=geta_VI;
% - FACE I -
for i=1:nn,
    for j=1:nn,
      xwk1=fun5(z_fI(i,j));
      xwk2=fun7(z_fI(i,j));
      getaVb_I(i,j,1)= -xwk1; 
      getaVb_I(i,j,2)= 0;
      getaVb_I(i,j,3)= x_fI(i,j)*xwk1*xwk1;
      %
      getaVb_I(i,j,1:3)=getaVb_I(i,j,1:3)/(1+(x_fI(i,j)*xwk2)^2);
    end
end

for i=1:nn,
    for j=1:nn,
      funfV(i,j)=dot(mfunfV(i,j,:),getaVb_V(i,j,:));
      funfV(i,j)=funfV(i,j)*gtVb_V(i,j);
     end
end
for i=1:nn,
    for j=1:nn,
      funfIII(i,j)=dot(mfunfIII(i,j,:),getaVb_III(i,j,:));
      funfIII(i,j)=funfIII(i,j)*gtVb_III(i,j);
     end
end
for i=1:nn,
    for j=1:nn,
      funfVI(i,j)=dot(mfunfVI(i,j,:),getaVb_VI(i,j,:));
      funfVI(i,j)=funfVI(i,j)*gtVb_VI(i,j);
     end
end
for i=1:nn,
    for j=1:nn,
      funfI(i,j)=dot(mfunfI(i,j,:),getaVb_I(i,j,:));
      funfI(i,j)=funfI(i,j)*gtVb_I(i,j);
     end
end

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
for i=1:nn,
    for j=1:nn,
        dg_alfa(i,j,5)=vad_fV(i,j);
        dg_beta(i,j,5)=vbd_fV(i,j);
    end
end
% FACE VI
for i=1:nn,
    for j=1:nn,
        dg_alfa(i,j,6)=vad_fV(2*(nn-1)+nn-i+1,j);
        dg_beta(i,j,6)=vbd_fV(nn-i+1,2*(nn-1)+j);
    end
end

div_fV=zeros(nn,nn);
gt_V=gtVb_V;
for i=1:nn,
    for j=1:nn,
        div_fV(i,j)=(dg_alfa(i,j,5) + dg_beta(i,j,5))/gt_V(i,j);
    end
end

div_fVI=zeros(nn,nn);
gt_VI=gtVa_VI; 
for i=1:nn,
    for j=1:nn,
        div_fVI(i,j)=(-dg_alfa(i,j,6) + dg_beta(i,j,6))/gt_VI(i,j);
    end
end

%% *** demi-somme + 1/3 somme de la divergence ****************************
uwk_I=div_fI(1:nn,1:nn,1);uwk_II=div_fII(1:nn,1:nn,1);uwk_III=div_fIII(1:nn,1:nn,1);
uwk_IV=div_fIV(1:nn,1:nn,1);uwk_V=div_fV(1:nn,1:nn,1);uwk_VI=div_fVI(1:nn,1:nn,1);
[div_fI(1:nn,1:nn),div_fII(1:nn,1:nn),div_fIII(1:nn,1:nn),div_fIV(1:nn,1:nn),div_fV(1:nn,1:nn),div_fVI(1:nn,1:nn,1)]=...
    ds101(uwk_I,uwk_II,uwk_III,uwk_IV,uwk_V,uwk_VI,n,nn);
