function [div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
    div72(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI,n,nn)
global na nb;
global radius;
%global xx yy deltab;
global alfa beta;
global betacr;
global alfa1;
%global alfag betag;
global p k;
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;
global gxi_I gxi_II gxi_III gxi_IV gxi_V gxi_VI;
global geta_I geta_II geta_III geta_IV geta_V geta_VI;

%% ************************************************************************1

xwk1=fun5(x_fII);
xwk2=fun7(x_fII);
xxtIa_I=y_fI./x_fI;
yytIa_I=z_fI./x_fI;
xxtIa_II=y_fII.*xwk1;
yytIa_II=z_fII.*xwk2;
xxtIa_III=y_fIII./x_fIII;
yytIa_III=-z_fIII./x_fIII;
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

gxiIa_II=zeros(nn,nn,3);
gxiIa_IV=zeros(nn,nn,3);

gxiIa_I=gxi_I;

for i=1:nn,
    for j=1:nn,
      gxiIa_II(i,j,1)= -y_fII(i,j)*xwk1(i,j)*xwk1(i,j);
      gxiIa_II(i,j,2)= xwk1(i,j);
      gxiIa_II(i,j,3)= 0;
      gxiIa_II(i,j,1:3)=gxiIa_II(i,j,1:3)/(1+(y_fII(i,j)*xwk2(i,j))^2);
    end
end

gxiIa_III=gxi_III;

xwk1=fun5(x_fIV);
xwk2=fun7(x_fIV);
for i=1:nn,
    for j=1:nn,
      gxiIa_IV(i,j,1)= -y_fIV(i,j).*xwk1(i,j).*xwk1(i,j);
      gxiIa_IV(i,j,2)= xwk1(i,j);
      gxiIa_IV(i,j,3)= 0;
      gxiIa_IV(i,j,1:3)=gxiIa_IV(i,j,1:3)/(1+(y_fIV(i,j)*xwk2(i,j))^2);
    end
end

funfI=zeros(nn,nn);funfII=zeros(nn,nn);funfIII=zeros(nn,nn);
funfIV=zeros(nn,nn);funfV=zeros(nn,nn);funfVI=zeros(nn,nn);

for i=1:nn,
    for j=1:nn,
      funfI(i,j)=dot(mfunfI(i,j,:),gxiIa_I(i,j,:));
      funfI(i,j)=funfI(i,j)*gtIa_I(i,j);
      funfII(i,j)=dot(mfunfII(i,j,:),gxiIa_II(i,j,:));
      funfII(i,j)=funfII(i,j)*gtIa_II(i,j);
      funfIII(i,j)=dot(mfunfIII(i,j,:),gxiIa_III(i,j,:));
      funfIII(i,j)=funfIII(i,j)*gtIa_III(i,j);
      funfIV(i,j)=dot(mfunfIV(i,j,:),gxiIa_IV(i,j,:));
      funfIV(i,j)=funfIV(i,j)*gtIa_IV(i,j);
     end
end

va_fI=zeros(4*(nn-1),nn);
funbII1=zeros(nn,1);
funbIV1=zeros(nn,1);

for jline1=1:nn, 
    va_fI(1:nn-1,jline1)=funfI(1:nn-1,jline1);
end

for i=1:nn-1 
    betaspline(1:nn)=beta(i,1:nn);
    funspline(1:nn)=funfII(i,1:nn);
    ppspline=spline(betaspline,funspline);
    funbII1(1:nn)=ppval(ppspline,betacr(i,1:nn));
    va_fI(nn-1+i,1:nn)=funbII1(1:nn);
end

for jline1=1:nn,
    va_fI(2*nn-1:3*nn-3,jline1)=funfIII(1:nn-1,nn-jline1+1);
end

for i=1:nn-1, 
   betaspline(1:nn)=beta(i,1:nn);
   funspline(1:nn)=funfIV(i,1:nn);
   ppspline=spline(betaspline,funspline);
   funbIV1(1:nn)=ppval(ppspline,betacr(i,nn+1-1:nn));
   va_fI(3*nn-3+i,1:nn)=funbIV1(1:nn);
end

vad_fI=zeros(na,nn);
for jline1=1:nn,
%     alfag4=alfag(:,jline1);
%     dalfag4=k*alfag4;
%     dalfag4(1)=dalfag4(1)+2*pi;
%     dalfag4(na)=dalfag4(na)+2*pi;
%     alfag4d=(p\dalfag4);
    funa7=va_fI(:,jline1); 
    funad7=p\(k*funa7); 
    funad8=funad7;%./alfag4d; 
    vad_fI(:,jline1)=funad8; 
end

%% ************************************************************************


xxtIb_I=y_fI./x_fI;
yytIb_I=z_fI./x_fI;
xwk1=fun5(x_fV);
xwk2=fun7(x_fV);
xxtIb_V=y_fV.*xwk1;
yytIb_V=z_fV.*xwk2;
xxtIb_III=y_fIII./x_fIII;
yytIb_III=-z_fIII./x_fIII;
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

getaIb_V=zeros(nn,nn,3);
getaIb_VI=zeros(nn,nn,3);

getaIb_I=geta_I;

xwk1=fun5(x_fV);
xwk2=fun7(x_fV);
for i=1:nn,
    for j=1:nn,
      getaIb_V(i,j,1)= -z_fV(i,j)*xwk1(i,j)*xwk2(i,j);
      getaIb_V(i,j,2)= 0;
      getaIb_V(i,j,3)= xwk2(i,j);
      getaIb_V(i,j,1:3)=getaIb_V(i,j,1:3)/(1+(z_fV(i,j)*xwk2(i,j))^2);
    end
end

getaIb_III=geta_III;

xwk1=fun5(x_fVI);
xwk2=fun7(x_fVI);
for i=1:nn,
    for j=1:nn,
      getaIb_VI(i,j,1)= -z_fVI(i,j)*xwk1(i,j)*xwk2(i,j);
      getaIb_VI(i,j,2)= 0;
      getaIb_VI(i,j,3)= xwk2(i,j);
      getaIb_VI(i,j,1:3)=getaIb_VI(i,j,1:3)/(1+(z_fVI(i,j)*xwk2(i,j))^2);
      funfI(i,j)=dot(mfunfI(i,j,:),getaIb_I(i,j,:));
      funfI(i,j)=funfI(i,j)*gtIb_I(i,j);
      funfV(i,j)=dot(mfunfV(i,j,:),getaIb_V(i,j,:));
      funfV(i,j)=funfV(i,j)*gtIb_V(i,j);
      funfIII(i,j)=dot(mfunfIII(i,j,:),getaIb_III(i,j,:));
      funfIII(i,j)=funfIII(i,j)*gtIb_III(i,j);
      funfVI(i,j)=dot(mfunfVI(i,j,:),getaIb_VI(i,j,:));
      funfVI(i,j)=funfVI(i,j)*gtIb_VI(i,j);
     end
end

vb_fI=zeros(nn,4*(nn-1)); 
funaV1=zeros(nn,1);
funaVI1=zeros(nn,1);

for iline1=1:nn, 
    vb_fI(iline1,1:nn-1)=funfI(iline1,1:nn-1);
end

for j=1:nn-1 
    alfaspline(1:nn)=alfa(1:nn,j);
    funspline(1:nn)=funfV(1:nn,j);
    ppspline=spline(alfaspline,funspline);
    funaV1(1:nn)=ppval(ppspline,alfa1(1:nn,j));
    vb_fI(1:nn,nn-1+j)=funaV1(1:nn);
end

 for iline1=1:nn,
     vb_fI(iline1,2*nn-1:3*nn-3)=funfIII(iline1,nn:-1:2); 
 end

 for j=1:nn-1
    alfaspline(1:nn)=alfa(1:nn,j);
    funspline(1:nn)=funfVI(1:nn,j);
    ppspline=spline(alfaspline,funspline);
    funaVI1(1:nn)=ppval(ppspline,alfa1(nn+1-1:nn,j));
    vb_fI(1:nn,3*nn-3+j)=funaVI1(1:nn);
 end

vbd_fI=zeros(nn,nb);
for iline1=1:nn,
%     betag1=betag(iline1,:);
%     dbetag1=k*betag1';
%     dbetag1(1)=dbetag1(1)+2*pi;
%     dbetag1(na)=dbetag1(na)+2*pi;
%     betag1d=(p\dbetag1);
    funb1=vb_fI(iline1,:);
    funbd1=p\(k*funb1');
    funbd2=funbd1;%./betag1d;
    vbd_fI(iline1,:)=funbd2; 
end

dg_alfa=zeros(nn,nn,6);dg_beta=zeros(nn,nn,6);

for i=1:nn,
    for j=1:nn,
        dg_alfa(i,j,1) = vad_fI(i,j);
        dg_beta(i,j,1) = vbd_fI(i,j);
        dg_alfa(i,j,3)=vad_fI(2*(nn-1)+i,nn-j+1);
        dg_beta(i,j,3)=vbd_fI(i,2*(nn-1)+nn-j+1); 
    end
end

div_fI=zeros(nn,nn);
gt_I=gtIa_I;
for i=1:nn,
    for j=1:nn,
%         xwk=(1+xx(i,j)^2)*(1+yy(i,j)^2)/(deltab(i,j)^2);
%         xwk_xi  = (1/sqrt(1+yy(i,j)^2))*xwk;
%         xwk_eta = (1/sqrt(1+xx(i,j)^2))*xwk;
%         div_fI(i,j)=(xwk_xi*dg_alfa(i,j,1) + xwk_eta*dg_beta(i,j,1))/gt_I(i,j);
        div_fI(i,j)=(dg_alfa(i,j,1) + dg_beta(i,j,1))/gt_I(i,j);
    end
end

div_fIII=zeros(nn,nn);

gt_III=gtIa_III; 

for i=1:nn,
    for j=1:nn,
%         xwk=(1+xx(i,j)^2)*(1+yy(i,j)^2)/(deltab(i,j)^2);
%         xwk_xi  = (1/sqrt(1+yy(i,j)^2))*xwk;
%         xwk_eta = (1/sqrt(1+xx(i,j)^2))*xwk;
%         div_fIII(i,j)=(xwk_xi*dg_alfa(i,j,3) - xwk_eta*dg_beta(i,j,3))/gt_III(i,j);
        div_fIII(i,j)=(dg_alfa(i,j,3) + dg_beta(i,j,3))/gt_III(i,j);
    end
end

%% ************************************************************************

xxtIIa_II=zeros(nn,nn);xxtIIa_III=zeros(nn,nn);
xxtIIa_IV=zeros(nn,nn);xxtIIa_I=zeros(nn,nn);
yytIIa_II=zeros(nn,nn);yytIIa_III=zeros(nn,nn);
yytIIa_IV=zeros(nn,nn);yytIIa_I=zeros(nn,nn);


xwk1=fun5(y_fIII);
xwk2=fun7(y_fIII);
for i=1:nn,
    for j=1:nn,
        xxtIIa_II(i,j)=-x_fII(i,j)/y_fII(i,j);
        yytIIa_II(i,j)=z_fII(i,j)/y_fII(i,j);
        xxtIIa_III(i,j)=-x_fIII(i,j)*xwk1(i,j);
        yytIIa_III(i,j)=z_fIII(i,j)*xwk2(i,j);
        xxtIIa_IV(i,j)=-x_fIV(i,j)/y_fIV(i,j);
        yytIIa_IV(i,j)=-z_fIV(i,j)/y_fIV(i,j);
    end
end

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

gxiIIa_III=zeros(nn,nn,3);
gxiIIa_I=zeros(nn,nn,3);

gxiIIa_II=gxi_II;

xwk1=fun5(y_fIII);
xwk2=fun7(y_fIII);
for i=1:nn,
    for j=1:nn,
      gxiIIa_III(i,j,1)= -xwk1(i,j);
      gxiIIa_III(i,j,2)= x_fIII(i,j)*xwk1(i,j)*xwk1(i,j);
      gxiIIa_III(i,j,3)= 0;
      gxiIIa_III(i,j,1:3)=gxiIIa_III(i,j,1:3)/(1+(x_fIII(i,j)*xwk2(i,j))^2);
    end
end

gxiIIa_IV=gxi_IV;

xwk1=fun5(y_fI);
xwk2=fun7(y_fI);
for i=1:nn,
    for j=1:nn,
      gxiIIa_I(i,j,1)= -xwk1(i,j);
      gxiIIa_I(i,j,2)= x_fI(i,j)*xwk1(i,j)*xwk1(i,j);
      gxiIIa_I(i,j,3)= 0;
      gxiIIa_I(i,j,1:3)=gxiIIa_I(i,j,1:3)/(1+(x_fI(i,j)*xwk2(i,j))^2);
    end
end

funfI=zeros(nn,nn);funfII=zeros(nn,nn);funfIII=zeros(nn,nn);
funfIV=zeros(nn,nn);funfV=zeros(nn,nn);funfVI=zeros(nn,nn);

for i=1:nn,
    for j=1:nn,
      funfII(i,j)=dot(mfunfII(i,j,:),gxiIIa_II(i,j,:));
      funfII(i,j)=funfII(i,j)*gtIIa_II(i,j);
      funfIII(i,j)=dot(mfunfIII(i,j,:),gxiIIa_III(i,j,:));
      funfIII(i,j)=funfIII(i,j)*gtIIa_III(i,j);
      funfIV(i,j)=dot(mfunfIV(i,j,:),gxiIIa_IV(i,j,:));
      funfIV(i,j)=funfIV(i,j)*gtIIa_IV(i,j);
      funfI(i,j)=dot(mfunfI(i,j,:),gxiIIa_I(i,j,:));
      funfI(i,j)=funfI(i,j)*gtIIa_I(i,j);
     end
end

va_fII=zeros(4*(nn-1),nn); 
funbIII1=zeros(nn,1);
funbI1=zeros(nn,1);

for jline1=1:nn,
    va_fII(1:nn-1,jline1)=funfII(1:nn-1,jline1);
end

for i=1:nn-1 
    betaspline(1:nn)=beta(i,1:nn);
    funspline(1:nn)=funfIII(i,1:nn);
    ppspline=spline(betaspline,funspline);
    funbIII1(1:nn)=ppval(ppspline,betacr(i,1:nn));
    va_fII(nn-1+i,1:nn)=funbIII1(1:nn);
end

for jline1=1:nn,
    va_fII(2*nn-1:3*nn-3,jline1)=funfIV(1:nn-1,nn-jline1+1); 
end

for i=1:nn-1, 
   betaspline(1:nn)=beta(i,1:nn);
   funspline(1:nn)=funfI(i,1:nn);
   ppspline=spline(betaspline,funspline);
   funbI1(1:nn)=ppval(ppspline,betacr(i,nn+1-1:nn));
   va_fII(3*nn-3+i,1:nn)=funbI1(1:nn);
end

vad_fII=zeros(na,nn);
for jline1=1:nn,
%     alfag5=alfag(:,jline1);
%     dalfag5=k*alfag5;
%     dalfag5(1)=dalfag5(1)+2*pi;
%     dalfag5(na)=dalfag5(na)+2*pi;
%     alfag5d=(p\dalfag5);
    funa9=va_fII(:,jline1);
    funad9=p\(k*funa9);
    funad10=funad9;%./alfag5d;
    vad_fII(:,jline1)=funad10; 
end

%% ************************************************************************

xwk=fun7(y_fV);
xxtIIb_II=-x_fII./y_fII;
yytIIb_II=z_fII./y_fII;
xxtIIb_V=-x_fV.*xwk;
yytIIb_V=z_fV.*xwk;
xxtIIb_IV=-x_fIV./y_fIV;
yytIIb_IV=-z_fIV./y_fIV;
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

getaIIb_V=zeros(nn,nn,3);
getaIIb_VI=zeros(nn,nn,3);

getaIIb_II=geta_II;

xwk1=fun5(y_fV);
xwk2=fun7(y_fV);
for i=1:nn,
    for j=1:nn,
      getaIIb_V(i,j,1)=0;
      getaIIb_V(i,j,2)=-z_fV(i,j)*xwk1(i,j);
      getaIIb_V(i,j,3)= 1;
      getaIIb_V(i,j,1:3)=xwk2(i,j)*getaIIb_V(i,j,1:3)/(1+(z_fV(i,j)*xwk2(i,j))^2);
    end
end

getaIIb_IV=geta_IV;

xwk1=fun5(y_fVI);
xwk2=fun7(y_fVI);
for i=1:nn,
    for j=1:nn,
      getaIIb_VI(i,j,1)= 0;
      getaIIb_VI(i,j,2)=-z_fVI(i,j)*xwk1(i,j)*xwk2(i,j);
      getaIIb_VI(i,j,3)= xwk2(i,j);
      getaIIb_VI(i,j,1:3)=getaIIb_VI(i,j,1:3)/(1+(z_fVI(i,j)*xwk2(i,j))^2);
    end
end

for i=1:nn,
    for j=1:nn,
      funfII(i,j)=dot(mfunfII(i,j,:),getaIIb_II(i,j,:));
      funfII(i,j)=funfII(i,j)*gtIIb_II(i,j);
      funfV(i,j)=dot(mfunfV(i,j,:),getaIIb_V(i,j,:));
      funfV(i,j)=funfV(i,j)*gtIIb_V(i,j);
      funfIV(i,j)=dot(mfunfIV(i,j,:),getaIIb_IV(i,j,:));
      funfIV(i,j)=funfIV(i,j)*gtIIb_IV(i,j);
      funfVI(i,j)=dot(mfunfVI(i,j,:),getaIIb_VI(i,j,:));
      funfVI(i,j)=funfVI(i,j)*gtIIb_VI(i,j);
     end
end

vb_fII=zeros(nn,4*(nn-1));
funbV1=zeros(nn,1);
funbVI1=zeros(nn,1);

for iline1=1:nn,
    vb_fII(iline1,1:nn-1)=funfII(iline1,1:nn-1);
end

for i=1:nn-1 
    betaspline(1:nn)=beta(nn-i+1,1:nn);
    funspline(1:nn)=funfV(nn-i+1,1:nn);
    ppspline=spline(betaspline,funspline);
    funbV1(1:nn)=ppval(ppspline,betacr(i,1:nn));
    vb_fII(1:nn,nn-1+i)=funbV1(1:nn);
end

for iline1=1:nn,
    vb_fII(iline1,2*nn-1:3*nn-3)=funfIV(iline1,nn:-1:2); 
end

for i=1:nn-1 
    betaspline(1:nn)=beta(i,1:nn);
    funspline(1:nn)=funfVI(i,1:nn);
    ppspline=spline(betaspline,funspline);
    funbVI1(1:nn)=ppval(ppspline,betacr(i,1:nn));
    vb_fII(1:nn,3*nn-3+i)=funbVI1(1:nn);
end

vbd_fII=zeros(nn,nb);
for iline1=1:nn,
%     betag2=betag(iline1,:);
%     dbetag2=k*betag2';
%     dbetag2(1)=dbetag2(1)+2*pi;
%     dbetag2(na)=dbetag2(na)+2*pi;
%     betag2d=(p\dbetag2);
    funb2=vb_fII(iline1,:);
    funbd2=p\(k*funb2');
    funbd3=funbd2;%./betag2d;
    vbd_fII(iline1,:)=funbd3; 
end

for i=1:nn,
    for j=1:nn,
        dg_alfa(i,j,2)=vad_fII(i,j);
        dg_beta(i,j,2)=vbd_fII(i,j);
        dg_alfa(i,j,4)=vad_fII(2*(nn-1)+i,nn-j+1);
        dg_beta(i,j,4)=vbd_fII(i,2*(nn-1)+nn-j+1);
    end
end

div_fII=zeros(nn,nn);
gt_II=gtIIa_II; 
for i=1:nn,
    for j=1:nn,
%         xwk=(1+xx(i,j)^2)*(1+yy(i,j)^2)/(deltab(i,j)^2);
%         xwk_xi  = (1/sqrt(1+yy(i,j)^2))*xwk;
%         xwk_eta = (1/sqrt(1+xx(i,j)^2))*xwk;
%         div_fII(i,j)=(xwk_xi*dg_alfa(i,j,2) + xwk_eta*dg_beta(i,j,2))/gt_II(i,j);
        div_fII(i,j)=(dg_alfa(i,j,2) + dg_beta(i,j,2))/gt_II(i,j);
    end
end

div_fIV=zeros(nn,nn);
gt_IV=gtIIa_IV;
for i=1:nn,
    for j=1:nn,
%         xwk=(1+xx(i,j)^2)*(1+yy(i,j)^2)/(deltab(i,j)^2);
%         xwk_xi  = (1/sqrt(1+yy(i,j)^2))*xwk;
%         xwk_eta = (1/sqrt(1+xx(i,j)^2))*xwk;
%         div_fIV(i,j)=(xwk_xi*dg_alfa(i,j,4) - xwk_eta*dg_beta(i,j,4))/gt_IV(i,j);
        div_fIV(i,j)=(dg_alfa(i,j,4) + dg_beta(i,j,4))/gt_IV(i,j);
    end
end

%% ************************************************************************

xxtVa_V=y_fV./z_fV;
yytVa_V=-x_fV./z_fV;
xwk1=fun5(z_fII);
xwk2=fun7(z_fII);
xxtVa_II=y_fII.*xwk2;
yytVa_II=-x_fII.*xwk1;
xxtVa_VI=-y_fVI./z_fVI;
yytVa_VI=-x_fVI./z_fVI;
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

gxiVa_II=zeros(nn,nn,3);
gxiVa_IV=zeros(nn,nn,3);

gxiVa_V=gxi_V;

xwk1=fun5(z_fII);
xwk2=fun7(z_fII);
for i=1:nn,
    for j=1:nn,
      gxiVa_II(i,j,1)= 0;
      gxiVa_II(i,j,2)= xwk2(i,j);
      gxiVa_II(i,j,3)= -y_fII(i,j)*xwk1(i,j)*xwk2(i,j);
      gxiVa_II(i,j,1:3)=gxiVa_II(i,j,1:3)/(1+(y_fII(i,j)*xwk2(i,j))^2);
    end
end

gxiVa_VI=gxi_VI;

xwk1=fun5(z_fIV);
xwk2=fun7(z_fIV);
for i=1:nn,
    for j=1:nn,
      gxiVa_IV(i,j,1)= 0;
      gxiVa_IV(i,j,2)= xwk2(i,j);
      gxiVa_IV(i,j,3)= -y_fIV(i,j)*xwk1(i,j)*xwk2(i,j);
      gxiVa_IV(i,j,1:3)=gxiVa_IV(i,j,1:3)/(1+(y_fIV(i,j)*xwk2(i,j))^2);
    end
end

funfI=zeros(nn,nn);funfII=zeros(nn,nn);funfIII=zeros(nn,nn);
funfIV=zeros(nn,nn);funfV=zeros(nn,nn);funfVI=zeros(nn,nn);

for i=1:nn,
    for j=1:nn,
      funfV(i,j)=dot(mfunfV(i,j,:),gxiVa_V(i,j,:));
      funfV(i,j)=funfV(i,j)*gtVa_V(i,j);
      funfII(i,j)=dot(mfunfII(i,j,:),gxiVa_II(i,j,:));
      funfII(i,j)=funfII(i,j)*gtVa_II(i,j);
      funfVI(i,j)=dot(mfunfVI(i,j,:),gxiVa_VI(i,j,:));
      funfVI(i,j)=funfVI(i,j)*gtVa_VI(i,j);
      funfIV(i,j)=dot(mfunfIV(i,j,:),gxiVa_IV(i,j,:));
      funfIV(i,j)=funfIV(i,j)*gtVa_IV(i,j);
     end
end

va_fV=zeros(4*(nn-1),nn); 
funaII1=zeros(nn,1);
funaIV1=zeros(nn,1);
for jline1=1:nn, 
    va_fV(1:nn-1,jline1)=funfV(1:nn-1,jline1);
end

for j=1:nn-1, 
    alfaspline(1:nn)=alfa(1:nn,nn+1-j);
    funspline(1:nn)=funfII(1:nn,nn+1-j);
    ppspline=spline(alfaspline,funspline);
    funaII1(1:nn)=ppval(ppspline,alfa1(1:nn,j)); 
    va_fV(nn-1+j,1:nn)=funaII1(1:nn);
end 

for jline1=1:nn,
    va_fV(2*nn-1:3*nn-3,jline1)=funfVI(nn:-1:2,jline1); 
end

for j=1:nn-1, 
    alfaspline(1:nn)=alfa(1:nn,j);
    funspline(1:nn)=funfIV(1:nn,j);
    ppspline=spline(alfaspline,funspline);
    funaIV1(1:nn)=ppval(ppspline,alfa1(1:nn,j));
    va_fV(3*nn-3+j,1:nn)=funaIV1(1:nn);
end 

vad_fV=zeros(na,nn);
for jline1=1:nn,
%     alfag6=alfag(1:na,jline1);
%     dalfag6=k*alfag6;
%     dalfag6(1)=dalfag6(1)+2*pi;
%     dalfag6(na)=dalfag6(na)+2*pi;
%     alfag6d=(p\dalfag6);
    funa11=va_fV(:,jline1);
    funad11=p\(k*funa11);
    funad12=funad11;%./alfag6d;
    vad_fV(:,jline1)=funad12; 
end

%% ************************************************************************

xxtVb_V=y_fV./z_fV;
yytVb_V=-x_fV./z_fV;
xwk1=fun5(z_fIII);
xwk2=fun7(z_fIII);
xxtVb_III=y_fIII.*xwk2;
yytVb_III=-x_fIII.*xwk1;
xxtVb_VI=-y_fVI./z_fVI;
yytVb_VI=-x_fVI./z_fVI;
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

getaVb_III=zeros(nn,nn,3);
getaVb_I=zeros(nn,nn,3);

getaVb_V=geta_V;

xwk1=fun5(z_fIII);
xwk2=fun7(z_fIII);
for i=1:nn,
    for j=1:nn,
      getaVb_III(i,j,1)= -xwk1(i,j); 
      getaVb_III(i,j,2)= 0;
      getaVb_III(i,j,3)= x_fIII(i,j)*xwk1(i,j)*xwk1(i,j);
      getaVb_III(i,j,1:3)=getaVb_III(i,j,1:3)/(1+(x_fIII(i,j)*xwk2(i,j))^2);
    end
end

getaVb_VI=geta_VI;

xwk1=fun5(z_fI);
xwk2=fun7(z_fI);
for i=1:nn,
    for j=1:nn,
      getaVb_I(i,j,1)= -xwk1(i,j); 
      getaVb_I(i,j,2)= 0;
      getaVb_I(i,j,3)= x_fI(i,j)*xwk1(i,j)*xwk1(i,j);
      getaVb_I(i,j,1:3)=getaVb_I(i,j,1:3)/(1+(x_fI(i,j)*xwk2(i,j))^2);
      funfV(i,j)=dot(mfunfV(i,j,:),getaVb_V(i,j,:));
      funfV(i,j)=funfV(i,j)*gtVb_V(i,j);
      funfIII(i,j)=dot(mfunfIII(i,j,:),getaVb_III(i,j,:));
      funfIII(i,j)=funfIII(i,j)*gtVb_III(i,j);
      funfVI(i,j)=dot(mfunfVI(i,j,:),getaVb_VI(i,j,:));
      funfVI(i,j)=funfVI(i,j)*gtVb_VI(i,j);
      funfI(i,j)=dot(mfunfI(i,j,:),getaVb_I(i,j,:));
      funfI(i,j)=funfI(i,j)*gtVb_I(i,j);
     end
end

vb_fV=zeros(nn,4*(nn-1)); 
funaIII1=zeros(nn,1);
funaI1=zeros(nn,1);

for iline1=1:nn, 
    vb_fV(iline1,1:nn-1)=funfV(iline1,1:nn-1);
end

for j=1:nn-1, 
    alfaspline(1:nn)=alfa(1:nn,nn+1-j);
    funspline(1:nn)=funfIII(1:nn,nn+1-j);
    ppspline=spline(alfaspline,funspline);
    funaIII1(1:nn)=ppval(ppspline,alfa1(nn+1-1:nn,j)); 
    vb_fV(1:nn,nn-1+j)=funaIII1(1:nn);
end 

for iline1=1:nn,
    vb_fV(iline1,2*nn-1:3*nn-3)=funfVI(nn-iline1+1,1:nn-1); %
end
for j=1:nn-1,
    alfaspline(1:nn)=alfa(1:nn,j);
    funspline(1:nn)=funfI(1:nn,j);
    ppspline=spline(alfaspline,funspline);
    funaI1(1:nn)=ppval(ppspline,alfa1(nn+1-1:nn,j)); 
    vb_fV(1:nn,3*nn-3+j)=funaI1(1:nn);
end 

vbd_fV=zeros(nn,nb);
for iline1=1:nn,
%     betag5=betag(iline1,:);
%     dbetag5=k*betag5';
%     dbetag5(1)=dbetag5(1)+2*pi;
%     dbetag5(na)=dbetag5(na)+2*pi;
%     betag5d=(p\dbetag5);
    funb5=vb_fV(iline1,:);
    funbd5=p\(k*funb5');
    funbd6=funbd5;%./betag5d;
    vbd_fV(iline1,:)=funbd6; 
end


for i=1:nn,
    for j=1:nn,
        dg_alfa(i,j,5)=vad_fV(i,j);
        dg_beta(i,j,5)=vbd_fV(i,j);
        dg_alfa(i,j,6)=vad_fV(2*(nn-1)+nn-i+1,j);
        dg_beta(i,j,6)=vbd_fV(nn-i+1,2*(nn-1)+j);
    end
end

div_fV=zeros(nn,nn);

%gt_V=gtVb_V; 
for i=1:nn,
    for j=1:nn,
%         xwk=(1+xx(i,j)^2)*(1+yy(i,j)^2)/(deltab(i,j)^2);
%         xwk_xi  = (1/sqrt(1+yy(i,j)^2))*xwk;
%         xwk_eta = (1/sqrt(1+xx(i,j)^2))*xwk;
%        div_fV(i,j)=(xwk_xi*dg_alfa(i,j,5) + xwk_eta*dg_beta(i,j,5))/gt_V(i,j);
        div_fV(i,j)=dg_alfa(i,j,5) + dg_beta(i,j,5);
    end
end

div_fVI=zeros(nn,nn);

gt_VI=gtVa_VI;
for i=1:nn,
    for j=1:nn,
%         xwk=(1+xx(i,j)^2)*(1+yy(i,j)^2)/(deltab(i,j)^2);
%         xwk_xi  = (1/sqrt(1+yy(i,j)^2))*xwk;
%         xwk_eta = (1/sqrt(1+xx(i,j)^2))*xwk;
%        div_fVI(i,j)=(-xwk_xi*dg_alfa(i,j,6) + xwk_eta*dg_beta(i,j,6))/gt_VI(i,j);
        div_fVI(i,j)=(dg_alfa(i,j,6) + dg_beta(i,j,6))/gt_VI(i,j);
    end
end

uwk_I=div_fI(1:nn,1:nn,1);uwk_II=div_fII(1:nn,1:nn,1);uwk_III=div_fIII(1:nn,1:nn,1);
uwk_IV=div_fIV(1:nn,1:nn,1);uwk_V=div_fV(1:nn,1:nn,1);uwk_VI=div_fVI(1:nn,1:nn,1);
[div_fI(1:nn,1:nn),div_fII(1:nn,1:nn),div_fIII(1:nn,1:nn),div_fIV(1:nn,1:nn),div_fV(1:nn,1:nn),div_fVI(1:nn,1:nn,1)]=...
    ds72(uwk_I,uwk_II,uwk_III,uwk_IV,uwk_V,uwk_VI,n,nn);
