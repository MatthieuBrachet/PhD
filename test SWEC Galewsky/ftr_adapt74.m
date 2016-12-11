function [funftI,funftII,funftIII,funftIV,funftV,funftVI]=...
    ftr_adapt74(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn)
global na nb;
global alfa beta;
global  betacr;
global alfa1;
global B_adap ftr1 ftr;
global dxi


epsilon=10^-16;
rth=10^-5;

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
funbII1=zeros(nn,1);
funbIV1=zeros(nn,1);
% Face I : 
for jline1=1:nn, 
    va_fI(1:nn-1,jline1)=funftI(1:nn-1,jline1);
end
% FACE II: 
for i=1:nn-1 
    betaspline(1:nn)=beta(i,1:nn);
    funspline(1:nn)=funftII(i,1:nn);
    ppspline=spline(betaspline,funspline);
    funbII1(1:nn)=ppval(ppspline,betacr(i,1:nn));
    va_fI(nn-1+i,1:nn)=funbII1(1:nn);
end
% FACE III:
for jline1=1:nn,
    va_fI(2*nn-1:3*nn-3,jline1)=funftIII(1:nn-1,nn-jline1+1); 
end
% FACE IV: 
for i=1:nn-1, 
   betaspline(1:nn)=beta(i,1:nn);
   funspline(1:nn)=funftIV(i,1:nn);
   ppspline=spline(betaspline,funspline);
   funbIV1(1:nn)=ppval(ppspline,betacr(i,nn+1-[1:nn]));
   va_fI(3*nn-3+i,1:nn)=funbIV1(1:nn);
end
% FILTRAGE
for jline1=1:nn,
    u=va_fI(1:na,jline1);
    uph=u-ftr*u;
    r=0.5*((B_adap*uph).^2+((B_adap')*uph).^2)./dxi^2+epsilon;
    sigma=0.5*(1-rth./r+abs(1-rth./r));
    va_fI(1:na,jline1)=sigma.*(ftr1*u)+(1-sigma).*u;
end
% TRANSFERT FACE I+III
for jline1=1:nn,
   funftI(1:nn,jline1)= va_fI(1:nn,jline1);
   funftIII(1:nn,nn-jline1+1)= va_fI(2*nn-1:3*nn-2,jline1);
end

%% RESEAU 2: 

vb_fI=zeros(nn,4*(nn-1)); 
funaV1=zeros(nn,1);
funaVI1=zeros(nn,1);
% Face I :
for iline1=1:nn,
    vb_fI(iline1,1:nn-1)=funftI(iline1,1:nn-1);
end
% FACE V: 
for j=1:nn-1
    alfaspline(1:nn)=alfa(1:nn,j);
    funspline(1:nn)=funftV(1:nn,j);
    ppspline=spline(alfaspline,funspline);
    funaV1(1:nn)=ppval(ppspline,alfa1(1:nn,j));
    vb_fI(1:nn,nn-1+j)=funaV1(1:nn);
end
% FACE III:
 for iline1=1:nn,
     vb_fI(iline1,2*nn-1:3*nn-3)=funftIII(iline1,nn:-1:2); 
 end
% FACE VI:
 for j=1:nn-1, 
    alfaspline(1:nn)=alfa(1:nn,j);
    funspline(1:nn)=funftVI(1:nn,j);
    ppspline=spline(alfaspline,funspline);
    funaVI1(1:nn)=ppval(ppspline,alfa1(nn+1-[1:nn],j));
    vb_fI(1:nn,3*nn-3+j)=funaVI1(1:nn);
 end
 % FILTRAGE
for iline1=1:nn,
    u=vb_fI(iline1,1:nb)';
    uph=u-ftr*u;
    r=0.5*((B_adap*uph).^2+((B_adap')*uph).^2)./dxi^2+epsilon;
    sigma=0.5*(1-rth./r+abs(1-rth./r));
    vb_fI(iline1,1:nb)=sigma.*(ftr1*u)+(1-sigma).*u;
end
% TRANSFERT FACE I+III
for iline1=1:nn,
   funftI(iline1,1:nn) = vb_fI(iline1,1:nn);
   funftIII(iline1,nn:-1:1) = vb_fI(iline1,2*nn-1:3*nn-2);
end

%% RESEAU 3: 
va_fII=zeros(4*(nn-1),nn);
funbIII1=zeros(nn,1);
funbI1=zeros(nn,1);
% Face II :
for jline1=1:nn, 
    va_fII(1:nn-1,jline1)=funftII(1:nn-1,jline1);
end
% FACE III: 
for i=1:nn-1 
    betaspline(1:nn)=beta(i,1:nn);
    funspline(1:nn)=funftIII(i,1:nn);
    ppspline=spline(betaspline,funspline);
    funbIII1(1:nn)=ppval(ppspline,betacr(i,1:nn));
    va_fII(nn-1+i,1:nn)=funbIII1(1:nn);
end
% FACE IV: 
for jline1=1:nn,
    va_fII(2*nn-1:3*nn-3,jline1)=funftIV(1:nn-1,nn-jline1+1);  
end
% FACE I: 
for i=1:nn-1,
   betaspline(1:nn)=beta(i,1:nn);
   funspline(1:nn)=funftI(i,1:nn);
   ppspline=spline(betaspline,funspline);
   funbI1(1:nn)=ppval(ppspline,betacr(i,nn+1-[1:nn]));
   va_fII(3*nn-3+i,1:nn)=funbI1(1:nn);
end
% FILTRAGE
for jline1=1:nn,
    u=va_fII(1:na,jline1);
    uph=u-ftr*u;
    r=0.5*((B_adap*uph).^2+((B_adap')*uph).^2)./dxi^2+epsilon;
    sigma=0.5*(1-rth./r+abs(1-rth./r));
    va_fII(1:na,jline1)=sigma.*(ftr1*u)+(1-sigma).*u;
end
% TRANSFERT FACE II+IV
for jline1=1:nn,
   funftII(1:nn,jline1) = va_fII(1:nn,jline1);
   funftIV(1:nn,nn-jline1+1) = va_fII(2*nn-1:3*nn-2,jline1);
end

%% RESEAU 4:
vb_fII=zeros(nn,4*(nn-1));
funbV1=zeros(nn,1);
funbVI1=zeros(nn,1);
% Face II :
for iline1=1:nn, 
    vb_fII(iline1,1:nn-1)=funftII(iline1,1:nn-1);
end
% FACE V: 
for i=1:nn-1 
    betaspline(1:nn)=beta(nn-i+1,1:nn);
    funspline(1:nn)=funftV(nn-i+1,1:nn);
    ppspline=spline(betaspline,funspline);
    funbV1(1:nn)=ppval(ppspline,betacr(i,1:nn));
    vb_fII(1:nn,nn-1+i)=funbV1(1:nn);
end
% FACE III: 
for iline1=1:nn,
    vb_fII(iline1,2*nn-1:3*nn-3)=funftIV(iline1,nn:-1:2); 
end
% FACE VI: 
for i=1:nn-1 
    betaspline(1:nn)=beta(i,1:nn);
    funspline(1:nn)=funftVI(i,1:nn);
    ppspline=spline(betaspline,funspline);
    funbVI1(1:nn)=ppval(ppspline,betacr(i,1:nn));
    vb_fII(1:nn,3*nn-3+i)=funbVI1(1:nn);
end
 % FILTRAGE
for iline1=1:nn,
    u=vb_fII(iline1,1:nb)';
    uph=u-ftr*u;
    r=0.5*((B_adap*uph).^2+((B_adap')*uph).^2)./dxi^2+epsilon;
    sigma=0.5*(1-rth./r+abs(1-rth./r));
    vb_fII(iline1,1:nb)=sigma.*(ftr1*u)+(1-sigma).*u;
end
% TRANSFERT FACE II+IV
for iline1=1:nn,
   funftII(iline1,1:nn) = vb_fII(iline1,1:nn);
   funftIV(iline1,nn:-1:1) = vb_fII(iline1,2*nn-1:3*nn-2);
end

%% RESEAU 5: 

va_fV=zeros(4*(nn-1),nn); 
funaII1=zeros(nn,1);
funaIV1=zeros(nn,1);
% FACE V
for jline1=1:nn,
    va_fV(1:nn-1,jline1)=funftV(1:nn-1,jline1);
end
% FACE II: 
for j=1:nn-1, 
    alfaspline(1:nn)=alfa(1:nn,nn+1-j);
    funspline(1:nn)=funftII(1:nn,nn+1-j);
    ppspline=spline(alfaspline,funspline);
    funaII1(1:nn)=ppval(ppspline,alfa1(1:nn,j)); 
    va_fV(nn-1+j,1:nn)=funaII1(1:nn);
end 
% FACE VI: 
for jline1=1:nn
    va_fV(2*nn-1:3*nn-3,jline1)=funftVI(nn:-1:2,jline1); 
end
% FACE IV:
for j=1:nn-1,
    alfaspline(1:nn)=alfa(1:nn,j);
    funspline(1:nn)=funftIV(1:nn,j);
    ppspline=spline(alfaspline,funspline);
    funaIV1(1:nn)=ppval(ppspline,alfa1(1:nn,j));
    va_fV(3*nn-3+j,1:nn)=funaIV1(1:nn);
end 
% FILTRAGE
for jline1=1:nn,
    u=va_fV(1:na,jline1);
    uph=u-ftr*u;
    r=0.5*((B_adap*uph).^2+((B_adap')*uph).^2)./dxi^2+epsilon;
    sigma=0.5*(1-rth./r+abs(1-rth./r));
    va_fV(1:na,jline1)=sigma.*(ftr1*u)+(1-sigma).*u;
end
% TRANSFERT FACE V+VI
for jline1=1:nn,
   funftV(1:nn,jline1) = va_fV(1:nn,jline1);
   funftVI(nn:-1:1,jline1) = va_fV(2*nn-1:3*nn-2,jline1);
end

%% RESEAU 6: ASSEMBLAGE DES DONNEES SUR LE RESEAU V-BETA
vb_fV=zeros(nn,4*(nn-1)); 
funaIII1=zeros(nn,1);
funaI1=zeros(nn,1);
% Face V : 
for iline1=1:nn, 
    vb_fV(iline1,1:nn-1)=funfV(iline1,1:nn-1);
end
% FACE III
for j=1:nn-1,
    alfaspline(1:nn)=alfa(1:nn,nn+1-j);
    funspline(1:nn)=funfIII(1:nn,nn+1-j);
    ppspline=spline(alfaspline,funspline);
    funaIII1(1:nn)=ppval(ppspline,alfa1(nn+1-[1:nn],j));
    vb_fV(1:nn,nn-1+j)=funaIII1(1:nn);
end 
% FACE IV
for iline1=1:nn,
    vb_fV(iline1,2*nn-1:3*nn-3)=funfVI(nn-iline1+1,1:nn-1); %
end
% FACE I
for j=1:nn-1,
    alfaspline(1:nn)=alfa(1:nn,j);
    funspline(1:nn)=funfI(1:nn,j);
    ppspline=spline(alfaspline,funspline);
    funaI1(1:nn)=ppval(ppspline,alfa1(nn+1-[1:nn],j));
    vb_fV(1:nn,3*nn-3+j)=funaI1(1:nn);
end 
 % FILTRAGE
for iline1=1:nn,
    u=vb_fV(iline1,1:nb)';
    uph=u-ftr*u;
    r=0.5*((B_adap*uph).^2+((B_adap')*uph).^2)./dxi^2+epsilon;
    sigma=0.5*(1-rth./r+abs(1-rth./r));
    vb_fV(iline1,1:nb)=sigma.*(ftr1*u)+(1-sigma).*u;
end
% TRANSFERT FACE V+VI
for iline1=1:nn,
   funftV(iline1,1:nn) = vb_fV(iline1,1:nn);
   funftVI(nn-iline1+1,1:nn) = vb_fV(iline1,2*nn-1:3*nn-2);
end