function [grxi_I,grxi_II,grxi_III,grxi_IV,grxi_V,grxi_VI...
    ,greta_I,greta_II,greta_III,greta_IV,greta_V,greta_VI]=...
    grxieta72(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn)
global na nb;
global alfa beta;
global betacr;
global alfa1;
global p_div k_div;
grxi_I=zeros(nn,nn); grxi_II=zeros(nn,nn);grxi_III=zeros(nn,nn);
grxi_IV=zeros(nn,nn); grxi_V=zeros(nn,nn);grxi_VI=zeros(nn,nn);
greta_I=zeros(nn,nn); greta_II=zeros(nn,nn);greta_III=zeros(nn,nn);
greta_IV=zeros(nn,nn); greta_V=zeros(nn,nn);greta_VI=zeros(nn,nn);

betaspline=zeros(nn,1);
alfaspline=zeros(nn,1);
funspline=zeros(nn,1);

%% RESEAU 1 : ASSEMBLAGE DES DONNEES SUR LE RESEAU I-ALPHA

va_fI=zeros(4*(nn-1),nn); 
funbII1=zeros(nn,1);
funbIV1=zeros(nn,1);
% Face I : 
for jline1=1:nn, 
    va_fI(1:nn-1,jline1)=funfI(1:nn-1,jline1);
end
% FACE II: 
for i=1:nn-1 
    betaspline(1:nn)=beta(i,1:nn);
    funspline(1:nn)=funfII(i,1:nn);
    ppspline=spline(betaspline,funspline);
    funbII1(1:nn)=ppval(ppspline,betacr(i,1:nn));
    va_fI(nn-1+i,1:nn)=funbII1(1:nn);
end
% FACE III: 
for jline1=1:nn,
    va_fI(2*nn-1:3*nn-3,jline1)=funfIII(1:nn-1,nn-jline1+1); 
end
% FACE IV: 
for i=1:nn-1, 
   betaspline(1:nn)=beta(i,1:nn);
   funspline(1:nn)=funfIV(i,1:nn);
   ppspline=spline(betaspline,funspline);
   funbIV1(1:nn)=ppval(ppspline,betacr(i,nn+1-[1:nn]));
   va_fI(3*nn-3+i,1:nn)=funbIV1(1:nn);
end

vad_fI=zeros(na,nn);
for jline1=1:nn, 
    funa7=va_fI(:,jline1); 
    funad7=p_div\(k_div*funa7);
    vad_fI(:,jline1)=funad7;
end

%% RESEAU 2: ASSEMBLAGE DES DONNEES SUR LE RESEAU I-BETA

vb_fI=zeros(nn,4*(nn-1)); 
funaV1=zeros(nn,1);
funaVI1=zeros(nn,1);
% Face I : 
for iline1=1:nn, 
    vb_fI(iline1,1:nn-1)=funfI(iline1,1:nn-1);
end
% FACE V: 
for j=1:nn-1 
    alfaspline(1:nn)=alfa(1:nn,j);
    funspline(1:nn)=funfV(1:nn,j);
    ppspline=spline(alfaspline,funspline);
    funaV1(1:nn)=ppval(ppspline,alfa1(1:nn,j));
    vb_fI(1:nn,nn-1+j)=funaV1(1:nn);
end
% FACE III:
 for iline1=1:nn,
     vb_fI(iline1,2*nn-1:3*nn-3)=funfIII(iline1,nn:-1:2); 
 end
% FACE VI: 
 for j=1:nn-1, 
    alfaspline(1:nn)=alfa(1:nn,j);
    funspline(1:nn)=funfVI(1:nn,j);
    ppspline=spline(alfaspline,funspline);
    funaVI1(1:nn)=ppval(ppspline,alfa1(nn+1-[1:nn],j));
    vb_fI(1:nn,3*nn-3+j)=funaVI1(1:nn);
 end

vbd_fI=zeros(nn,nb);
for iline1=1:nn,
    funb1=vb_fI(iline1,:);
    funbd1=p_div\(k_div*funb1');
    vbd_fI(iline1,:)=funbd1; 
end

%% RESEAU 3: ASSEMBLAGE DES DONNEES SUR LE RESEAU II-ALPHA
va_fII=zeros(4*(nn-1),nn); 
funbIII1=zeros(nn,1);
funbI1=zeros(nn,1);
% Face II : 
for jline1=1:nn,
    va_fII(1:nn-1,jline1)=funfII(1:nn-1,jline1);
end
% FACE III: 
for i=1:nn-1 
    betaspline(1:nn)=beta(i,1:nn);
    funspline(1:nn)=funfIII(i,1:nn);
    ppspline=spline(betaspline,funspline);
    funbIII1(1:nn)=ppval(ppspline,betacr(i,1:nn));
    va_fII(nn-1+i,1:nn)=funbIII1(1:nn);
end
% FACE IV:
for jline1=1:nn,
    va_fII(2*nn-1:3*nn-3,jline1)=funfIV(1:nn-1,nn-jline1+1); 
end
% FACE I:
for i=1:nn-1, 
   betaspline(1:nn)=beta(i,1:nn);
   funspline(1:nn)=funfI(i,1:nn);
   ppspline=spline(betaspline,funspline);
   funbI1(1:nn)=ppval(ppspline,betacr(i,nn+1-[1:nn]));
   va_fII(3*nn-3+i,1:nn)=funbI1(1:nn);
end

vad_fII=zeros(na,nn);
for jline1=1:nn,
    funa9=va_fII(:,jline1);
    funad9=p_div\(k_div*funa9);
    vad_fII(:,jline1)=funad9; 
end

%% RESEAU 4: ASSEMBLAGE DES DONNEES SUR LE RESEAU II-BETA

vb_fII=zeros(nn,4*(nn-1)); 
funbV1=zeros(nn,1);
funbVI1=zeros(nn,1);
% Face II :
for iline1=1:nn,
    vb_fII(iline1,1:nn-1)=funfII(iline1,1:nn-1);
end
% FACE V: 
for i=1:nn-1 
    betaspline(1:nn)=beta(nn-i+1,1:nn);
    funspline(1:nn)=funfV(nn-i+1,1:nn);
    ppspline=spline(betaspline,funspline);
    funbV1(1:nn)=ppval(ppspline,betacr(i,1:nn));
    vb_fII(1:nn,nn-1+i)=funbV1(1:nn);
end
% FACE III: 
for iline1=1:nn,
    vb_fII(iline1,2*nn-1:3*nn-3)=funfIV(iline1,nn:-1:2);
end
% FACE VI: 
for i=1:nn-1 
    betaspline(1:nn)=beta(i,1:nn);
    funspline(1:nn)=funfVI(i,1:nn);
    ppspline=spline(betaspline,funspline);
    funbVI1(1:nn)=ppval(ppspline,betacr(i,1:nn));
    vb_fII(1:nn,3*nn-3+i)=funbVI1(1:nn);
end

vbd_fII=zeros(nn,nb);
for iline1=1:nn,
    funb2=vb_fII(iline1,:);
    funbd2=p_div\(k_div*funb2');
    vbd_fII(iline1,:)=funbd2; 
end

%% RESEAU 5: 

va_fV=zeros(4*(nn-1),nn); 
funaII1=zeros(nn,1);
funaIV1=zeros(nn,1);
% Face V
for jline1=1:nn,
    va_fV(1:nn-1,jline1)=funfV(1:nn-1,jline1);
end
% FACE II: 
for j=1:nn-1,
    alfaspline(1:nn)=alfa(1:nn,nn+1-j);
    funspline(1:nn)=funfII(1:nn,nn+1-j);
    ppspline=spline(alfaspline,funspline);
    funaII1(1:nn)=ppval(ppspline,alfa1(1:nn,j)); 
    va_fV(nn-1+j,1:nn)=funaII1(1:nn);
end 
% FACE VI: 
for jline1=1:nn,
    va_fV(2*nn-1:3*nn-3,jline1)=funfVI(nn:-1:2,jline1); 
end
% FACE IV:
for j=1:nn-1,
    alfaspline(1:nn)=alfa(1:nn,j);
    funspline(1:nn)=funfIV(1:nn,j);
    ppspline=spline(alfaspline,funspline);
    funaIV1(1:nn)=ppval(ppspline,alfa1(1:nn,j)); 
    va_fV(3*nn-3+j,1:nn)=funaIV1(1:nn);
end 

vad_fV=zeros(na,nn);
for jline1=1:nn,
    funa11=va_fV(:,jline1);
    funad11=p_div\(k_div*funa11);
    vad_fV(:,jline1)=funad11; 
end

%% RESEAU 6:

vb_fV=zeros(nn,4*(nn-1)); 
funaIII1=zeros(nn,1);
funaI1=zeros(nn,1);
% Face V : 
for iline1=1:nn, 
    vb_fV(iline1,1:nn-1)=funfV(iline1,1:nn-1);
end
% Face III
for j=1:nn-1, 
    alfaspline(1:nn)=alfa(1:nn,nn+1-j);
    funspline(1:nn)=funfIII(1:nn,nn+1-j);
    ppspline=spline(alfaspline,funspline);
    funaIII1(1:nn)=ppval(ppspline,alfa1(nn+1-[1:nn],j));
    vb_fV(1:nn,nn-1+j)=funaIII1(1:nn);
end 
% Face VI
for iline1=1:nn,
    vb_fV(iline1,2*nn-1:3*nn-3)=funfVI(nn-iline1+1,1:nn-1); %
end
% Face I
for j=1:nn-1, 
    alfaspline(1:nn)=alfa(1:nn,j);
    funspline(1:nn)=funfI(1:nn,j);
    ppspline=spline(alfaspline,funspline);
    funaI1(1:nn)=ppval(ppspline,alfa1(nn+1-[1:nn],j));
    vb_fV(1:nn,3*nn-3+j)=funaI1(1:nn);
end 

vbd_fV=zeros(nn,nb);
for iline1=1:nn,
    funb5=vb_fV(iline1,:);
    funbd5=p_div\(k_div*funb5');
    vbd_fV(iline1,:)=funbd5;
end

%% *** Assemblage *********************************************************
dg_alfa=zeros(nn,nn,6);dg_beta=zeros(nn,nn,6);
% FACE I
for i=1:nn,
    for j=1:nn,
        dg_alfa(i,j,1) = vad_fI(i,j);
        dg_beta(i,j,1) = vbd_fI(i,j);
    end
end
% FACE II
for i=1:nn,
    for j=1:nn,
        dg_alfa(i,j,2)=vad_fII(i,j);
        dg_beta(i,j,2)=vbd_fII(i,j);
    end
end
% FACE III
for i=1:nn,
    for j=1:nn,
        dg_alfa(i,j,3)=vad_fI(2*(nn-1)+i,nn-j+1);
        dg_beta(i,j,3)=vbd_fI(i,2*(nn-1)+nn-j+1); 
    end
end
% FACE IV
for i=1:nn,
    for j=1:nn,
        dg_alfa(i,j,4)=vad_fII(2*(nn-1)+i,nn-j+1);
        dg_beta(i,j,4)=vbd_fII(i,2*(nn-1)+nn-j+1);
    end
end
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
%% ETAPE 3 SUIVANTE DIFFERENTE DU CALCUL DU GRADIENT:

% FACE I
for i=1:nn,
    for j=1:nn,
        grxi_I(i,j)=dg_alfa(i,j,1);
        greta_I(i,j)=dg_beta(i,j,1);
    end
end

% Face II
for i=1:nn,
    for j=1:nn,
        grxi_II(i,j)=dg_alfa(i,j,2);
        greta_II(i,j)=dg_beta(i,j,2);

    end
end
% FACE III: 
for i=1:nn,
    for j=1:nn,
        grxi_III(i,j)=dg_alfa(i,j,3);
        greta_III(i,j)=-dg_beta(i,j,3);
    end
end
% Face IV
for i=1:nn,
    for j=1:nn,
        grxi_IV(i,j)=dg_alfa(i,j,4);
        greta_IV(i,j)=-dg_beta(i,j,4);
    end
end
% FACE V:
for i=1:nn,
    for j=1:nn,
        grxi_V(i,j)=dg_alfa(i,j,5);
        greta_V(i,j)=dg_beta(i,j,5);
    end
end
% FACE VI: 
for i=1:nn,
    for j=1:nn,
        grxi_VI(i,j)=-dg_alfa(i,j,6);
        greta_VI(i,j)=dg_beta(i,j,6);
    end
end