function [div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
    div72(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI,n,nn)
global gxi_I gxi_II gxi_III gxi_IV gxi_V gxi_VI;
global geta_I geta_II geta_III geta_IV geta_V geta_VI;
%% ************************************************************************
% METHODE 2 POUR CALCUL DE LA DIVERGENCE SPHEREIQUE D'UN CHAMP DE VECTEURS
% TANGENTIEL. vOIR FORMULE (34) DANS PAPIER CS5 (HERMITIAN APPROXIMATION
% OF THE SPHERICAL DIVERGENCE ON THE CUBED SPHERE).
%% *** ETAPE 1: ***********************************************************
% --------
% CALCUL DE D/DXI (F) ET D/DETA (F)AVEC F=CHAMP DE VECTEURS.
% 3 APPELS POUR LES 3 COMPOSANTES DE F
grxivec_I=zeros(nn,nn,3);grxivec_II=zeros(nn,nn,3);grxivec_III=zeros(nn,nn,3);
grxivec_IV=zeros(nn,nn,3);grxivec_V=zeros(nn,nn,3);grxivec_VI=zeros(nn,nn,3);

gretavec_I=zeros(nn,nn,3);gretavec_II=zeros(nn,nn,3);gretavec_III=zeros(nn,nn,3);
gretavec_IV=zeros(nn,nn,3);gretavec_V=zeros(nn,nn,3);gretavec_VI=zeros(nn,nn,3);

%% 1- COMPOSANTE X DU CHAMP DE VECTEURS MFUN:
funfI=mfunfI(:,:,1);
funfII=mfunfII(:,:,1);
funfIII=mfunfIII(:,:,1);
funfIV=mfunfIV(:,:,1);
funfV=mfunfV(:,:,1);
funfVI=mfunfVI(:,:,1);

[grxi_I,grxi_II,grxi_III,grxi_IV,grxi_V,grxi_VI...
    ,greta_I,greta_II,greta_III,greta_IV,greta_V,greta_VI]=...
    grxieta72(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn);

grxivec_I(:,:,1)=grxi_I;
grxivec_II(:,:,1)=grxi_II;
grxivec_III(:,:,1)=grxi_III;
grxivec_IV(:,:,1)=grxi_IV;
grxivec_V(:,:,1)=grxi_V;
grxivec_VI(:,:,1)=grxi_VI;

gretavec_I(:,:,1)=greta_I;
gretavec_II(:,:,1)=greta_II;
gretavec_III(:,:,1)=greta_III;
gretavec_IV(:,:,1)=greta_IV;
gretavec_V(:,:,1)=greta_V;
gretavec_VI(:,:,1)=greta_VI;

%% 2- COMPOSANTE Y DU CHAMP DE VECTEURS MFUN:
funfI=mfunfI(:,:,2);
funfII=mfunfII(:,:,2);
funfIII=mfunfIII(:,:,2);
funfIV=mfunfIV(:,:,2);
funfV=mfunfV(:,:,2);
funfVI=mfunfVI(:,:,2);

[grxi_I,grxi_II,grxi_III,grxi_IV,grxi_V,grxi_VI...
    ,greta_I,greta_II,greta_III,greta_IV,greta_V,greta_VI]=...
    grxieta72(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn);

grxivec_I(:,:,2)=grxi_I;
grxivec_II(:,:,2)=grxi_II;
grxivec_III(:,:,2)=grxi_III;
grxivec_IV(:,:,2)=grxi_IV;
grxivec_V(:,:,2)=grxi_V;
grxivec_VI(:,:,2)=grxi_VI;

gretavec_I(:,:,2)=greta_I;
gretavec_II(:,:,2)=greta_II;
gretavec_III(:,:,2)=greta_III;
gretavec_IV(:,:,2)=greta_IV;
gretavec_V(:,:,2)=greta_V;
gretavec_VI(:,:,2)=greta_VI;

%% 3- COMPOSANTE Z DU CHAMP DE VECTEURS MFUN:
funfI=mfunfI(:,:,3);
funfII=mfunfII(:,:,3);
funfIII=mfunfIII(:,:,3);
funfIV=mfunfIV(:,:,3);
funfV=mfunfV(:,:,3);
funfVI=mfunfVI(:,:,3);

[grxi_I,grxi_II,grxi_III,grxi_IV,grxi_V,grxi_VI...
    ,greta_I,greta_II,greta_III,greta_IV,greta_V,greta_VI]=...
    grxieta72(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn);

grxivec_I(:,:,3)=grxi_I;
grxivec_II(:,:,3)=grxi_II;
grxivec_III(:,:,3)=grxi_III;
grxivec_IV(:,:,3)=grxi_IV;
grxivec_V(:,:,3)=grxi_V;
grxivec_VI(:,:,3)=grxi_VI;

gretavec_I(:,:,3)=greta_I;
gretavec_II(:,:,3)=greta_II;
gretavec_III(:,:,3)=greta_III;
gretavec_IV(:,:,3)=greta_IV;
gretavec_V(:,:,3)=greta_V;
gretavec_VI(:,:,3)=greta_VI;
%% *** ETAPE 2: ***********************************************************
div_fI=zeros(nn,nn);
for i=1:nn,
    for j=1:nn,
        div_fI(i,j)=grxivec_I(i,j,1)*gxi_I(i,j,1)+grxivec_I(i,j,2)*gxi_I(i,j,2)+grxivec_I(i,j,3)*gxi_I(i,j,3)...
            + gretavec_I(i,j,1)*geta_I(i,j,1)+gretavec_I(i,j,2)*geta_I(i,j,2)+gretavec_I(i,j,3)*geta_I(i,j,3);
    end
end
div_fII=zeros(nn,nn);
for i=1:nn,
    for j=1:nn,
        div_fII(i,j)=grxivec_II(i,j,1)*gxi_II(i,j,1)+grxivec_II(i,j,2)*gxi_II(i,j,2)+grxivec_II(i,j,3)*gxi_II(i,j,3)...
            + gretavec_II(i,j,1)*geta_II(i,j,1)+gretavec_II(i,j,2)*geta_II(i,j,2)+gretavec_II(i,j,3)*geta_II(i,j,3);
    end
end
div_fIII=zeros(nn,nn);
for i=1:nn,
    for j=1:nn,
        div_fIII(i,j)=grxivec_III(i,j,1)*gxi_III(i,j,1)+grxivec_III(i,j,2)*gxi_III(i,j,2)+grxivec_III(i,j,3)*gxi_III(i,j,3)...
            + gretavec_III(i,j,1)*geta_III(i,j,1)+gretavec_III(i,j,2)*geta_III(i,j,2)+gretavec_III(i,j,3)*geta_III(i,j,3);
    end
end
div_fIV=zeros(nn,nn);
for i=1:nn,
    for j=1:nn,
        div_fIV(i,j)=grxivec_IV(i,j,1)*gxi_IV(i,j,1)+grxivec_IV(i,j,2)*gxi_IV(i,j,2)+grxivec_IV(i,j,3)*gxi_IV(i,j,3)...
            + gretavec_IV(i,j,1)*geta_IV(i,j,1)+gretavec_IV(i,j,2)*geta_IV(i,j,2)+gretavec_IV(i,j,3)*geta_IV(i,j,3);
    end
end
div_fV=zeros(nn,nn);
for i=1:nn,
    for j=1:nn,
        div_fV(i,j)=grxivec_V(i,j,1)*gxi_V(i,j,1)+grxivec_V(i,j,2)*gxi_V(i,j,2)+grxivec_V(i,j,3)*gxi_V(i,j,3)...
            + gretavec_V(i,j,1)*geta_V(i,j,1)+gretavec_V(i,j,2)*geta_V(i,j,2)+gretavec_V(i,j,3)*geta_V(i,j,3);
    end
end
div_fVI=zeros(nn,nn);
for i=1:nn,
    for j=1:nn,
        div_fVI(i,j)=grxivec_VI(i,j,1)*gxi_VI(i,j,1)+grxivec_VI(i,j,2)*gxi_VI(i,j,2)+grxivec_VI(i,j,3)*gxi_VI(i,j,3)...
            + gretavec_VI(i,j,1)*geta_VI(i,j,1)+gretavec_VI(i,j,2)*geta_VI(i,j,2)+gretavec_VI(i,j,3)*geta_VI(i,j,3);
    end
end
%% ************************************************************************
% % demi-somme + 1/3 somme de la divergence
uwk_I=div_fI(1:nn,1:nn,1);uwk_II=div_fII(1:nn,1:nn,1);uwk_III=div_fIII(1:nn,1:nn,1);
uwk_IV=div_fIV(1:nn,1:nn,1);uwk_V=div_fV(1:nn,1:nn,1);uwk_VI=div_fVI(1:nn,1:nn,1);
[div_fI(1:nn,1:nn),div_fII(1:nn,1:nn),div_fIII(1:nn,1:nn),div_fIV(1:nn,1:nn),div_fV(1:nn,1:nn),div_fVI(1:nn,1:nn,1)]=...
    ds72(uwk_I,uwk_II,uwk_III,uwk_IV,uwk_V,uwk_VI,n,nn);

ee=1; ite=0;
while ee>1e-14 & ite<50
    [divf_fI,divf_fII,divf_fIII,divf_fIV,divf_fV,divf_fVI]=ftr_div72(div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI,n,nn);
    ite=ite+1;
    
    e_fI=abs(divf_fI-div_fI);
    e_fII=abs(divf_fII-div_fII);
    e_fIII=abs(divf_fIII-div_fIII);
    e_fIV=abs(divf_fIV-div_fIV);
    e_fV=abs(divf_fV-div_fV);
    e_fVI=abs(divf_fVI-div_fVI);
    ee=max(max([e_fI, e_fII, e_fIII, e_fIV, e_fV, e_fVI]));
    
    div_fI=divf_fI;
    div_fII=divf_fII;
    div_fIII=divf_fIII;
    div_fIV=divf_fIV;
    div_fV=divf_fV;
    div_fVI=divf_fVI;
end

