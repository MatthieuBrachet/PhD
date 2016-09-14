function [rot_fI,rot_fII,rot_fIII,rot_fIV,rot_fV,rot_fVI]=...
    rot75(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI,n,nn)
% COMPUTATION SPHERICAL CURL (ON THE CUBED SPHERE GRID)

global gxi_I gxi_II gxi_III gxi_IV gxi_V gxi_VI;
global geta_I geta_II geta_III geta_IV geta_V geta_VI;


rot_fI=zeros(nn,nn,3);
rot_fII=zeros(nn,nn,3);
rot_fIII=zeros(nn,nn,3);
rot_fIV=zeros(nn,nn,3);
rot_fV=zeros(nn,nn,3);
rot_fVI=zeros(nn,nn,3);

for kk=1:3
    funfI=mfunfI(:,:,kk);
    funfII=mfunfII(:,:,kk);
    funfIII=mfunfIII(:,:,kk);
    funfIV=mfunfIV(:,:,kk);
    funfV=mfunfV(:,:,kk);
    funfVI=mfunfVI(:,:,kk);
    
    
    [grxi_I,grxi_II,grxi_III,grxi_IV,grxi_V,grxi_VI...
        ,greta_I,greta_II,greta_III,greta_IV,greta_V,greta_VI]=...
        grxieta75(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn);
    
    dg_xi_fI(:,:,kk)=grxi_I;
    dg_xi_fII(:,:,kk)=grxi_II;
    dg_xi_fIII(:,:,kk)=grxi_III;
    dg_xi_fIV(:,:,kk)=grxi_IV;
    dg_xi_fV(:,:,kk)=grxi_V;
    dg_xi_fVI(:,:,kk)=grxi_VI;
    
    dg_eta_fI(:,:,kk)=greta_I;
    dg_eta_fII(:,:,kk)=greta_II;
    dg_eta_fIII(:,:,kk)=greta_III;
    dg_eta_fIV(:,:,kk)=greta_IV;
    dg_eta_fV(:,:,kk)=greta_V;
    dg_eta_fVI(:,:,kk)=greta_VI;
    
end


%% ETAPE 3 - ASSEMBLAGE ROTATIONNEL EN FONCTION DES DERIVEES XI/ETA

% 8.1 - FACE I: CALCUL DU ROT. EN FONCTION DES DERIVEES ALPHA ET BETA
for i=1:nn,
    for j=1:nn,
        rot_fI(i,j,1:3)=cross(gxi_I(i,j,1:3),dg_xi_fI(i,j,1:3)) + cross(geta_I(i,j,1:3), dg_eta_fI(i,j,1:3));
    end
end

% 8.2  FACE II: CALCUL DU ROT. EN FONCTION DES DERIVEES ALPHA ET BETA FACE II
for i=1:nn,
    for j=1:nn,
        rot_fII(i,j,1:3)=cross(gxi_II(i,j,1:3),dg_xi_fII(i,j,1:3)) + cross(geta_II(i,j,1:3),dg_eta_fII(i,j,1:3));
    end
end

% 8.3 -  FACE III: CALCUL DU ROT. EN FONCTION DES DERIVEES ALPHA ET BETA FACE III
for i=1:nn
    for j=1:nn
        rot_fIII(i,j,1:3)=cross(gxi_III(i,j,1:3),dg_xi_fIII(i,j,1:3)) + cross(geta_III(i,j,1:3), dg_eta_fIII(i,j,1:3));
    end
end

% 8.4 - FACE IV: CALCUL DU ROT. EN FONCTION DES DERIVEES ALPHA ET BETA FACE IV
for i=1:nn,
    for j=1:nn,
        rot_fIV(i,j,1:3)=cross(gxi_IV(i,j,1:3),dg_xi_fIV(i,j,1:3)) + cross(geta_IV(i,j,1:3),dg_eta_fIV(i,j,1:3));
    end
end

% 8.5 - FACE V: CALCUL DU ROT. EN FONCTION DES DERIVEES ALPHA ET BETA
for i=1:nn,
    for j=1:nn,
        rot_fV(i,j,1:3)=cross(gxi_V(i,j,1:3),dg_xi_fV(i,j,1:3)) + cross(geta_V(i,j,1:3),dg_eta_fV(i,j,1:3));
    end
end

% 8.6 - FACE VI : CALCUL DU ROT. EN FONCTION DES DERIVEES ALPHA ET BETA 
for i=1:nn,
    for j=1:nn,
        rot_fVI(i,j,1:3)=cross(gxi_VI(i,j,1:3),dg_xi_fVI(i,j,1:3)) + cross(geta_VI(i,j,1:3),dg_eta_fVI(i,j,1:3));
    end
end

%% ETAPE 4 - GESTION DES BORDS
% Moyennage sur les bords des panels
% 1/2 SOMME ARRETES, 1/3 SOMME SOMMETS
% COMPONENT 1
uwk_I=rot_fI(1:nn,1:nn,1);uwk_II=rot_fII(1:nn,1:nn,1);uwk_III=rot_fIII(1:nn,1:nn,1);
uwk_IV=rot_fIV(1:nn,1:nn,1);uwk_V=rot_fV(1:nn,1:nn,1);uwk_VI=rot_fVI(1:nn,1:nn,1);
[rot_fI(1:nn,1:nn,1),rot_fII(1:nn,1:nn,1),rot_fIII(1:nn,1:nn,1),rot_fIV(1:nn,1:nn,1),rot_fV(1:nn,1:nn,1),rot_fVI(1:nn,1:nn,1)]=...
    ds72(uwk_I,uwk_II,uwk_III,uwk_IV,uwk_V,uwk_VI,n,nn);
% COMPONENT 2
uwk_I=rot_fI(1:nn,1:nn,2);uwk_II=rot_fII(1:nn,1:nn,2);uwk_III=rot_fIII(1:nn,1:nn,2);
uwk_IV=rot_fIV(1:nn,1:nn,2);uwk_V=rot_fV(1:nn,1:nn,2);uwk_VI=rot_fVI(1:nn,1:nn,2);
[rot_fI(1:nn,1:nn,2),rot_fII(1:nn,1:nn,2),rot_fIII(1:nn,1:nn,2),rot_fIV(1:nn,1:nn,2),rot_fV(1:nn,1:nn,2),rot_fVI(1:nn,1:nn,2)]=...
    ds72(uwk_I,uwk_II,uwk_III,uwk_IV,uwk_V,uwk_VI,n,nn);
% COMPONENT 3
uwk_I=rot_fI(1:nn,1:nn,3);uwk_II=rot_fII(1:nn,1:nn,3);uwk_III=rot_fIII(1:nn,1:nn,3);
uwk_IV=rot_fIV(1:nn,1:nn,3);uwk_V=rot_fV(1:nn,1:nn,3);uwk_VI=rot_fVI(1:nn,1:nn,3);
[rot_fI(1:nn,1:nn,3),rot_fII(1:nn,1:nn,3),rot_fIII(1:nn,1:nn,3),rot_fIV(1:nn,1:nn,3),rot_fV(1:nn,1:nn,3),rot_fVI(1:nn,1:nn,3)]=...
    ds72(uwk_I,uwk_II,uwk_III,uwk_IV,uwk_V,uwk_VI,n,nn);


end