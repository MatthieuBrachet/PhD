function [rot_fI,rot_fII,rot_fIII,rot_fIV,rot_fV,rot_fVI]=...
    rot74(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI,n,nn)
% COMPUTATION SPHERICAL CURL (ON THE CUBED SPHERE GRID)

global alfa beta;
global betacr;
global alfa1;
global p kxi keta;
global gxi_I gxi_II gxi_III gxi_IV gxi_V gxi_VI;
global geta_I geta_II geta_III geta_IV geta_V geta_VI;


rot_fI=zeros(nn,nn,3);
rot_fII=zeros(nn,nn,3);
rot_fIII=zeros(nn,nn,3);
rot_fIV=zeros(nn,nn,3);
rot_fV=zeros(nn,nn,3);
rot_fVI=zeros(nn,nn,3);

betaspline=zeros(nn,1);
alfaspline=zeros(nn,1);
funspline=zeros(nn,1);

%%  PHASE 1: ASSEMBLAGE DES 6 RESEAUX FRONT/BOTTOM
%% RESEAU 1 : ASSEMBLAGE DES DONNEES SUR LE RESEAU I-ALPHA

for kk=1:3
    % Face I : 
    for jline1=1:nn, 
        va_fI(1:nn-1,jline1)=mfunfI(1:nn-1,jline1,kk);
    end

    % FaceII: 
    for i=1:nn-1 
        betaspline(1:nn)=beta(i,1:nn);
        funspline(1:nn)=mfunfII(i,1:nn,kk);
        ppspline=spline(betaspline,funspline);
        funbII1(1:nn)=ppval(ppspline,betacr(i,1:nn));
        va_fI(nn-1+i,1:nn)=funbII1(1:nn);
    end

    % Face III : 
    for jline1=1:nn,
        va_fI(2*nn-1:3*nn-3,jline1)=mfunfIII(1:nn-1,nn-jline1+1,kk);
    end

    % Face IV: 
    for i=1:nn-1, 
       betaspline(1:nn)=beta(i,1:nn);
       funspline(1:nn)=mfunfIV(i,1:nn,kk);
       ppspline=spline(betaspline,funspline);
       funbIV1(1:nn)=ppval(ppspline,betacr(i,nn+1-[1:nn]));
       va_fI(3*nn-3+i,1:nn)=funbIV1(1:nn);
    end


    % CALCUL DES DERIVEES ALFA SUR LE RESEAU DE GRANDS CERCLES I-ALFA
    for jline1=1:nn, 
     funa7=va_fI(:,jline1); 
     test=kxi*funa7;
     funad8=p\test;
     vad_fI(:,jline1,kk)=funad8;
    end
end


%% RESEAU 2: ASSEMBLAGE DES DONNEES SUR LE RESEAU I-BETA

for kk=1:3
    % Face I :
    for iline1=1:nn,
        vb_fI(iline1,1:nn-1)=mfunfI(iline1,1:nn-1,kk);
    end

    % FACE V: 
    for j=1:nn-1
        alfaspline(1:nn)=alfa(1:nn,j);
        funspline(1:nn)=mfunfV(1:nn,j,kk);
        ppspline=spline(alfaspline,funspline);
        funaV1(1:nn)=ppval(ppspline,alfa1(1:nn,j));
        vb_fI(1:nn,nn-1+j)=funaV1(1:nn);
    end

    % FACE III:
     for iline1=1:nn,
         vb_fI(iline1,2*nn-1:3*nn-3)=mfunfIII(iline1,nn:-1:2,kk); 
     end

    % Face VI : 
     for j=1:nn-1,
        alfaspline(1:nn)=alfa(1:nn,j);
        funspline(1:nn)=mfunfVI(1:nn,j,kk);
        ppspline=spline(alfaspline,funspline);
        funaVI1(1:nn)=ppval(ppspline,alfa1(nn+1-[1:nn],j));
        vb_fI(1:nn,3*nn-3+j)=funaVI1(1:nn);
     end

    % CALCUL DES DERIVEES BETA SUR LE RESEAU DE GRANDS CERCLES I-BETA
    for iline1=1:nn,
        funb1=vb_fI(iline1,:);
        test=keta*(funb1');
        funbd2=p\test;
        vbd_fI(iline1,:,kk)=funbd2; 
    end
end

%% RESEAU 3: ASSEMBLAGE DES DONNEES SUR LE RESEAU II-ALPHA

for kk=1:3
    % Face II :
    for jline1=1:nn,
        va_fII(1:nn-1,jline1)=mfunfII(1:nn-1,jline1,kk);
    end

    % FACE III: 
    for i=1:nn-1
        betaspline(1:nn)=beta(i,1:nn);
        funspline(1:nn)=mfunfIII(i,1:nn,kk);
        ppspline=spline(betaspline,funspline);
        funbIII1(1:nn)=ppval(ppspline,betacr(i,1:nn));
        va_fII(nn-1+i,1:nn)=funbIII1(1:nn);
    end

    % FACE IV: 
    for jline1=1:nn,
        va_fII(2*nn-1:3*nn-3,jline1)=mfunfIV(1:nn-1,nn-jline1+1,kk); 
    end

    % FACE I: 
    for i=1:nn-1, 
       betaspline(1:nn)=beta(i,1:nn);
       funspline(1:nn)=mfunfI(i,1:nn,kk);
       ppspline=spline(betaspline,funspline);
       funbI1(1:nn)=ppval(ppspline,betacr(i,nn+1-[1:nn]));
       va_fII(3*nn-3+i,1:nn)=funbI1(1:nn);
    end

    % CALCUL DES DERIVEES
    for jline1=1:nn,
     funa9=va_fII(:,jline1);
     test=kxi*funa9;
     funad10=p\test;
     vad_fII(:,jline1,kk)=funad10;
    end
end

%% RESEAU 4: ASSEMBLAGE DES DONNEES SUR LE RESEAU II-BETA

for kk=1:3
    % Face II : 
    for iline1=1:nn, 
        vb_fII(iline1,1:nn-1)=mfunfII(iline1,1:nn-1,kk);
    end

    % FACE V:
    for i=1:nn-1
        betaspline(1:nn)=beta(nn-i+1,1:nn);
        funspline(1:nn)=mfunfV(nn-i+1,1:nn,kk);
        ppspline=spline(betaspline,funspline);
        funbV1(1:nn)=ppval(ppspline,betacr(i,1:nn));
        vb_fII(1:nn,nn-1+i)=funbV1(1:nn);
    end

    % FACE IV:
    for iline1=1:nn,
        vb_fII(iline1,2*nn-1:3*nn-3)=mfunfIV(iline1,nn:-1:2,kk);
    end

    % FACE VI:
    for i=1:nn-1 
        betaspline(1:nn)=beta(i,1:nn);
        funspline(1:nn)=mfunfVI(i,1:nn,kk);
        ppspline=spline(betaspline,funspline);
        funbVI1(1:nn)=ppval(ppspline,betacr(i,1:nn));
        vb_fII(1:nn,3*nn-3+i)=funbVI1(1:nn);
    end

    % CALCUL DES DERIVEES BETA SUR LE RESEAU DE GRANDS CERCLES II-BETA
    for iline1=1:nn,
        funb2=vb_fII(iline1,:);
        test=keta*(funb2');
        funbd3=p\test;
        vbd_fII(iline1,:,kk)=funbd3; 
    end
end

%% RESEAU 5: ASSEMBLAGE DES DONNEES SUR LE RESEAU V-ALPHA

for kk=1:3
    % Face V
    for jline1=1:nn,
        va_fV(1:nn-1,jline1)=mfunfV(1:nn-1,jline1,kk);
    end

    % FACE II: 
    for j=1:nn-1, 
        alfaspline(1:nn)=alfa(1:nn,nn+1-j);
        funspline(1:nn)=mfunfII(1:nn,nn+1-j,kk);
        ppspline=spline(alfaspline,funspline);
        funaII1(1:nn)=ppval(ppspline,alfa1(1:nn,j));
        va_fV(nn-1+j,1:nn)=funaII1(1:nn);
    end 

    % FACE VI: 
    for jline1=1:nn,
        va_fV(2*nn-1:3*nn-3,jline1)=mfunfVI(nn:-1:2,jline1,kk); 
    end

    % FACE IV: 
    for j=1:nn-1, 
        alfaspline(1:nn)=alfa(1:nn,j);
        funspline(1:nn)=mfunfIV(1:nn,j,kk);
        ppspline=spline(alfaspline,funspline);
        funaIV1(1:nn)=ppval(ppspline,alfa1(1:nn,j));
        va_fV(3*nn-3+j,1:nn)=funaIV1(1:nn);
    end 

    % CALCUL DES DERIVEES ALPHA SUR LE RESEAU DE GRANDS CERCLES I-ALPHA
    for jline1=1:nn,
        funa11=va_fV(:,jline1);
        test=kxi*funa11;
        funad12=p\test;
        vad_fV(:,jline1,kk)=funad12; 
    end
end

%% RESEAU 6: ASSEMBLAGE DES DONNEES SUR LE RESEAU V-BETA

for kk=1:3
    % Face V : 
    for iline1=1:nn, 
        vb_fV(iline1,1:nn-1)=mfunfV(iline1,1:nn-1,kk);
    end

    % Face III :
    for j=1:nn-1, 
        alfaspline(1:nn)=alfa(1:nn,nn+1-j);
        funspline(1:nn)=mfunfIII(1:nn,nn+1-j,kk);
        ppspline=spline(alfaspline,funspline);
        funaIII1(1:nn)=ppval(ppspline,alfa1(nn+1-[1:nn],j)); 
        vb_fV(1:nn,nn-1+j)=funaIII1(1:nn);
    end 

    % Face VI
    for iline1=1:nn,
        vb_fV(iline1,2*nn-1:3*nn-3)=mfunfVI(nn-iline1+1,1:nn-1,kk);
    end

    % Face I :
    for j=1:nn-1,
        alfaspline(1:nn)=alfa(1:nn,j);
        funspline(1:nn)=mfunfI(1:nn,j,kk);
        ppspline=spline(alfaspline,funspline);
        funaI1(1:nn)=ppval(ppspline,alfa1(nn+1-[1:nn],j)); 
        vb_fV(1:nn,3*nn-3+j)=funaI1(1:nn);
    end 

    % CALCUL DES DERIVEES BETA SUR LE RESEAU DE GRANDS CERCLES V-BETA
    for iline1=1:nn,
        funb5=vb_fV(iline1,:);
        test=keta*(funb5');
        funbd6=p\test;
        vbd_fV(iline1,:,kk)=funbd6; 
    end
end

%% FIN ASSEMBLAGE DES DERIVEES HERMITIENNES SUR LES 6 RESEAUX DE CERCLES %%%%%%%%%%%%%%





%% ETAPE 2 - ASSEMBLAGE DES DERIVEES ALPHA/BETA SUR CHACUNE DES 6 FACES 
% DEDUITES DES DERIVEES HERMITIENNES SUR LES 6 RESEAUX DE GRANDS CERCLES.

% FACE I
for i=1:nn,
    for j=1:nn,
        dg_xi_fI(i,j,1:3) = vad_fI(i,j,1:3);
        dg_eta_fI(i,j,1:3) = vbd_fI(i,j,1:3);
    end
end
% FACE II
for i=1:nn,
    for j=1:nn,
        dg_xi_fII(i,j,1:3)=vad_fII(i,j,1:3);
        dg_eta_fII(i,j,1:3)=vbd_fII(i,j,1:3);
    end
end
% FACE III
for i=1:nn,
    for j=1:nn,
        dg_xi_fIII(i,j,1:3)=vad_fI(2*(nn-1)+i,nn-j+1,1:3);
        % -
        dg_eta_fIII(i,j,1:3)=-vbd_fI(i,2*(nn-1)+nn-j+1,1:3); 
    end
end
% FACE IV
for i=1:nn,
    for j=1:nn,
        dg_xi_fIV(i,j,1:3)=vad_fII(2*(nn-1)+i,nn-j+1,1:3);
        % -
        dg_eta_fIV(i,j,1:3)=-vbd_fII(i,2*(nn-1)+nn-j+1,1:3);
    end
end
% FACE V
for i=1:nn,
    for j=1:nn,
        dg_xi_fV(i,j,1:3)=vad_fV(i,j,1:3);
        dg_eta_fV(i,j,1:3)=vbd_fV(i,j,1:3);
    end
end
% FACE VI
for i=1:nn,
    for j=1:nn,
        % -
        dg_xi_fVI(i,j,1:3)=-vad_fV(2*(nn-1)+nn-i+1,j,1:3);
        dg_eta_fVI(i,j,1:3)=vbd_fV(nn-i+1,2*(nn-1)+j,1:3);
    end
end


%% ETAPE 3 - ASSEMBLAGE GRADIENTS EN FONCTION DES DERIVEES ALPHA/BETA

% 8.1 - FACE I: CALCUL DU GRADIENT EN FONCTION DES DERIVEES ALPHA ET BETA
% FACE

for i=1:nn,
    for j=1:nn,
        rot_fI(i,j,1:3)=cross(gxi_I(i,j,1:3),dg_xi_fI(i,j,1:3)) + cross(geta_I(i,j,1:3), dg_eta_fI(i,j,1:3));
        % rot_fI(i,j,1:3)=pdt_vect(gxi_I(i,j,1:3),dg_xi1(i,j,1:3)) + pdt_vect(geta_I(i,j,1:3), dg_eta1(i,j,1:3));
    end
end

% 8.2  FACE II: CALCUL DU GRADIENT EN FONCTION DES DERIVEES ALPHA ET BETA FACE II

for i=1:nn,
    for j=1:nn,
        rot_fII(i,j,1:3)=cross(gxi_II(i,j,1:3),dg_xi_fII(i,j,1:3)) + cross(geta_II(i,j,1:3),dg_eta_fII(i,j,1:3));
        % rot_fII(i,j,1:3)=pdt_vect(gxi_II(i,j,1:3),dg_xi2(i,j,1:3)) + pdt_vect(geta_II(i,j,1:3),dg_eta2(i,j,1:3));
    end
end

% 8.3 -  FACE III: CALCUL DU GRADIENT EN FONCTION DES DERIVEES ALPHA ET BETA FACE III

% FACE III: CALCUL DU GRADIENT EN FONCTION DES DERIVEES ALPHA ET BETA FACE
for i=1:nn
    for j=1:nn
        rot_fIII(i,j,1:3)=cross(gxi_III(i,j,1:3),dg_xi_fIII(i,j,1:3)) + cross(geta_III(i,j,1:3), dg_eta_fIII(i,j,1:3));
        % rot_fIII(i,j,1:3)=pdt_vect(gxi_III(i,j,1:3),dg_xi3(i,j,1:3)) + pdt_vect(geta_III(i,j,1:3), dg_eta3(i,j,1:3));
    end
end

% 8.4 - FACE IV: CALCUL DU GRADIENT EN FONCTION DES DERIVEES ALPHA ET BETA FACE IV

for i=1:nn,
    for j=1:nn,
        rot_fIV(i,j,1:3)=cross(gxi_IV(i,j,1:3),dg_xi_fIV(i,j,1:3)) + cross(geta_IV(i,j,1:3),dg_eta_fIV(i,j,1:3));
        % rot_fIV(i,j,1:3)=pdt_vect(gxi_IV(i,j,1:3),dg_xi4(i,j,1:3)) + pdt_vect(geta_IV(i,j,1:3),dg_eta4(i,j,1:3));
    end
end

% 8.5 - FACE V: CALCUL DU GRADIENT EN FONCTION DES DERIVEES ALPHA ET BETA

% FACE V: CALCUL DU GRADIENT EN FONCTION DES DERIVEES ALPHA ET BETA FACE V
for i=1:nn,
    for j=1:nn,
        rot_fV(i,j,1:3)=cross(gxi_V(i,j,1:3),dg_xi_fV(i,j,1:3)) + cross(geta_V(i,j,1:3),dg_eta_fV(i,j,1:3));
        % rot_fV(i,j,1:3)=pdt_vect(gxi_V(i,j,1:3),dg_xi5(i,j,1:3)) + pdt_vect(geta_V(i,j,1:3),dg_eta5(i,j,1:3));
    end
end

% 8.6 - FACE VI : CALCUL DU GRADIENT EN FONCTION DES DERIVEES ALPHA ET BETA 

% FACE VI: CALCUL DU GRADIENT EN FONCTION DES DERIVEES ALPHA ET BETA FACE
for i=1:nn,
    for j=1:nn,
        rot_fVI(i,j,1:3)=cross(gxi_VI(i,j,1:3),dg_xi_fVI(i,j,1:3)) + cross(geta_VI(i,j,1:3),dg_eta_fVI(i,j,1:3));
        % rot_fVI(i,j,1:3)=pdt_vect(gxi_VI(i,j,1:3),dg_xi6(i,j,1:3)) + pdt_vect(geta_VI(i,j,1:3),dg_eta6(i,j,1:3));
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