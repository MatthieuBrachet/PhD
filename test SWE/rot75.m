function [rot_fI,rot_fII,rot_fIII,rot_fIV,rot_fV,rot_fVI]=...
    rot75(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI,n,nn)
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
%% x-comp
% RESEAU 1 : ASSEMBLAGE DES DONNEES SUR LE RESEAU I-ALPHA

funfI=mfunfI(:,:,1);
funfII=mfunfII(:,:,1);
funfIII=mfunfIII(:,:,1);
funfIV=mfunfIV(:,:,1);
funfV=mfunfV(:,:,1);
funfVI=mfunfVI(:,:,1);

    % Face I : 
    for jline1=1:nn, 
        va_fI(1:nn-1,jline1)=funfI(1:nn-1,jline1);
    end

    % FaceII: 
    for i=1:nn-1 
        betaspline(1:nn)=beta(i,1:nn);
        funspline(1:nn)=funfII(i,1:nn);
        ppspline=spline(betaspline,funspline);
        funbII1(1:nn)=ppval(ppspline,betacr(i,1:nn));
        va_fI(nn-1+i,1:nn)=funbII1(1:nn);
    end

    % Face III : 
    for jline1=1:nn,
        va_fI(2*nn-1:3*nn-3,jline1)=funfIII(1:nn-1,nn-jline1+1);
    end

    % Face IV: 
    for i=1:nn-1, 
       betaspline(1:nn)=beta(i,1:nn);
       funspline(1:nn)=funfIV(i,1:nn);
       ppspline=spline(betaspline,funspline);
       funbIV1(1:nn)=ppval(ppspline,betacr(i,nn+1-[1:nn]));
       va_fI(3*nn-3+i,1:nn)=funbIV1(1:nn);
    end


    % CALCUL DES DERIVEES ALFA SUR LE RESEAU DE GRANDS CERCLES I-ALFA
    for jline1=1:nn, 
     funa7=va_fI(:,jline1); 
     test=kxi*funa7;
     funad8=p\test;
     vad_fI(:,jline1,1)=funad8;
    end

% RESEAU 2: ASSEMBLAGE DES DONNEES SUR LE RESEAU I-BETA

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

    % Face VI : 
     for j=1:nn-1,
        alfaspline(1:nn)=alfa(1:nn,j);
        funspline(1:nn)=funfVI(1:nn,j);
        ppspline=spline(alfaspline,funspline);
        funaVI1(1:nn)=ppval(ppspline,alfa1(nn+1-[1:nn],j));
        vb_fI(1:nn,3*nn-3+j)=funaVI1(1:nn);
     end

    % CALCUL DES DERIVEES BETA SUR LE RESEAU DE GRANDS CERCLES I-BETA
    for iline1=1:nn,
        funb1=vb_fI(iline1,:);
        test=keta*(funb1');
        funbd2=p\test;
        vbd_fI(iline1,:,1)=funbd2; 
    end
    
    
%% y-comp
% RESEAU 1 : ASSEMBLAGE DES DONNEES SUR LE RESEAU I-ALPHA

funfI=mfunfI(:,:,2);
funfII=mfunfII(:,:,2);
funfIII=mfunfIII(:,:,2);
funfIV=mfunfIV(:,:,2);
funfV=mfunfV(:,:,2);
funfVI=mfunfVI(:,:,2);

    % Face I : 
    for jline1=1:nn, 
        va_fI(1:nn-1,jline1)=funfI(1:nn-1,jline1);
    end

    % FaceII: 
    for i=1:nn-1 
        betaspline(1:nn)=beta(i,1:nn);
        funspline(1:nn)=funfII(i,1:nn);
        ppspline=spline(betaspline,funspline);
        funbII1(1:nn)=ppval(ppspline,betacr(i,1:nn));
        va_fI(nn-1+i,1:nn)=funbII1(1:nn);
    end

    % Face III : 
    for jline1=1:nn,
        va_fI(2*nn-1:3*nn-3,jline1)=funfIII(1:nn-1,nn-jline1+1);
    end

    % Face IV: 
    for i=1:nn-1, 
       betaspline(1:nn)=beta(i,1:nn);
       funspline(1:nn)=funfIV(i,1:nn);
       ppspline=spline(betaspline,funspline);
       funbIV1(1:nn)=ppval(ppspline,betacr(i,nn+1-[1:nn]));
       va_fI(3*nn-3+i,1:nn)=funbIV1(1:nn);
    end


    % CALCUL DES DERIVEES ALFA SUR LE RESEAU DE GRANDS CERCLES I-ALFA
    for jline1=1:nn, 
     funa7=va_fI(:,jline1); 
     test=kxi*funa7;
     funad8=p\test;
     vad_fI(:,jline1,2)=funad8;
    end

% RESEAU 2: ASSEMBLAGE DES DONNEES SUR LE RESEAU I-BETA

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

    % Face VI : 
     for j=1:nn-1,
        alfaspline(1:nn)=alfa(1:nn,j);
        funspline(1:nn)=funfVI(1:nn,j);
        ppspline=spline(alfaspline,funspline);
        funaVI1(1:nn)=ppval(ppspline,alfa1(nn+1-[1:nn],j));
        vb_fI(1:nn,3*nn-3+j)=funaVI1(1:nn);
     end

    % CALCUL DES DERIVEES BETA SUR LE RESEAU DE GRANDS CERCLES I-BETA
    for iline1=1:nn,
        funb1=vb_fI(iline1,:);
        test=keta*(funb1');
        funbd2=p\test;
        vbd_fI(iline1,:,2)=funbd2; 
    end

%% z-comp
% RESEAU 1 : ASSEMBLAGE DES DONNEES SUR LE RESEAU I-ALPHA

funfI=mfunfI(:,:,3);
funfII=mfunfII(:,:,3);
funfIII=mfunfIII(:,:,3);
funfIV=mfunfIV(:,:,3);
funfV=mfunfV(:,:,3);
funfVI=mfunfVI(:,:,3);

    % Face I : 
    for jline1=1:nn, 
        va_fI(1:nn-1,jline1)=funfI(1:nn-1,jline1);
    end

    % FaceII: 
    for i=1:nn-1 
        betaspline(1:nn)=beta(i,1:nn);
        funspline(1:nn)=funfII(i,1:nn);
        ppspline=spline(betaspline,funspline);
        funbII1(1:nn)=ppval(ppspline,betacr(i,1:nn));
        va_fI(nn-1+i,1:nn)=funbII1(1:nn);
    end

    % Face III : 
    for jline1=1:nn,
        va_fI(2*nn-1:3*nn-3,jline1)=funfIII(1:nn-1,nn-jline1+1);
    end

    % Face IV: 
    for i=1:nn-1, 
       betaspline(1:nn)=beta(i,1:nn);
       funspline(1:nn)=funfIV(i,1:nn);
       ppspline=spline(betaspline,funspline);
       funbIV1(1:nn)=ppval(ppspline,betacr(i,nn+1-[1:nn]));
       va_fI(3*nn-3+i,1:nn)=funbIV1(1:nn);
    end


    % CALCUL DES DERIVEES ALFA SUR LE RESEAU DE GRANDS CERCLES I-ALFA
    for jline1=1:nn, 
     funa7=va_fI(:,jline1); 
     test=kxi*funa7;
     funad8=p\test;
     vad_fI(:,jline1,3)=funad8;
    end

% RESEAU 2: ASSEMBLAGE DES DONNEES SUR LE RESEAU I-BETA

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

    % Face VI : 
     for j=1:nn-1,
        alfaspline(1:nn)=alfa(1:nn,j);
        funspline(1:nn)=funfVI(1:nn,j);
        ppspline=spline(alfaspline,funspline);
        funaVI1(1:nn)=ppval(ppspline,alfa1(nn+1-[1:nn],j));
        vb_fI(1:nn,3*nn-3+j)=funaVI1(1:nn);
     end

    % CALCUL DES DERIVEES BETA SUR LE RESEAU DE GRANDS CERCLES I-BETA
    for iline1=1:nn,
        funb1=vb_fI(iline1,:);
        test=keta*(funb1');
        funbd2=p\test;
        vbd_fI(iline1,:,3)=funbd2; 
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




%% ETAPE 3 - ASSEMBLAGE GRADIENTS EN FONCTION DES DERIVEES ALPHA/BETA

% 8.1 - FACE I: CALCUL DU GRADIENT EN FONCTION DES DERIVEES ALPHA ET BETA
% FACE

for i=1:nn,
    for j=1:nn,
        rot_fI(i,j,1:3)=cross(gxi_I(i,j,1:3),dg_xi_fI(i,j,1:3)) + cross(geta_I(i,j,1:3), dg_eta_fI(i,j,1:3));
    end
end



end