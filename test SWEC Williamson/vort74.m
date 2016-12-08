function [vort_fI,vort_fII,vort_fIII,vort_fIV,vort_fV,vort_fVI]=...
    vort74(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI,n,nn)
% COMPUTATION SPHERICAL CURL (ON THE CUBED SPHERE GRID)

global alfa beta;
global betacr;
global alfa1;
global p kxi keta;
global gxi_I gxi_II gxi_III gxi_IV gxi_V gxi_VI;
global geta_I geta_II geta_III geta_IV geta_V geta_VI;
global gr_I gr_II gr_III gr_IV gr_V gr_VI;

vort_fI=zeros(nn,nn);
vort_fII=zeros(nn,nn);
vort_fIII=zeros(nn,nn);
vort_fIV=zeros(nn,nn);
vort_fV=zeros(nn,nn);
vort_fVI=zeros(nn,nn);

betaspline=zeros(nn,1);
alfaspline=zeros(nn,1);
funspline=zeros(nn,1);

%%  PHASE 1: ASSEMBLAGE DES 6 RESEAUX FRONT/BOTTOM
%% RESEAU 1 : ASSEMBLAGE DES DONNEES SUR LE RESEAU I-XI

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


    % CALCUL DES DERIVEES XI SUR LE RESEAU DE GRANDS CERCLES I-XI
    for jline1=1:nn, 
     funa7=va_fI(:,jline1); 
     test=kxi*funa7;
     funad8=p\test;
     vad_fI(:,jline1,kk)=funad8;
    end

end


%% RESEAU 2: ASSEMBLAGE DES DONNEES SUR LE RESEAU I-ETA

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

    % CALCUL DES DERIVEES ETA SUR LE RESEAU DE GRANDS CERCLES I-ETA
    for iline1=1:nn,
        funb1=vb_fI(iline1,:);
        test=keta*(funb1');
        funbd2=p\test;
        vbd_fI(iline1,:,kk)=funbd2; 
    end
    
end

%% RESEAU 3: ASSEMBLAGE DES DONNEES SUR LE RESEAU II-XI

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

%% RESEAU 4: ASSEMBLAGE DES DONNEES SUR LE RESEAU II-ETA

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

    % CALCUL DES DERIVEES ETA SUR LE RESEAU DE GRANDS CERCLES II-ETA
    for iline1=1:nn,
        funb2=vb_fII(iline1,:);
        test=keta*(funb2');
        funbd3=p\test;
        vbd_fII(iline1,:,kk)=funbd3; 
    end
end

%% RESEAU 5: ASSEMBLAGE DES DONNEES SUR LE RESEAU V-XI

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

    % CALCUL DES DERIVEES XI SUR LE RESEAU DE GRANDS CERCLES I-XI
    for jline1=1:nn,
        funa11=va_fV(:,jline1);
        test=kxi*funa11;
        funad12=p\test;
        vad_fV(:,jline1,kk)=funad12; 
    end
end

%% RESEAU 6: ASSEMBLAGE DES DONNEES SUR LE RESEAU V-ETA

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

    % CALCUL DES DERIVEES ETA SUR LE RESEAU DE GRANDS CERCLES V-ETA
    for iline1=1:nn,
        funb5=vb_fV(iline1,:);
        test=keta*(funb5');
        funbd6=p\test;
        vbd_fV(iline1,:,kk)=funbd6; 
    end
end

%% FIN ASSEMBLAGE DES DERIVEES HERMITIENNES SUR LES 6 RESEAUX DE CERCLES %%%%%%%%%%%%%%





%% ETAPE 2 - ASSEMBLAGE DES DERIVEES XI/ETA SUR CHACUNE DES 6 FACES 
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


%% ETAPE 3 - ASSEMBLAGE vortATIONNEL EN FONCTION DES DERIVEES XI/ETA

% 8.1 - FACE I: CALCUL DU vort. EN FONCTION DES DERIVEES ALPHA ET BETA
for i=1:nn,
    for j=1:nn,
        vort_fI(i,j)=dot(cross(gxi_I(i,j,1:3),dg_xi_fI(i,j,1:3)) + cross(geta_I(i,j,1:3), dg_eta_fI(i,j,1:3)), gr_I(i,j,1:3));
    end
end

% 8.2  FACE II: CALCUL DU vort. EN FONCTION DES DERIVEES ALPHA ET BETA FACE II
for i=1:nn,
    for j=1:nn,
        vort_fII(i,j)=dot(cross(gxi_II(i,j,1:3),dg_xi_fII(i,j,1:3)) + cross(geta_II(i,j,1:3),dg_eta_fII(i,j,1:3)), gr_II(i,j,1:3));
    end
end

% 8.3 -  FACE III: CALCUL DU vort. EN FONCTION DES DERIVEES ALPHA ET BETA FACE III
for i=1:nn
    for j=1:nn
        vort_fIII(i,j)=dot(cross(gxi_III(i,j,1:3),dg_xi_fIII(i,j,1:3)) + cross(geta_III(i,j,1:3), dg_eta_fIII(i,j,1:3)), gr_III(i,j,1:3));
    end
end

% 8.4 - FACE IV: CALCUL DU vort. EN FONCTION DES DERIVEES ALPHA ET BETA FACE IV
for i=1:nn,
    for j=1:nn,
        vort_fIV(i,j)=dot(cross(gxi_IV(i,j,1:3),dg_xi_fIV(i,j,1:3)) + cross(geta_IV(i,j,1:3),dg_eta_fIV(i,j,1:3)), gr_IV(i,j,1:3));
    end
end

% 8.5 - FACE V: CALCUL DU vort. EN FONCTION DES DERIVEES ALPHA ET BETA
for i=1:nn,
    for j=1:nn,
        vort_fV(i,j)=dot(cross(gxi_V(i,j,1:3),dg_xi_fV(i,j,1:3)) + cross(geta_V(i,j,1:3),dg_eta_fV(i,j,1:3)), gr_V(i,j,1:3));
    end
end

% 8.6 - FACE VI : CALCUL DU vort. EN FONCTION DES DERIVEES ALPHA ET BETA 
for i=1:nn,
    for j=1:nn,
        vort_fVI(i,j)=dot(cross(gxi_VI(i,j,1:3),dg_xi_fVI(i,j,1:3)) + cross(geta_VI(i,j,1:3),dg_eta_fVI(i,j,1:3)), gr_VI(i,j,1:3));
    end
end

%% ETAPE 4 - demi-somme + 1/3 somme de la divergence 
uwk_I=vort_fI(1:nn,1:nn);uwk_II=vort_fII(1:nn,1:nn);uwk_III=vort_fIII(1:nn,1:nn);
uwk_IV=vort_fIV(1:nn,1:nn);uwk_V=vort_fV(1:nn,1:nn);uwk_VI=vort_fVI(1:nn,1:nn);
[vort_fI(1:nn,1:nn),vort_fII(1:nn,1:nn),vort_fIII(1:nn,1:nn),vort_fIV(1:nn,1:nn),vort_fV(1:nn,1:nn),vort_fVI(1:nn,1:nn,1)]=...
    ds74(uwk_I,uwk_II,uwk_III,uwk_IV,uwk_V,uwk_VI,n,nn);
end