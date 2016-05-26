function [grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
    gr72(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn)
global mm na nb;
global radius;
global xi eta dxi deta xx yy delta deltab;
global alfa beta;
global alfacr betacr;
global alfa1;
global alfag betag;
global p k;
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;
global gxi_I gxi_II gxi_III gxi_IV gxi_V gxi_VI;
global geta_I geta_II geta_III geta_IV geta_V geta_VI;
global ite aaa bbb itestop;
global ftr;
% COMPUTATION SPHERICAL GRADIENT
grad_fI=zeros(nn,nn,3);grad_fII=zeros(nn,nn,3);grad_fIII=zeros(nn,nn,3);
grad_fIV=zeros(nn,nn,3);grad_fV=zeros(nn,nn,3);grad_fVI=zeros(nn,nn,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
betaspline=zeros(nn,1);
alfaspline=zeros(nn,1);
funspline=zeros(nn,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHASE 1: ASSEMBLAGE DES 6 RESEAUX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESEAU 1 : ASSEMBLAGE DES DONNEES SUR LE RESEAU I-ALPHA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------   
% tic;
%
va_fI=zeros(4*(nn-1),nn); % VALEURS =DONNEES RESEAU FACE I SELON ALPHA
funbII1=zeros(nn,1);
funbIV1=zeros(nn,1);
% Face I : TRANSFERT OF DATA OF FACE I 
for jline1=1:nn, % boucle sur les lignes iso-eta du reseaux Ia
    va_fI(1:nn-1,jline1)=funfI(1:nn-1,jline1);
end
% FACE II: SPLINE INTERPOLATION A L'AIDE DES ANGLES BETA
% INVERSION OF THE LOOP ON THE ISO-ETA LINES AND ON THE XI OF FACE II
for i=1:nn-1 % LOOP ON THE XI OF FACE II
    betaspline(1:nn)=beta(i,1:nn);
    funspline(1:nn)=funfII(i,1:nn);
    ppspline=spline(betaspline,funspline);
    funbII1(1:nn)=ppval(ppspline,betacr(i,1:nn));
    va_fI(nn-1+i,1:nn)=funbII1(1:nn);
end
for jline1=1:nn,
    % FACE III: TRANSFERT OF DATA OF FACE III
    va_fI(2*nn-1:3*nn-3,jline1)=funfIII(1:nn-1,nn-jline1+1); % symetrie sur adresses en j !
end
% FACE IV: SPLINE INTERPOLATION A L'AIDE DES ANGLES BETA
% INVERSION OF THE LOOP ON THE ISO-ETA LINES AND ON THE XI OF FACE IV
for i=1:nn-1, % LOOP ON THE XI OF FACE IV
   betaspline(1:nn)=beta(i,1:nn);
   funspline(1:nn)=funfIV(i,1:nn);
   ppspline=spline(betaspline,funspline);
   funbIV1(1:nn)=ppval(ppspline,betacr(i,nn+1-[1:nn]));
   va_fI(3*nn-3+i,1:nn)=funbIV1(1:nn);
end
% FILTRAGE
% for jline1=1:nn,
%    va_fI(1:na,jline1)=ftr*va_fI(1:na,jline1);
% end
% --------------------------------------------------
alfag4=zeros(na,1); % ANGLE ALONG EACH GREAT CIRCLE 
%xa=zeros(na,1);ya=zeros(na,1);za=zeros(na,1); % CARTESIAN COORDINATES IN COORDINATES ALFA/ETA (FACE I)
funadex2=zeros(na,nn);
vad_fI=zeros(na,nn);
for jline1=1:nn, % EACH GREAT CIRCLE OF  NETWORK I-ALPHA
 alfag4=alfag(:,jline1);
 dalfag4=k*alfag4;
 dalfag4(1)=dalfag4(1)+2*pi;
 dalfag4(na)=dalfag4(na)+2*pi;
 alfag4d=3*(p\dalfag4);
 %
 funa7=va_fI(:,jline1); % TABLEAU DE TRAVAIL VALEURS
 funad7=3*((p\k)*funa7); % TABLEAU DE TRAVAIL DERIVEES
 funad8=funad7./alfag4d; % AUTRE TABLEAU DE TRAVAIL DERIVEES
 vad_fI(:,jline1)=funad8; % VALUE OF THE DERIVATIVE WITH RESPECT TO ANGLE ALFA.
end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESEAU 2: ASSEMBLAGE DES DONNEES SUR LE RESEAU I-BETA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
vb_fI=zeros(nn,4*(nn-1)); % VALEURS =DONNEES RESEAU FACE I SELON BETA
%
funaV1=zeros(nn,1);
funaVI1=zeros(nn,1);
% Face I : TRANSFERT OF DATA OF FACE I
for iline1=1:nn, % boucle sur les lignes iso-xi du reseau I-beta
    vb_fI(iline1,1:nn-1)=funfI(iline1,1:nn-1);
end
% FACE V: SPLINE INTERPOLATION A L'AIDE DES ANGLES ALFA
% INVERSION OF THE LOOP ON THE ISO-XI LINES AND ON THE ETA OF FACE V
for j=1:nn-1 % LOOP ON THE ETA OF FACE V
    alfaspline(1:nn)=alfa(1:nn,j);
    funspline(1:nn)=funfV(1:nn,j);
    ppspline=spline(alfaspline,funspline);
    funaV1(1:nn)=ppval(ppspline,alfa1(1:nn,j));
    vb_fI(1:nn,nn-1+j)=funaV1(1:nn);
end
    % FACE III: TRANSFERT OF DATA
    % ---------------------------
 for iline1=1:nn,
     vb_fI(iline1,2*nn-1:3*nn-3)=funfIII(iline1,nn:-1:2); % symetrie sur adresses en i + inversion en j !
 end
     % FACE VI: SPLINE INTERPOLATION A L'AIDE DES ANGLES ALFA
    % --------------------------------------------------------
    % LOADING DATA FOR THE SPLINE INTERPOLATION :
 for j=1:nn-1, % LOOP ON THE ETA OF FACE VI
    alfaspline(1:nn)=alfa(1:nn,j);
    funspline(1:nn)=funfVI(1:nn,j);
    ppspline=spline(alfaspline,funspline);
    funaVI1(1:nn)=ppval(ppspline,alfa1(nn+1-[1:nn],j));
    vb_fI(1:nn,3*nn-3+j)=funaVI1(1:nn);
 end
 % FILTRAGE
% for iline1=1:nn,
%    vb_fI(iline1,1:nb)=vb_fI(iline1,1:nb)*ftr;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCUL DES DERIVEES BETA SUR LE RESEAU DE GRANDS CERCLES I-BETA
betag1=zeros(1,nb);
xb=zeros(1,nb);yb=zeros(1,nb);zb=zeros(1,nb);
funbdex3=zeros(nn,nb);
vbd_fI=zeros(nn,nb);
for iline1=1:nn,
 % betag1=betag_fI(iline1,:);
 betag1=betag(iline1,:);
 dbetag1=k*betag1';
 dbetag1(1)=dbetag1(1)+2*pi;
 dbetag1(na)=dbetag1(na)+2*pi;
 betag1d=3*(p\dbetag1);
 %
 funb1=vb_fI(iline1,:);
 funbd1=3*((p\k)*funb1');
 funbd2=funbd1./betag1d;
 vbd_fI(iline1,:)=funbd2; % VALUE OF THE DERIVATIVE WITH RESPECT TO ANGLE ALFA.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESEAU 3: ASSEMBLAGE DES DONNEES SUR LE RESEAU II-ALPHA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tic
va_fII=zeros(4*(nn-1),nn); % VALEURS =DONNEES RESEAU FACE II SELON ALPHA
funbIII1=zeros(nn,1);
funbI1=zeros(nn,1);
% Face II : TRANSFERT OF DATA OF FACE II
for jline1=1:nn, % boucle sur les lignes iso-eta du reseaux IIa
%     jbar= jline1-mm;
%     etabar= jbar*deta;
    va_fII(1:nn-1,jline1)=funfII(1:nn-1,jline1);
end
% FACE III: SPLINE INTERPOLATION A L'AIDE DES ANGLES BETA
% INVERSION OF THE LOOP ON THE ISO-ETA LINES AND ON THE XI OF FACE II
for i=1:nn-1 % LOOP ON THE XI OF FACE II
    betaspline(1:nn)=beta(i,1:nn);
    funspline(1:nn)=funfIII(i,1:nn);
    ppspline=spline(betaspline,funspline);
    funbIII1(1:nn)=ppval(ppspline,betacr(i,1:nn));
    va_fII(nn-1+i,1:nn)=funbIII1(1:nn);
end
% FACE IV: TRANSFERT OF DATA
    % ---------------------------
for jline1=1:nn,
    va_fII(2*nn-1:3*nn-3,jline1)=funfIV(1:nn-1,nn-jline1+1); % symetrie sur adresses en j !    
end
% % FACE I: SPLINE INTERPOLATION A L'AIDE DES ANGLES BETA
% INVERSION OF THE LOOP ON THE ISO-ETA LINES AND ON THE XI OF FACE I
for i=1:nn-1, % LOOP ON THE XI OF FACE I============
   betaspline(1:nn)=beta(i,1:nn);
   funspline(1:nn)=funfI(i,1:nn);
   ppspline=spline(betaspline,funspline);
   funbI1(1:nn)=ppval(ppspline,betacr(i,nn+1-[1:nn]));
   va_fII(3*nn-3+i,1:nn)=funbI1(1:nn);
end
% FILTRAGE
% for jline1=1:nn,
%    va_fII(1:na,jline1)=ftr*va_fII(1:na,jline1);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
alfag5=zeros(na,1);
xa=zeros(na,1);ya=zeros(na,1);za=zeros(na,1);
funadex4=zeros(na,nn);
vad_fII=zeros(na,nn);
for jline1=1:nn,
 alfag5=alfag(:,jline1);
 dalfag5=k*alfag5;
 dalfag5(1)=dalfag5(1)+2*pi;
 dalfag5(na)=dalfag5(na)+2*pi;
 alfag5d=3*(p\dalfag5);
 %
 funa9=va_fII(:,jline1);
 funad9=3*((p\k)*funa9);
 funad10=funad9./alfag5d;
 vad_fII(:,jline1)=funad10; % VALUE OF THE DERIVATIVE WITH RESPECT TO ANGLE ALFA.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESEAU 4: ASSEMBLAGE DES DONNEES SUR LE RESEAU II-BETA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tic;
vb_fII=zeros(nn,4*(nn-1)); % VALEURS =DONNEES RESEAU FACE II SELON BETA
funbV1=zeros(nn,1);
funbVI1=zeros(nn,1);
% Face II : TRANSFERT OF DATA OF FACE II
for iline1=1:nn, % boucle sur les lignes iso-xi du reseau II-beta
    vb_fII(iline1,1:nn-1)=funfII(iline1,1:nn-1);
end
% FACE V: SPLINE INTERPOLATION A L'AIDE DES ANGLES BETA
% INVERSION OF THE LOOP ON THE ISO-ETA LINES AND ON THE XI OF FACE V
for i=1:nn-1 % LOOP ON THE XI OF FACE V. i= order of encountering !!!
    betaspline(1:nn)=beta(nn-i+1,1:nn);
    funspline(1:nn)=funfV(nn-i+1,1:nn);
    ppspline=spline(betaspline,funspline);
    funbV1(1:nn)=ppval(ppspline,betacr(i,1:nn));
    vb_fII(1:nn,nn-1+i)=funbV1(1:nn);
end
% FACE III: TRANSFERT OF DATA OF FACE III
for iline1=1:nn,
    vb_fII(iline1,2*nn-1:3*nn-3)=funfIV(iline1,nn:-1:2); % symetrie sur adresses en i + inversion en j !
end
% FACE VI: SPLINE INTERPOLATION A L'AIDE DES ANGLES BETA
% INVERSION OF THE LOOP ON THE ISO-ETA LINES AND ON THE XI OF FACE Vi
for i=1:nn-1 % LOOP ON THE XI OF FACE VI. i= order of encountering !!!
    betaspline(1:nn)=beta(i,1:nn);
    funspline(1:nn)=funfVI(i,1:nn);
    ppspline=spline(betaspline,funspline);
    funbVI1(1:nn)=ppval(ppspline,betacr(i,1:nn));
    vb_fII(1:nn,3*nn-3+i)=funbVI1(1:nn);
end
 % FILTRAGE
% for iline1=1:nn,
%    vb_fII(iline1,1:nb)=vb_fII(iline1,1:nb)*ftr;
% end
%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCUL DES DERIVEES BETA SUR LE RESEAU DE GRANDS CERCLES II-BETA
betag2=zeros(1,nb);
xb=zeros(1,nb);yb=zeros(1,nb);zb=zeros(1,nb);
funbdex4=zeros(nn,nb);
vbd_fII=zeros(nn,nb);
for iline1=1:nn,
 betag2=betag(iline1,:);
 dbetag2=k*betag2';
 dbetag2(1)=dbetag2(1)+2*pi;
 dbetag2(na)=dbetag2(na)+2*pi;
 betag2d=3*(p\dbetag2);
 %
 funb2=vb_fII(iline1,:);
 funbd2=3*((p\k)*funb2');
 funbd3=funbd2./betag2d;
 vbd_fII(iline1,:)=funbd3; % VALUE OF THE DERIVATIVE WITH RESPECT TO ANGLE ALFA.
end
if ite==itestop,
    aaa=vb_fII(9,1:nb)';
    bbb=vbd_fII(9,1:nb)';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESEAU 5: ASSEMBLAGE DES DONNEES SUR LE RESEAU V-ALPHA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
va_fV=zeros(4*(nn-1),nn); % VALEURS =DONNEES RESEAU FACE V SELON ALPHA
funaII1=zeros(nn,1);
funaIV1=zeros(nn,1);
for jline1=1:nn, %% Face V : TRANSFERT OF DATA OF FACE V
    va_fV(1:nn-1,jline1)=funfV(1:nn-1,jline1);
end
% FACE II: SPLINE INTERPOLATION  A L'AIDE DES ANGLES ALFA
% INVERSION OF THE LOOP ON THE ISO-ETA LINES OF NETWORK VA 
% AND THE ETA LINES OF FACE II
for j=1:nn-1, % LOOP ON THE ETA-LINES OF FACE II. j= order of encountering !!!
    alfaspline(1:nn)=alfa(1:nn,nn+1-j);
    funspline(1:nn)=funfII(1:nn,nn+1-j);
    ppspline=spline(alfaspline,funspline);
    funaII1(1:nn)=ppval(ppspline,alfa1(1:nn,j)); % VERIFIER SI C'EST BIEN "J" SUR LE PLAN LOGIQUE: je ne suis pas sur!!!!!
    va_fV(nn-1+j,1:nn)=funaII1(1:nn);
end 
% FACE VI: TRANSFERT OF DATA
for jline1=1:nn,
    jbar= jline1-mm;
    etabar= jbar*deta;
    va_fV(2*nn-1:3*nn-3,jline1)=funfVI(nn:-1:2,jline1); % symetrie a bien regarder!
end
% FACE IV: SPLINE INTERPOLATION  A L'AIDE DES ANGLES ALFA
% INVERSION OF THE LOOP ON THE ISO-ETA LINES OF NETWORK VA 
% AND THE ETA LINES OF FACE IV
for j=1:nn-1, % LOOP ON THE ETA-LINES OF FACE IV. j= order of encountering !!!
    alfaspline(1:nn)=alfa(1:nn,j);
    funspline(1:nn)=funfIV(1:nn,j);
    ppspline=spline(alfaspline,funspline);
    funaIV1(1:nn)=ppval(ppspline,alfa1(1:nn,j)); % VERIFIER SI C'EST BIEN "J" SUR LE PLAN LOGIQUE: je ne suis pas sur!!!!!
    va_fV(3*nn-3+j,1:nn)=funaIV1(1:nn);
end 
% FILTRAGE
% for jline1=1:nn,
%    va_fV(1:na,jline1)=ftr*va_fV(1:na,jline1);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCUL DES DERIVEES ALPHA SUR LE RESEAU DE GRANDS CERCLES I-ALPHA
%jline=mm;% TO BE FIXED HERE= SELECTION OF THE ISOXI LINE, jline in 1,nn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
alfag6=zeros(na,1);
xa=zeros(na,1);ya=zeros(na,1);za=zeros(na,1);
funadex5=zeros(na,nn);
vad_fV=zeros(na,nn);
for jline1=1:nn,
 alfag6=alfag(1:na,jline1);
 dalfag6=k*alfag6;
 dalfag6(1)=dalfag6(1)+2*pi;
 dalfag6(na)=dalfag6(na)+2*pi;
 alfag6d=3*(p\dalfag6);
 %
 funa11=va_fV(:,jline1);
 funad11=3*((p\k)*funa11);
 funad12=funad11./alfag6d;
 vad_fV(:,jline1)=funad12; % VALUE OF THE DERIVATIVE WITH RESPECT TO ANGLE ALFA.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESEAU 6: ASSEMBLAGE DES DONNEES SUR LE RESEAU V-BETA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tic;
vb_fV=zeros(nn,4*(nn-1)); % VALEURS =DONNEES RESEAU FACE V SELON BETA
funaIII1=zeros(nn,1);
funaI1=zeros(nn,1);
% Face V : TRANSFERT OF DATA selon beta
for iline1=1:nn, % boucle sur les lignes iso-xi du reseau V-beta
    vb_fV(iline1,1:nn-1)=funfV(iline1,1:nn-1);
end
for j=1:nn-1, % LOOP ON THE ETA-LINES OF FACE III. j= order of encountering !!!
    alfaspline(1:nn)=alfa(1:nn,nn+1-j);
    funspline(1:nn)=funfIII(1:nn,nn+1-j);
    ppspline=spline(alfaspline,funspline);
    funaIII1(1:nn)=ppval(ppspline,alfa1(nn+1-[1:nn],j)); % VERIFIER SI C'EST BIEN "NN+1-[1,NN]" SUR LE PLAN LOGIQUE: je ne suis pas sur!!!!!
    vb_fV(1:nn,nn-1+j)=funaIII1(1:nn);
end 
for iline1=1:nn,
%     ibar= iline1-mm;
%     xibar= ibar*dxi; 
% FACE VI: TRANSFERT OF DATA
    % ---------------------------
    vb_fV(iline1,2*nn-1:3*nn-3)=funfVI(nn-iline1+1,1:nn-1); %
end
for j=1:nn-1, % LOOP ON THE ETA-LINES OF FACE I. j= order of encountering !!!
    alfaspline(1:nn)=alfa(1:nn,j);
    funspline(1:nn)=funfI(1:nn,j);
    ppspline=spline(alfaspline,funspline);
    funaI1(1:nn)=ppval(ppspline,alfa1(nn+1-[1:nn],j)); % VERIFIER SI C'EST BIEN "NN+1-[1,NN]" SUR LE PLAN LOGIQUE: je ne suis pas sur!!!!!
    vb_fV(1:nn,3*nn-3+j)=funaI1(1:nn);
end 
 % FILTRAGE
% for iline1=1:nn,
%    vb_fV(iline1,1:nb)=vb_fV(iline1,1:nb)*ftr;
% end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCUL DES DERIVEES BETA SUR LE RESEAU DE GRANDS CERCLES V-BETA
betag5=zeros(1,nb);
xb=zeros(1,nb);yb=zeros(1,nb);zb=zeros(1,nb);
funbdex5=zeros(nn,nb);
vbd_fV=zeros(nn,nb);
for iline1=1:nn,
 betag5=betag(iline1,:);
 dbetag5=k*betag5';
 dbetag5(1)=dbetag5(1)+2*pi;
 dbetag5(na)=dbetag5(na)+2*pi;
 betag5d=3*(p\dbetag5);
 %
 funb5=vb_fV(iline1,:);
 funbd5=3*((p\k)*funb5');
 funbd6=funbd5./betag5d;
 vbd_fV(iline1,:)=funbd6; % VALUE OF THE DERIVATIVE WITH RESPECT TO ANGLE BETA.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%% FIN ASSEMBLAGE DES DERIVEES HERMITIENNES SUR LES 6 RESEAUX DE CERCLES %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ETAPE 2 - ASSEMBLAGE DES DERIVEES ALPHA/BETA SUR CHACUNE DES 6 FACES 
% DEDUITES DES DERIVEES HERMITIENNES SUR LES 6 RESEAUX DE GRANDS CERCLES.
% 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ETAPE 3 - ASSEMBLAGE GRADIENTS EN FONCTION DES DERIVEES ALPHA/BETA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8.1 - FACE I: CALCUL DU GRADIENT EN FONCTION DES DERIVEES ALPHA ET BETA
% FACE
alfa_fI=zeros(nn,nn);beta_fI=zeros(nn,nn);
alfa_fI=alfag(1:nn,1:nn);beta_fI=betag(1:nn,1:nn);
%
%grad_I=zeros(nn,nn,3);
% 
for i=1:nn,
    for j=1:nn,
        xwk=(1+xx(i,j)^2)*(1+yy(i,j)^2)/(deltab(i,j)^2);
        xwk_xi  = (1/sqrt(1+yy(i,j)^2))*xwk;
        xwk_eta = (1/sqrt(1+xx(i,j)^2))*xwk;
        grad_I(i,j,1:3)=xwk_xi*dg_alfa(i,j,1)*gxi_I(i,j,1:3) + xwk_eta*dg_beta(i,j,1)*geta_I(i,j,1:3);
    end
end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % 8.2  FACE II: CALCUL DU GRADIENT EN FONCTION DES DERIVEES ALPHA ET BETA FACE II
%% alfa_fII=zeros(nn,nn);beta_fII=zeros(nn,nn);
alfa_fII=zeros(nn,nn);beta_fII=zeros(nn,nn);
alfa_fII=alfag(1:nn,1:nn);beta_fII=betag(1:nn,1:nn);
%
%grad_II=zeros(nn,nn,3);
%
for i=1:nn,
    for j=1:nn,
        xwk=(1+xx(i,j)^2)*(1+yy(i,j)^2)/(deltab(i,j)^2);
        xwk_xi  = (1/sqrt(1+yy(i,j)^2))*xwk;
        xwk_eta = (1/sqrt(1+xx(i,j)^2))*xwk;
        grad_II(i,j,1:3)=xwk_xi*dg_alfa(i,j,2)*gxi_II(i,j,1:3) + xwk_eta*dg_beta(i,j,2)*geta_II(i,j,1:3);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8.3 -  FACE III: CALCUL DU GRADIENT EN FONCTION DES DERIVEES ALPHA ET BETA FACE III
%alfa_fIII=zeros(nn,nn);beta_fIII=zeros(nn,nn);
alfa_fIII=zeros(nn,nn);beta_fIII=zeros(nn,nn);
alfa_fIII=alfag(1:nn,1:nn);beta_fIII=betag(1:nn,1:nn);
%
%grad_III=zeros(nn,nn,3);
%
% FACE III: CALCUL DU GRADIENT EN FONCTION DES DERIVEES ALPHA ET BETA FACE
for i=1:nn,
    for j=1:nn,
        xwk=(1+xx(i,j)^2)*(1+yy(i,j)^2)/(deltab(i,j)^2);
        xwk_xi  = (1/sqrt(1+yy(i,j)^2))*xwk;
        xwk_eta = (1/sqrt(1+xx(i,j)^2))*xwk;
        grad_III(i,j,1:3)=xwk_xi*dg_alfa(i,j,3)*gxi_III(i,j,1:3) - xwk_eta*dg_beta(i,j,3)*geta_III(i,j,1:3);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8.4 - FACE IV : CALCUL DU GRADIENT EN FONCTION DES DERIVEES ALPHA ET BETA
% 
alfa_fIV=zeros(nn,nn);beta_fIV=zeros(nn,nn);
alfa_fIV=alfag(1:nn,1:nn);beta_fIV=betag(1:nn,1:nn);
%
%grad_IV=zeros(nn,nn,3);
%
for i=1:nn,
    for j=1:nn,
        xwk=(1+xx(i,j)^2)*(1+yy(i,j)^2)/(deltab(i,j)^2);
        xwk_xi  = (1/sqrt(1+yy(i,j)^2))*xwk;
        xwk_eta = (1/sqrt(1+xx(i,j)^2))*xwk;
        grad_IV(i,j,1:3)=xwk_xi*dg_alfa(i,j,4)*gxi_IV(i,j,1:3) - xwk_eta*dg_beta(i,j,4)*geta_IV(i,j,1:3);
    end
end
%
% 8.5 - FACE V: CALCUL DU GRADIENT EN FONCTION DES DERIVEES ALPHA ET BETA
% ------------------------------------------------------------------------
alfa_fV=zeros(nn,nn);beta_fV=zeros(nn,nn);
alfa_fV=alfag(1:nn,1:nn);beta_fV=betag(1:nn,1:nn);
%
%grad_V=zeros(nn,nn,3);
% 
% FACE V: CALCUL DU GRADIENT EN FONCTION DES DERIVEES ALPHA ET BETA FACE V
for i=1:nn,
    for j=1:nn,
        xwk=(1+xx(i,j)^2)*(1+yy(i,j)^2)/(deltab(i,j)^2);
        xwk_xi  = (1/sqrt(1+yy(i,j)^2))*xwk;
        xwk_eta = (1/sqrt(1+xx(i,j)^2))*xwk;
        grad_V(i,j,1:3)=xwk_xi*dg_alfa(i,j,5)*gxi_V(i,j,1:3) + xwk_eta*dg_beta(i,j,5)*geta_V(i,j,1:3);
    end
end
% 8.6 - FACE VI : CALCUL DU GRADIENT EN FONCTION DES DERIVEES ALPHA ET BETA 
alfa_fVI=zeros(nn,nn);beta_fVI=zeros(nn,nn);
alfa_fVI=alfag(1:nn,1:nn);beta_fVI=betag(1:nn,1:nn);
%
%grad_VI=zeros(nn,nn,3);
% FACE VI: CALCUL DU GRADIENT EN FONCTION DES DERIVEES ALPHA ET BETA FACE
for i=1:nn,
    for j=1:nn,
        xwk=(1+xx(i,j)^2)*(1+yy(i,j)^2)/(deltab(i,j)^2);
        xwk_xi  = (1/sqrt(1+yy(i,j)^2))*xwk;
        xwk_eta = (1/sqrt(1+xx(i,j)^2))*xwk;
        grad_VI(i,j,1:3)=-xwk_xi*dg_alfa(i,j,6)*gxi_VI(i,j,1:3) + xwk_eta*dg_beta(i,j,6)*geta_VI(i,j,1:3);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1/2 SOMME ARRETES, 1/3 SOMME SOMMETS
% COMPONENT 1
uwk_I=grad_I(1:nn,1:nn,1);uwk_II=grad_II(1:nn,1:nn,1);uwk_III=grad_III(1:nn,1:nn,1);
uwk_IV=grad_IV(1:nn,1:nn,1);uwk_V=grad_V(1:nn,1:nn,1);uwk_VI=grad_VI(1:nn,1:nn,1);
[grad_I(1:nn,1:nn,1),grad_II(1:nn,1:nn,1),grad_III(1:nn,1:nn,1),grad_IV(1:nn,1:nn,1),grad_V(1:nn,1:nn,1),grad_VI(1:nn,1:nn,1)]=...
    ds72(uwk_I,uwk_II,uwk_III,uwk_IV,uwk_V,uwk_VI,n,nn);
% COMPONENT 2
uwk_I=grad_I(1:nn,1:nn,2);uwk_II=grad_II(1:nn,1:nn,2);uwk_III=grad_III(1:nn,1:nn,2);
uwk_IV=grad_IV(1:nn,1:nn,2);uwk_V=grad_V(1:nn,1:nn,2);uwk_VI=grad_VI(1:nn,1:nn,2);
[grad_I(1:nn,1:nn,2),grad_II(1:nn,1:nn,2),grad_III(1:nn,1:nn,2),grad_IV(1:nn,1:nn,2),grad_V(1:nn,1:nn,2),grad_VI(1:nn,1:nn,2)]=...
    ds72(uwk_I,uwk_II,uwk_III,uwk_IV,uwk_V,uwk_VI,n,nn);
% COMPONENT 3
uwk_I=grad_I(1:nn,1:nn,3);uwk_II=grad_II(1:nn,1:nn,3);uwk_III=grad_III(1:nn,1:nn,3);
uwk_IV=grad_IV(1:nn,1:nn,3);uwk_V=grad_V(1:nn,1:nn,3);uwk_VI=grad_VI(1:nn,1:nn,3);
[grad_I(1:nn,1:nn,3),grad_II(1:nn,1:nn,3),grad_III(1:nn,1:nn,3),grad_IV(1:nn,1:nn,3),grad_V(1:nn,1:nn,3),grad_VI(1:nn,1:nn,3)]=...
    ds72(uwk_I,uwk_II,uwk_III,uwk_IV,uwk_V,uwk_VI,n,nn);













