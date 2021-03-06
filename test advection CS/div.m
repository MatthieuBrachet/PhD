function [div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
    div(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI,n,nn)
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
global p1 k1;
global  aaa bbb itestop
%%%%%%%%%%%%%%%%%%%
%
% ETAPES 1 A 5: CALCUL DIVERGENCE FACES I EY III
% ETAPES 6 A 10: CALCUL DIVERGENCE FACES II EY IV
% ETAPES 11 A 15: CALCUL DIVERGENCE FACES V EY VI
% -------------------------------------------------------------------
% ETAPE 1: CALCUL DE LA PARTIE D/DXI DE LA DIVERGENCE (FACE I ET III)
% -------------------------------------------------------------------
% RESEAU I-ALPHA
% --------------
%gxit_Ia=zeros(4*(nn-1),nn,1:3);
% X et Y modifies: I et III idem que X et Y.
% Interpoles sur II et IV.
xxtIa_I=zeros(nn,nn);xxtIa_II=zeros(nn,nn);
xxtIa_III=zeros(nn,nn);xxtIa_IV=zeros(nn,nn);
%
yytIa_I=zeros(nn,nn);yytIa_II=zeros(nn,nn);
yytIa_III=zeros(nn,nn);yytIa_IV=zeros(nn,nn);
%
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
        xwk1=fun5(x_fIV(i,j)); % CORRECTION P67.M
        xwk2=fun7(x_fIV(i,j));
        xxtIa_IV(i,j)=y_fIV(i,j)*xwk1;
        yytIa_IV(i,j)=z_fIV(i,j)*xwk2;
    end
end
%
deltatIa_I=zeros(nn,nn);deltatIa_II=zeros(nn,nn);
deltatIa_III=zeros(nn,nn);deltatIa_IV=zeros(nn,nn);
%
deltatIa_I=(1+xxtIa_I.^2+yytIa_I.^2); 
deltatIa_II=(1+xxtIa_II.^2+yytIa_II.^2); 
deltatIa_III=(1+xxtIa_III.^2+yytIa_III.^2); 
deltatIa_IV=(1+xxtIa_IV.^2+yytIa_IV.^2);
%
deltabtIa_I=zeros(nn,nn);deltabtIa_II=zeros(nn,nn);
deltabtIa_III=zeros(nn,nn);deltabtIa_IV=zeros(nn,nn);
%
deltabtIa_I=sqrt(deltatIa_I);
deltabtIa_II=sqrt(deltatIa_II);
deltabtIa_III=sqrt(deltatIa_III);
deltabtIa_IV=sqrt(deltatIa_IV);
%
% (sqrt(G):
gtIa_I=zeros(nn,nn);gtIa_II=zeros(nn,nn);
gtIa_III=zeros(nn,nn);gtIa_IV=zeros(nn,nn);
%
gtIa_I=(radius^2)*(1+xxtIa_I.^2).*(1+yytIa_I.^2)./(deltabtIa_I.^3);
gtIa_II=(radius^2)*(1+xxtIa_II.^2).*(1+yytIa_II.^2)./(deltabtIa_II.^3);
gtIa_III=(radius^2)*(1+xxtIa_III.^2).*(1+yytIa_III.^2)./(deltabtIa_III.^3);
gtIa_IV=(radius^2)*(1+xxtIa_IV.^2).*(1+yytIa_IV.^2)./(deltabtIa_IV.^3);
%
gxiIa_I=zeros(nn,nn,3);gxiIa_II=zeros(nn,nn,3);
gxiIa_III=zeros(nn,nn,3);gxiIa_IV=zeros(nn,nn,3);
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
      %
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESEAU 1 : ASSEMBLAGE DES DONNEES SUR LE RESEAU I-ALPHA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%---------debut section gr59.m pour derivee hermitienne reseau I-alpha
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
% break;
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
%%%%%%%%% fin section gr59.m------------------------
%
% ETAPE 2: CALCUL DE LA PARTIE D/DETA DE LA DIVERGENCE (FACE I ET III)
% -------------------------------------------------------------------
%
% RESEAU I-BETA
% --------------
%gxit_Ib=zeros(4*(nn-1),nn,1:3);
% X et Y modifies: I et III idem que X et Y.
% Interpoles sur V et VI
xxtIb_I=zeros(nn,nn);xxtIb_V=zeros(nn,nn);
xxtIb_III=zeros(nn,nn);xxtIb_VI=zeros(nn,nn);
%
yytIb_I=zeros(nn,nn);yytIb_V=zeros(nn,nn);
yytIb_III=zeros(nn,nn);yytIb_VI=zeros(nn,nn);
%
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
%
deltatIb_I=zeros(nn,nn);deltatIb_V=zeros(nn,nn);
deltatIb_III=zeros(nn,nn);deltatIb_VI=zeros(nn,nn);
%
deltatIb_I=(1+xxtIb_I.^2+yytIb_I.^2); 
deltatIb_V=(1+xxtIb_V.^2+yytIb_V.^2); 
deltatIb_III=(1+xxtIb_III.^2+yytIb_III.^2); 
deltatIb_VI=(1+xxtIb_VI.^2+yytIb_VI.^2);
%
deltabtIb_I=zeros(nn,nn);deltabtIb_V=zeros(nn,nn);
deltabtIb_III=zeros(nn,nn);deltabtIb_VI=zeros(nn,nn);
%
deltabtIb_I=sqrt(deltatIb_I);
deltabtIb_V=sqrt(deltatIb_V);
deltabtIb_III=sqrt(deltatIb_III);
deltabtIb_VI=sqrt(deltatIb_VI);
% (sqrt(G):
gtIb_I=zeros(nn,nn);gtIb_V=zeros(nn,nn);
gtIb_III=zeros(nn,nn);gtIb_VI=zeros(nn,nn);
% 
gtIb_I=(radius^2)*(1+xxtIb_I.^2).*(1+yytIb_I.^2)./(deltabtIb_I.^3);
gtIb_V=(radius^2)*(1+xxtIb_V.^2).*(1+yytIb_V.^2)./(deltabtIb_V.^3);
gtIb_III=(radius^2)*(1+xxtIb_III.^2).*(1+yytIb_III.^2)./(deltabtIb_III.^3);
gtIb_VI=(radius^2)*(1+xxtIb_VI.^2).*(1+yytIb_VI.^2)./(deltabtIb_VI.^3);
%
getaIb_I=zeros(nn,nn,3);getaIb_V=zeros(nn,nn,3);
getaIb_III=zeros(nn,nn,3);getaIb_VI=zeros(nn,nn,3);
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
      %
      % getaIb_V(i,j,1:3)=getaIb_V(i,j,1:3)/(1+(z_fV(i,j)*xwk1)^2);
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
      %
      %getaIb_VI(i,j,1:3)=getaIb_VI(i,j,1:3)/(1+(z_fVI(i,j)*xwk1)^2);
      getaIb_VI(i,j,1:3)=getaIb_VI(i,j,1:3)/(1+(z_fVI(i,j)*xwk2)^2);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESEAU 2 : ASSEMBLAGE DES DONNEES SUR LE RESEAU I-BETA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%---------debut section gr59.m pour derivee hermitienne reseau I-beta
%
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
%%%%%%%%% fin section gr59.m------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%% FIN ASSEMBLAGE DES DERIVEES HERMITIENNES SUR LES 2 RESEAUX DE CERCLES %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ETAPE 3 - ASSEMBLAGE DES DERIVEES ALPHA/BETA SUR 
% les faces I et III.
% DEDUITES DES DERIVEES HERMITIENNES SUR LES 2 RESEAUX DE GRANDS CERCLES
% reseaux I_alpha et I_beta
% 
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
%%% ASSEMBLAGE DIVERGENCE FACES I ET III
%%%--------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ETAPE 4 - ASSEMBLAGE DIVERGENCE FACES I ET III EN FONCTION DES DERIVEES ALPHA/BETA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8.1 - FACE I: CALCUL DE LA DIVERGENCE EN FONCTION DES DERIVEES ALPHA ET
% BETA FACE I
%
alfa_fI=zeros(nn,nn);beta_fI=zeros(nn,nn);
alfa_fI=alfag(1:nn,1:nn);beta_fI=betag(1:nn,1:nn);
%
div_fI=zeros(nn,nn);
% 
gt_I=zeros(nn,nn);
gt_I=gtIa_I; % gtIa_I=gtIb_I=gt_I
for i=1:nn,
    for j=1:nn,
        xwk=(1+xx(i,j)^2)*(1+yy(i,j)^2)/(deltab(i,j)^2);
        xwk_xi  = (1/sqrt(1+yy(i,j)^2))*xwk;
        xwk_eta = (1/sqrt(1+xx(i,j)^2))*xwk;
        div_fI(i,j)=(xwk_xi*dg_alfa(i,j,1) + xwk_eta*dg_beta(i,j,1))/gt_I(i,j);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8.3 -  FACE III: CALCUL DE LA DIVERGENCE EN FONCTION DES DERIVEES ALPHA
% ET BETA FACE III
alfa_fIII=zeros(nn,nn);beta_fIII=zeros(nn,nn);
alfa_fIII=alfag(1:nn,1:nn);beta_fIII=betag(1:nn,1:nn);
%
div_fIII=zeros(nn,nn);
%
gt_III=zeros(nn,nn);
gt_III=gtIa_III; % gtIa_III=gtIb_III=gt_III
for i=1:nn,
    for j=1:nn,
        xwk=(1+xx(i,j)^2)*(1+yy(i,j)^2)/(deltab(i,j)^2);
        xwk_xi  = (1/sqrt(1+yy(i,j)^2))*xwk;
        xwk_eta = (1/sqrt(1+xx(i,j)^2))*xwk;
        div_fIII(i,j)=(xwk_xi*dg_alfa(i,j,3) - xwk_eta*dg_beta(i,j,3))/gt_III(i,j);
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ETAPE 5: COMPARAISON DIVERGENCE EXACTE ET CALCULEE POUR FACES I ET III.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A - FACE I: 
% COMPARAISON DIVERGENCE EXACTE / DIVERGENCE CALCULEE SUR FACE I
% gradient exact
% [func,func_x,func_y,func_z]=fun3(x_fI,y_fI,z_fI);
% projection du gradient sur la sphere
% norm_I=zeros(3,1);
% gradse_I=zeros(nn,nn,3);gradset_I=zeros(nn,nn,3); % GRADSET_I=GRADIENT TANGENTIEL= SPHERIQUE
% gradse_I(1:nn,1:nn,1)=func_x ; 
% gradse_I(1:nn,1:nn,2)=func_y ; 
% gradse_I(1:nn,1:nn,3)=func_z;
% norm_xI=zeros(nn,nn);norm_yI=zeros(nn,nn);norm_zI=zeros(nn,nn);
% for i=1:nn,
%     for j=1:nn,
%         norm_xI(i,j)=x_fI(i,j) ; norm_yI(i,j)=y_fI(i,j) ; norm_zI(i,j)=z_fI(i,j);
%         xwk(1)=gradse_I(i,j,1) ; xwk(2)=gradse_I(i,j,2) ; xwk(3)=gradse_I(i,j,3);
%         gradset_I(i,j,1)=gradse_I(i,j,1)-(xwk(1)*norm_xI(i,j)+xwk(2)*norm_yI(i,j)+xwk(3)*norm_zI(i,j))*norm_xI(i,j);
%         gradset_I(i,j,2)=gradse_I(i,j,2)-(xwk(1)*norm_xI(i,j)+xwk(2)*norm_yI(i,j)+xwk(3)*norm_zI(i,j))*norm_yI(i,j);
%         gradset_I(i,j,3)=gradse_I(i,j,3)-(xwk(1)*norm_xI(i,j)+xwk(2)*norm_yI(i,j)+xwk(3)*norm_zI(i,j))*norm_zI(i,j);
%     end
% end
% calcul erreur sur le gradient en X:
% errg_fI=zeros(nn,nn);errg_fIII=zeros(nn,nn);
% for i=1:nn,
%     for j=1:nn,
%         errg_fI(i,j)=max(abs(div_fI(i,j)-dfunI(i,j)));
%     end
% end
% % for i=1:nn,
% %     for j=1:nn,
% %         errg_fII(i,j)=max(abs(grad_II(i,j,:)-gradset_II(i,j,:)));
% %     end
% % end
% % B - FACE III: 
% % COMPARAISON DIVERGENCE EXACTE / DIVERGENCE CALCULEE SUR FACE III
% for i=1:nn,
%     for j=1:nn,
%         errg_fIII(i,j)=max(abs(div_fIII(i,j)-dfunIII(i,j)));
%     end
% end
% % for i=1:nn,
% %     for j=1:nn,
% %         errg_fIV(i,j)=max(abs(grad_IV(i,j,:)-gradset_IV(i,j,:)));
% %     end
% % end
% % for i=1:nn,
% %     for j=1:nn,
% %         errg_fV(i,j)=max(abs(grad_V(i,j,:)-gradset_V(i,j,:)));
% %     end
% % end
% % for i=1:nn,
% %     for j=1:nn,
% %         errg_fVI(i,j)=max(abs(grad_VI(i,j,:)-gradset_VI(i,j,:)));
% %     end
% % end
% errI_25=max(max(abs(errg_fI)));
% %errII_26=max(max(abs(errg_fII)))
% errIII_27=max(max(abs(errg_fIII)));
% % errIV_28=max(max(abs(errg_fIV)))
% % errV_29=max(max(abs(errg_fV)))
% % errVI_30=max(max(abs(errg_fVI)))
% %err_grad_31=max([errVI_30,errV_29,errIV_28,errIII_27,errII_26,errI_25])
% err_grad_31=max([errIII_27,errI_25]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%break;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ETAPE 6: CALCUL DE LA PARTIE D/DXI DE LA DIVERGENCE (FACE II ET IV)
% -------------------------------------------------------------------
% RESEAU II-ALPHA
% --------------
%gxit_Ia=zeros(4*(nn-1),nn,1:3);
% X et Y modifies: II et IV idem que X et Y.
% Interpoles sur I et III.
xxtIIa_II=zeros(nn,nn);xxtIIa_III=zeros(nn,nn);
xxtIIa_IV=zeros(nn,nn);xxtIIa_I=zeros(nn,nn);
%
yytIIa_II=zeros(nn,nn);yytIIa_III=zeros(nn,nn);
yytIIa_IV=zeros(nn,nn);yytIIa_I=zeros(nn,nn);
%
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
%
deltatIIa_II=zeros(nn,nn);deltatIIa_III=zeros(nn,nn);
deltatIIa_IV=zeros(nn,nn);deltatIIa_I=zeros(nn,nn);
%
deltatIIa_II=(1+xxtIIa_II.^2+yytIIa_II.^2); 
deltatIIa_III=(1+xxtIIa_III.^2+yytIIa_III.^2); 
deltatIIa_IV=(1+xxtIIa_IV.^2+yytIIa_IV.^2); 
deltatIIa_I=(1+xxtIIa_I.^2+yytIIa_I.^2);
%
deltabtIIa_II=zeros(nn,nn);deltabtIIa_III=zeros(nn,nn);
deltabtIIa_IV=zeros(nn,nn);deltabtIIa_I=zeros(nn,nn);
%
deltabtIIa_II=sqrt(deltatIIa_II);
deltabtIIa_III=sqrt(deltatIIa_III);
deltabtIIa_IV=sqrt(deltatIIa_IV);
deltabtIIa_I=sqrt(deltatIIa_I);
%
% (sqrt(G):
gtIIa_II=zeros(nn,nn);gtIIa_III=zeros(nn,nn);
gtIIa_IV=zeros(nn,nn);gtIIa_I=zeros(nn,nn);
%
gtIIa_II=(radius^2)*(1+xxtIIa_II.^2).*(1+yytIIa_II.^2)./(deltabtIIa_II.^3);
gtIIa_III=(radius^2)*(1+xxtIIa_III.^2).*(1+yytIIa_III.^2)./(deltabtIIa_III.^3);
gtIIa_IV=(radius^2)*(1+xxtIIa_IV.^2).*(1+yytIIa_IV.^2)./(deltabtIIa_IV.^3);
gtIIa_I=(radius^2)*(1+xxtIIa_I.^2).*(1+yytIIa_I.^2)./(deltabtIIa_I.^3);
%
gxiIIa_II=zeros(nn,nn,3);gxiIIa_III=zeros(nn,nn,3);
gxiIIa_IV=zeros(nn,nn,3);gxiIIa_I=zeros(nn,nn,3);
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
      %
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
      %
      gxiIIa_I(i,j,1:3)=gxiIIa_I(i,j,1:3)/(1+(x_fI(i,j)*xwk2)^2);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESEAU 2 : ASSEMBLAGE DES DONNEES SUR LE RESEAU II-ALPHA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
funfI=zeros(nn,nn);funfII=zeros(nn,nn);funfIII=zeros(nn,nn);
funfIV=zeros(nn,nn);funfV=zeros(nn,nn);funfVI=zeros(nn,nn);
%
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
%---------debut section gr59.m pour derivee hermitienne reseau II-alpha
%
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
%%%%%%%%% fin section gr59.m------------------------
% break;
%
% ETAPE 7: CALCUL DE LA PARTIE D/DETA DE LA DIVERGENCE (FACE II ET IV)
% -------------------------------------------------------------------
%
% RESEAU II-BETA
% --------------
%gxit_Ib=zeros(4*(nn-1),nn,1:3);
% X et Y modifies: I et III idem que X et Y.
% Interpoles sur V et VI
xxtIIb_II=zeros(nn,nn);xxtIIb_V=zeros(nn,nn);
xxtIIb_IV=zeros(nn,nn);xxtIIb_VI=zeros(nn,nn);
%
yytIIb_II=zeros(nn,nn);yytIIb_V=zeros(nn,nn);
yytIIb_IV=zeros(nn,nn);yytIIb_VI=zeros(nn,nn);
%
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
%
deltatIIb_II=zeros(nn,nn);deltatIIb_V=zeros(nn,nn);
deltatIIb_IV=zeros(nn,nn);deltatIIb_VI=zeros(nn,nn);
%
deltatIIb_II=(1+xxtIIb_II.^2+yytIIb_II.^2); 
deltatIIb_V=(1+xxtIIb_V.^2+yytIIb_V.^2); 
deltatIIb_IV=(1+xxtIIb_IV.^2+yytIIb_IV.^2); 
deltatIIb_VI=(1+xxtIIb_VI.^2+yytIIb_VI.^2);
%
deltabtIIb_II=zeros(nn,nn);deltabtIIb_V=zeros(nn,nn);
deltabtIIb_IV=zeros(nn,nn);deltabtIIb_VI=zeros(nn,nn);
%
deltabtIIb_II=sqrt(deltatIIb_II);
deltabtIIb_V=sqrt(deltatIIb_V);
deltabtIIb_IV=sqrt(deltatIIb_IV);
deltabtIIb_VI=sqrt(deltatIIb_VI);
% (sqrt(G):
gtIIb_II=zeros(nn,nn);gtIIb_V=zeros(nn,nn);
gtIIb_IV=zeros(nn,nn);gtIIb_VI=zeros(nn,nn);
% 
gtIIb_II=(radius^2)*(1+xxtIIb_II.^2).*(1+yytIIb_II.^2)./(deltabtIIb_II.^3);
gtIIb_V=(radius^2)*(1+xxtIIb_V.^2).*(1+yytIIb_V.^2)./(deltabtIIb_V.^3);
gtIIb_IV=(radius^2)*(1+xxtIIb_IV.^2).*(1+yytIIb_IV.^2)./(deltabtIIb_IV.^3);
gtIIb_VI=(radius^2)*(1+xxtIIb_VI.^2).*(1+yytIIb_VI.^2)./(deltabtIIb_VI.^3);
%
getaIIb_II=zeros(nn,nn,3);getaIIb_V=zeros(nn,nn,3);
getaIIb_IV=zeros(nn,nn,3);getaIIb_VI=zeros(nn,nn,3);
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
      %
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
      getaIIb_VI(i,j,2)=-z_fVI(i,j)*xwk1*xwk2; ;
      getaIIb_VI(i,j,3)= xwk2;
      %
      getaIIb_VI(i,j,1:3)=getaIIb_VI(i,j,1:3)/(1+(z_fVI(i,j)*xwk2)^2);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESEAU 2 : ASSEMBLAGE DES DONNEES SUR LE RESEAU II-BETA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%
%---------debut section gr59.m pour derivee hermitienne reseau II-beta
%
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
% FACE IV: TRANSFERT OF DATA OF FACE IV
for iline1=1:nn,
    vb_fII(iline1,2*nn-1:3*nn-3)=funfIV(iline1,nn:-1:2); % symetrie sur adresses en i + inversion en j !
end
% FACE VI: SPLINE INTERPOLATION A L'AIDE DES ANGLES BETA
% INVERSION OF THE LOOP ON THE ISO-ETA LINES AND ON THE XI OF FACE VI
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
%%%%%%%%% fin section gr59.m------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%% FIN ASSEMBLAGE DES DERIVEES HERMITIENNES SUR LES 2 RESEAUX DE CERCLES %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ETAPE 8 - ASSEMBLAGE DES DERIVEES ALPHA/BETA SUR 
% les faces I et III.
% DEDUITES DES DERIVEES HERMITIENNES SUR LES 2 RESEAUX DE GRANDS CERCLES
% reseaux I_alpha et I_beta
% 
%% dg_alfa=zeros(nn,nn,6);dg_beta=zeros(nn,nn,6);

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
%%% ASSEMBLAGE DIVERGENCE FACES II ET IV
%%%--------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ETAPE 9 - ASSEMBLAGE DIVERGENCE FACES II ET IV EN FONCTION DES DERIVEES ALPHA/BETA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8.1 - FACE II: CALCUL DE LA DIVERGENCE EN FONCTION DES DERIVEES ALPHA ET
% BETA FACE II
%
alfa_fII=zeros(nn,nn);beta_fII=zeros(nn,nn);
alfa_fII=alfag(1:nn,1:nn);beta_fII=betag(1:nn,1:nn);
%
div_fII=zeros(nn,nn);
% 
gt_II=zeros(nn,nn);
gt_II=gtIIa_II; % gtIIa_II=gtIIb_II=gt_II
for i=1:nn,
    for j=1:nn,
        xwk=(1+xx(i,j)^2)*(1+yy(i,j)^2)/(deltab(i,j)^2);
        xwk_xi  = (1/sqrt(1+yy(i,j)^2))*xwk;
        xwk_eta = (1/sqrt(1+xx(i,j)^2))*xwk;
        div_fII(i,j)=(xwk_xi*dg_alfa(i,j,2) + xwk_eta*dg_beta(i,j,2))/gt_II(i,j);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8.3 -  FACE IV: CALCUL DE LA DIVERGENCE EN FONCTION DES DERIVEES ALPHA
% ET BETA FACE IV
alfa_fIV=zeros(nn,nn);beta_fIV=zeros(nn,nn);
alfa_fIV=alfag(1:nn,1:nn);beta_fIV=betag(1:nn,1:nn);
%
div_fIV=zeros(nn,nn);
%
gt_IV=zeros(nn,nn);
gt_IV=gtIIa_IV; % gtIIa_IV=gtIIb_IV=gt_IV
for i=1:nn,
    for j=1:nn,
        xwk=(1+xx(i,j)^2)*(1+yy(i,j)^2)/(deltab(i,j)^2);
        xwk_xi  = (1/sqrt(1+yy(i,j)^2))*xwk;
        xwk_eta = (1/sqrt(1+xx(i,j)^2))*xwk;
        div_fIV(i,j)=(xwk_xi*dg_alfa(i,j,4) - xwk_eta*dg_beta(i,j,4))/gt_IV(i,j);
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ETAPE 10: COMPARAISON DIVERGENCE EXACTE ET CALCULEE POUR FACES II ET IV.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A - FACE I: 
% COMPARAISON DIVERGENCE EXACTE / DIVERGENCE CALCULEE SUR FACE I
% gradient exact
% [func,func_x,func_y,func_z]=fun3(x_fI,y_fI,z_fI);
% projection du gradient sur la sphere
% norm_I=zeros(3,1);
% gradse_I=zeros(nn,nn,3);gradset_I=zeros(nn,nn,3); % GRADSET_I=GRADIENT TANGENTIEL= SPHERIQUE
% gradse_I(1:nn,1:nn,1)=func_x ; 
% gradse_I(1:nn,1:nn,2)=func_y ; 
% gradse_I(1:nn,1:nn,3)=func_z;
% norm_xI=zeros(nn,nn);norm_yI=zeros(nn,nn);norm_zI=zeros(nn,nn);
% for i=1:nn,
%     for j=1:nn,
%         norm_xI(i,j)=x_fI(i,j) ; norm_yI(i,j)=y_fI(i,j) ; norm_zI(i,j)=z_fI(i,j);
%         xwk(1)=gradse_I(i,j,1) ; xwk(2)=gradse_I(i,j,2) ; xwk(3)=gradse_I(i,j,3);
%         gradset_I(i,j,1)=gradse_I(i,j,1)-(xwk(1)*norm_xI(i,j)+xwk(2)*norm_yI(i,j)+xwk(3)*norm_zI(i,j))*norm_xI(i,j);
%         gradset_I(i,j,2)=gradse_I(i,j,2)-(xwk(1)*norm_xI(i,j)+xwk(2)*norm_yI(i,j)+xwk(3)*norm_zI(i,j))*norm_yI(i,j);
%         gradset_I(i,j,3)=gradse_I(i,j,3)-(xwk(1)*norm_xI(i,j)+xwk(2)*norm_yI(i,j)+xwk(3)*norm_zI(i,j))*norm_zI(i,j);
%     end
% end
% calcul erreur sur le gradient en X:
% errg_fII=zeros(nn,nn);errg_fIV=zeros(nn,nn);
% for i=1:nn,
%     for j=1:nn,
%         errg_fII(i,j)=max(abs(div_fII(i,j)-dfunII(i,j)));
%     end
% end
% % B - FACE IV: 
% % COMPARAISON DIVERGENCE EXACTE / DIVERGENCE CALCULEE SUR FACE IV
% for i=1:nn,
%     for j=1:nn,
%         errg_fIV(i,j)=max(abs(div_fIV(i,j)-dfunIV(i,j)));
%     end
% end
% %errI_25=max(max(abs(errg_fI)))
% errII_26=max(max(abs(errg_fII)));
% %errIII_27=max(max(abs(errg_fIII)))
% errIV_28=max(max(abs(errg_fIV)));
% % errV_29=max(max(abs(errg_fV)))
% % errVI_30=max(max(abs(errg_fVI)))
% %err_grad_31=max([errVI_30,errV_29,errIV_28,errIII_27,errII_26,errI_25])
% err_grad_32=max([errIV_28,errII_26]);
% break;
% -------------------------------------------------------------------
% ETAPE 11: CALCUL DE LA PARTIE D/DXI DE LA DIVERGENCE (FACE V ET VI)
% -------------------------------------------------------------------
% RESEAU V-ALPHA
% --------------
%gxit_Ia=zeros(4*(nn-1),nn,1:3);
% X et Y modifies: I et III idem que X et Y.
% Interpoles sur II et IV.
xxtVa_V=zeros(nn,nn);xxtVa_II=zeros(nn,nn);
xxtVa_VI=zeros(nn,nn);xxtVa_IV=zeros(nn,nn);
%
yytVa_V=zeros(nn,nn);yytVa_II=zeros(nn,nn);
yytVa_VI=zeros(nn,nn);yytVa_IV=zeros(nn,nn);
%
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
%
deltatVa_V=zeros(nn,nn);deltatVa_II=zeros(nn,nn);
deltatVa_VI=zeros(nn,nn);deltatVa_IV=zeros(nn,nn);
%
deltatVa_V=(1+xxtVa_V.^2+yytVa_V.^2); 
deltatVa_II=(1+xxtVa_II.^2+yytVa_II.^2); 
deltatVa_VI=(1+xxtVa_VI.^2+yytVa_VI.^2); 
deltatVa_IV=(1+xxtVa_IV.^2+yytVa_IV.^2);
%
deltabtVa_V=zeros(nn,nn);deltabtVa_II=zeros(nn,nn);
deltabtVa_VI=zeros(nn,nn);deltabtVa_IV=zeros(nn,nn);
%
deltabtVa_V=sqrt(deltatVa_V);
deltabtVa_II=sqrt(deltatVa_II);
deltabtVa_VI=sqrt(deltatVa_VI);
deltabtVa_IV=sqrt(deltatVa_IV);
%
% (sqrt(G):
gtVa_V=zeros(nn,nn);gtVa_II=zeros(nn,nn);
gtVa_VI=zeros(nn,nn);gtVa_IV=zeros(nn,nn);
%
gtVa_V =(radius^2)*(1+xxtVa_V.^2).*(1+yytVa_V.^2)./(deltabtVa_V.^3);
gtVa_II=(radius^2)*(1+xxtVa_II.^2).*(1+yytVa_II.^2)./(deltabtVa_II.^3);
gtVa_VI=(radius^2)*(1+xxtVa_VI.^2).*(1+yytVa_VI.^2)./(deltabtVa_VI.^3);
gtVa_IV=(radius^2)*(1+xxtVa_IV.^2).*(1+yytVa_IV.^2)./(deltabtVa_IV.^3);
%
gxiVa_V=zeros(nn,nn,3);gxiVa_II=zeros(nn,nn,3);
gxiVa_VI=zeros(nn,nn,3);gxiVa_IV=zeros(nn,nn,3);
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
      %
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
      %
      gxiVa_IV(i,j,1:3)=gxiVa_IV(i,j,1:3)/(1+(y_fIV(i,j)*xwk2)^2);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESEAU 5 : ASSEMBLAGE DES DONNEES SUR LE RESEAU V-ALPHA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
funfI=zeros(nn,nn);funfII=zeros(nn,nn);funfIII=zeros(nn,nn);
funfIV=zeros(nn,nn);funfV=zeros(nn,nn);funfVI=zeros(nn,nn);
%
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
%
%---------debut section gr59.m pour derivee hermitienne reseau V-alpha
%
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
%%%%%%%%% fin section gr59.m------------------------
%
% ETAPE 12: CALCUL DE LA PARTIE D/DETA DE LA DIVERGENCE (FACE V ET VI)
% -------------------------------------------------------------------
%
% RESEAU V-BETA
% --------------
%gxit_Ib=zeros(4*(nn-1),nn,1:3);
% X et Y modifies: I et III idem que X et Y.
% Interpoles sur V et VI
xxtVb_V=zeros(nn,nn);xxtVb_III=zeros(nn,nn);
xxtVb_VI=zeros(nn,nn);xxtVb_I=zeros(nn,nn);
%
yytVb_V=zeros(nn,nn);yytVb_III=zeros(nn,nn);
yytVb_VI=zeros(nn,nn);yytVb_I=zeros(nn,nn);
%
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
%
deltatVb_V=zeros(nn,nn);deltatVb_III=zeros(nn,nn);
deltatVb_VI=zeros(nn,nn);deltatVb_I=zeros(nn,nn);
%
deltatVb_V=(1+xxtVb_V.^2+yytVb_V.^2); 
deltatVb_III=(1+xxtVb_III.^2+yytVb_III.^2); 
deltatVb_VI=(1+xxtVb_VI.^2+yytVb_VI.^2); 
deltatVb_I=(1+xxtVb_I.^2+yytVb_I.^2);
%
deltabtVb_V=zeros(nn,nn);deltabtVb_III=zeros(nn,nn);
deltabtVb_VI=zeros(nn,nn);deltabtVb_IV=zeros(nn,nn);
%
deltabtVb_V=sqrt(deltatVb_V);
deltabtVb_III=sqrt(deltatVb_III);
deltabtVb_VI=sqrt(deltatVb_VI);
deltabtVb_I=sqrt(deltatVb_I);
% (sqrt(G):
gtVb_V=zeros(nn,nn);gtVb_III=zeros(nn,nn);
gtVb_VI=zeros(nn,nn);gtVb_I=zeros(nn,nn);
% 
gtVb_V=(radius^2)*(1+xxtVb_V.^2).*(1+yytVb_V.^2)./(deltabtVb_V.^3);
gtVb_III=(radius^2)*(1+xxtVb_III.^2).*(1+yytVb_III.^2)./(deltabtVb_III.^3);
gtVb_VI=(radius^2)*(1+xxtVb_VI.^2).*(1+yytVb_VI.^2)./(deltabtVb_VI.^3);
gtVb_I=(radius^2)*(1+xxtVb_I.^2).*(1+yytVb_I.^2)./(deltabtVb_I.^3);
%
getaVb_V=zeros(nn,nn,3);getaVb_III=zeros(nn,nn,3);
getaVb_VI=zeros(nn,nn,3);getaVb_I=zeros(nn,nn,3);
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
      %
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESEAU 6 : ASSEMBLAGE DES DONNEES SUR LE RESEAU V-BETA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%
%---------debut section gr59.m pour derivee hermitienne reseau V-beta
%
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
%%%%%%%%% fin section gr59.m------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%% FIN ASSEMBLAGE DES DERIVEES HERMITIENNES SUR LES 2 RESEAUX DE CERCLES %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ETAPE 13 - ASSEMBLAGE DES DERIVEES ALPHA/BETA SUR 
% les faces I et III.
% DEDUITES DES DERIVEES HERMITIENNES SUR LES 2 RESEAUX DE GRANDS CERCLES
% reseaux I_alpha et I_beta
% 
%dg_alfa=zeros(nn,nn,6);dg_beta=zeros(nn,nn,6);
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
%%% ASSEMBLAGE DIVERGENCE FACES V ET VI
%%%--------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ETAPE 14 - ASSEMBLAGE DIVERGENCE FACES V ET VI EN FONCTION DES DERIVEES ALPHA/BETA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8.1 - FACE V: CALCUL DE LA DIVERGENCE EN FONCTION DES DERIVEES ALPHA ET BETA FACE V
%
alfa_fV=zeros(nn,nn);beta_fV=zeros(nn,nn);
alfa_fV=alfag(1:nn,1:nn);beta_fV=betag(1:nn,1:nn);
%
div_fV=zeros(nn,nn);
% 
gt_V=zeros(nn,nn);
gt_V=gtVb_V; % gtVb_V=gtVb_V=gt_V
for i=1:nn,
    for j=1:nn,
        xwk=(1+xx(i,j)^2)*(1+yy(i,j)^2)/(deltab(i,j)^2);
        xwk_xi  = (1/sqrt(1+yy(i,j)^2))*xwk;
        xwk_eta = (1/sqrt(1+xx(i,j)^2))*xwk;
        div_fV(i,j)=(xwk_xi*dg_alfa(i,j,5) + xwk_eta*dg_beta(i,j,5))/gt_V(i,j);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8.3 -  FACE VI: CALCUL DE LA DIVERGENCE EN FONCTION DES DERIVEES ALPHA
% ET BETA FACE VI
alfa_fVI=zeros(nn,nn);beta_fVI=zeros(nn,nn);
alfa_fVI=alfag(1:nn,1:nn);beta_fVI=betag(1:nn,1:nn);
%
div_fVI=zeros(nn,nn);
%
gt_VI=zeros(nn,nn);
gt_VI=gtVa_VI; % gtVa_VI=gtVb_VI=gt_VI
for i=1:nn,
    for j=1:nn,
        xwk=(1+xx(i,j)^2)*(1+yy(i,j)^2)/(deltab(i,j)^2);
        xwk_xi  = (1/sqrt(1+yy(i,j)^2))*xwk;
        xwk_eta = (1/sqrt(1+xx(i,j)^2))*xwk;
        div_fVI(i,j)=(-xwk_xi*dg_alfa(i,j,6) + xwk_eta*dg_beta(i,j,6))/gt_VI(i,j);
    end
end
%
% demi-somme + 1/3 somme de la divergence
uwk_I=div_fI(1:nn,1:nn,1);uwk_II=div_fII(1:nn,1:nn,1);uwk_III=div_fIII(1:nn,1:nn,1);
uwk_IV=div_fIV(1:nn,1:nn,1);uwk_V=div_fV(1:nn,1:nn,1);uwk_VI=div_fVI(1:nn,1:nn,1);
[div_fI(1:nn,1:nn),div_fII(1:nn,1:nn),div_fIII(1:nn,1:nn),div_fIV(1:nn,1:nn),div_fV(1:nn,1:nn),div_fVI(1:nn,1:nn,1)]=...
    ds_1b(uwk_I,uwk_II,uwk_III,uwk_IV,uwk_V,uwk_VI,n,nn);
