function [funftI,funftII,funftIII,funftIV,funftV,funftVI]=...
    ftr72(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
betaspline=zeros(nn,1);
alfaspline=zeros(nn,1);
funspline=zeros(nn,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
funftI(1:nn,1:nn)=funfI(1:nn,1:nn);
funftII(1:nn,1:nn)=funfII(1:nn,1:nn);
funftIII(1:nn,1:nn)=funfIII(1:nn,1:nn);
funftIV(1:nn,1:nn)=funfIV(1:nn,1:nn);
funftV(1:nn,1:nn)=funfV(1:nn,1:nn);
funftVI(1:nn,1:nn)=funfVI(1:nn,1:nn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESEAU 1 : ASSEMBLAGE DES DONNEES SUR LE RESEAU I-ALPHA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
va_fI=zeros(4*(nn-1),nn); % VALEURS =DONNEES RESEAU FACE I SELON ALPHA
funbII1=zeros(nn,1);
funbIV1=zeros(nn,1);
% Face I : TRANSFERT OF DATA OF FACE I 
for jline1=1:nn, % boucle sur les lignes iso-eta du reseaux Ia
    va_fI(1:nn-1,jline1)=funftI(1:nn-1,jline1);
end
% FACE II: SPLINE INTERPOLATION A L'AIDE DES ANGLES BETA
% INVERSION OF THE LOOP ON THE ISO-ETA LINES AND ON THE XI OF FACE II
for i=1:nn-1 % LOOP ON THE XI OF FACE II
    betaspline(1:nn)=beta(i,1:nn);
    funspline(1:nn)=funftII(i,1:nn);
    ppspline=spline(betaspline,funspline);
    funbII1(1:nn)=ppval(ppspline,betacr(i,1:nn));
    va_fI(nn-1+i,1:nn)=funbII1(1:nn);
end
for jline1=1:nn,
    % FACE III: TRANSFERT OF DATA OF FACE III
    va_fI(2*nn-1:3*nn-3,jline1)=funftIII(1:nn-1,nn-jline1+1); % symetrie sur adresses en j !
end
% FACE IV: SPLINE INTERPOLATION A L'AIDE DES ANGLES BETA
% INVERSION OF THE LOOP ON THE ISO-ETA LINES AND ON THE XI OF FACE IV
for i=1:nn-1, % LOOP ON THE XI OF FACE IV
   betaspline(1:nn)=beta(i,1:nn);
   funspline(1:nn)=funftIV(i,1:nn);
   ppspline=spline(betaspline,funspline);
   funbIV1(1:nn)=ppval(ppspline,betacr(i,nn+1-[1:nn]));
   va_fI(3*nn-3+i,1:nn)=funbIV1(1:nn);
end
% FILTRAGE
for jline1=1:nn,
   va_fI(1:na,jline1)=ftr*va_fI(1:na,jline1);
end
% TRANSFERT FACE I+III
for jline1=1:nn,
   funftI(1:nn,jline1)   = va_fI(1:nn,jline1);
   funftIII(1:nn,nn-jline1+1)    = va_fI(2*nn-1:3*nn-2,jline1);
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
    vb_fI(iline1,1:nn-1)=funftI(iline1,1:nn-1);
end
% FACE V: SPLINE INTERPOLATION A L'AIDE DES ANGLES ALFA
% INVERSION OF THE LOOP ON THE ISO-XI LINES AND ON THE ETA OF FACE V
for j=1:nn-1 % LOOP ON THE ETA OF FACE V
    alfaspline(1:nn)=alfa(1:nn,j);
    funspline(1:nn)=funftV(1:nn,j);
    ppspline=spline(alfaspline,funspline);
    funaV1(1:nn)=ppval(ppspline,alfa1(1:nn,j));
    vb_fI(1:nn,nn-1+j)=funaV1(1:nn);
end
    % FACE III: TRANSFERT OF DATA
    % ---------------------------
 for iline1=1:nn,
     vb_fI(iline1,2*nn-1:3*nn-3)=funftIII(iline1,nn:-1:2); % symetrie sur adresses en i + inversion en j !
 end
     % FACE VI: SPLINE INTERPOLATION A L'AIDE DES ANGLES ALFA
    % --------------------------------------------------------
    % LOADING DATA FOR THE SPLINE INTERPOLATION :
 for j=1:nn-1, % LOOP ON THE ETA OF FACE VI
    alfaspline(1:nn)=alfa(1:nn,j);
    funspline(1:nn)=funftVI(1:nn,j);
    ppspline=spline(alfaspline,funspline);
    funaVI1(1:nn)=ppval(ppspline,alfa1(nn+1-[1:nn],j));
    vb_fI(1:nn,3*nn-3+j)=funaVI1(1:nn);
 end
 % FILTRAGE
for iline1=1:nn,
   vb_fI(iline1,1:nb)=vb_fI(iline1,1:nb)*ftr;
end
% TRANSFERT FACE I+III
for iline1=1:nn,
   funftI(iline1,1:nn)   = vb_fI(iline1,1:nn);
   funftIII(iline1,nn:-1:1)    = vb_fI(iline1,2*nn-1:3*nn-2);
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
    va_fII(1:nn-1,jline1)=funftII(1:nn-1,jline1);
end
% FACE III: SPLINE INTERPOLATION A L'AIDE DES ANGLES BETA
% INVERSION OF THE LOOP ON THE ISO-ETA LINES AND ON THE XI OF FACE II
for i=1:nn-1 % LOOP ON THE XI OF FACE II
    betaspline(1:nn)=beta(i,1:nn);
    funspline(1:nn)=funftIII(i,1:nn);
    ppspline=spline(betaspline,funspline);
    funbIII1(1:nn)=ppval(ppspline,betacr(i,1:nn));
    va_fII(nn-1+i,1:nn)=funbIII1(1:nn);
end
% FACE IV: TRANSFERT OF DATA
    % ---------------------------
for jline1=1:nn,
    va_fII(2*nn-1:3*nn-3,jline1)=funftIV(1:nn-1,nn-jline1+1); % symetrie sur adresses en j !    
end
% % FACE I: SPLINE INTERPOLATION A L'AIDE DES ANGLES BETA
% INVERSION OF THE LOOP ON THE ISO-ETA LINES AND ON THE XI OF FACE I
for i=1:nn-1, % LOOP ON THE XI OF FACE I
   betaspline(1:nn)=beta(i,1:nn);
   funspline(1:nn)=funftI(i,1:nn);
   ppspline=spline(betaspline,funspline);
   funbI1(1:nn)=ppval(ppspline,betacr(i,nn+1-[1:nn]));
   va_fII(3*nn-3+i,1:nn)=funbI1(1:nn);
end
% FILTRAGE
for jline1=1:nn,
   va_fII(1:na,jline1)=ftr*va_fII(1:na,jline1);
end
% TRANSFERT FACE II+IV
for jline1=1:nn,
   funftII(1:nn,jline1)   = va_fII(1:nn,jline1);
   funftIV(1:nn,nn-jline1+1)    = va_fII(2*nn-1:3*nn-2,jline1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESEAU 4: ASSEMBLAGE DES DONNEES SUR LE RESEAU II-BETA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tic;
vb_fII=zeros(nn,4*(nn-1)); % VALEURS =DONNEES RESEAU FACE II SELON BETA
funbV1=zeros(nn,1);
funbVI1=zeros(nn,1);
% Face II : TRANSFERT OF DATA OF FACE II
for iline1=1:nn, % boucle sur les lignes iso-xi du reseau II-beta
    vb_fII(iline1,1:nn-1)=funftII(iline1,1:nn-1);
end
% FACE V: SPLINE INTERPOLATION A L'AIDE DES ANGLES BETA
% INVERSION OF THE LOOP ON THE ISO-ETA LINES AND ON THE XI OF FACE V
for i=1:nn-1 % LOOP ON THE XI OF FACE V. i= order of encountering !!!
    betaspline(1:nn)=beta(nn-i+1,1:nn);
    funspline(1:nn)=funftV(nn-i+1,1:nn);
    ppspline=spline(betaspline,funspline);
    funbV1(1:nn)=ppval(ppspline,betacr(i,1:nn));
    vb_fII(1:nn,nn-1+i)=funbV1(1:nn);
end
% FACE III: TRANSFERT OF DATA OF FACE III
for iline1=1:nn,
    vb_fII(iline1,2*nn-1:3*nn-3)=funftIV(iline1,nn:-1:2); % symetrie sur adresses en i + inversion en j !
end
% FACE VI: SPLINE INTERPOLATION A L'AIDE DES ANGLES BETA
% INVERSION OF THE LOOP ON THE ISO-ETA LINES AND ON THE XI OF FACE Vi
for i=1:nn-1 % LOOP ON THE XI OF FACE VI. i= order of encountering !!!
    betaspline(1:nn)=beta(i,1:nn);
    funspline(1:nn)=funftVI(i,1:nn);
    ppspline=spline(betaspline,funspline);
    funbVI1(1:nn)=ppval(ppspline,betacr(i,1:nn));
    vb_fII(1:nn,3*nn-3+i)=funbVI1(1:nn);
end
 % FILTRAGE
for iline1=1:nn,
   vb_fII(iline1,1:nb)=vb_fII(iline1,1:nb)*ftr;
end
% TRANSFERT FACE II+IV
for iline1=1:nn,
   funftII(iline1,1:nn)   = vb_fII(iline1,1:nn);
   funftIV(iline1,nn:-1:1)    = vb_fII(iline1,2*nn-1:3*nn-2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESEAU 5: ASSEMBLAGE DES DONNEES SUR LE RESEAU V-ALPHA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
va_fV=zeros(4*(nn-1),nn); % VALEURS =DONNEES RESEAU FACE V SELON ALPHA
funaII1=zeros(nn,1);
funaIV1=zeros(nn,1);
for jline1=1:nn, %% Face V : TRANSFERT OF DATA OF FACE V
    va_fV(1:nn-1,jline1)=funftV(1:nn-1,jline1);
end
% FACE II: SPLINE INTERPOLATION  A L'AIDE DES ANGLES ALFA
% INVERSION OF THE LOOP ON THE ISO-ETA LINES OF NETWORK VA 
% AND THE ETA LINES OF FACE II
for j=1:nn-1, % LOOP ON THE ETA-LINES OF FACE II. j= order of encountering !!!
    alfaspline(1:nn)=alfa(1:nn,nn+1-j);
    funspline(1:nn)=funftII(1:nn,nn+1-j);
    ppspline=spline(alfaspline,funspline);
    funaII1(1:nn)=ppval(ppspline,alfa1(1:nn,j)); % VERIFIER SI C'EST BIEN "J" SUR LE PLAN LOGIQUE: je ne suis pas sur!!!!!
    va_fV(nn-1+j,1:nn)=funaII1(1:nn);
end 
% FACE VI: TRANSFERT OF DATA
for jline1=1:nn,
    jbar= jline1-mm;
    etabar= jbar*deta;
    va_fV(2*nn-1:3*nn-3,jline1)=funftVI(nn:-1:2,jline1); % symetrie a bien regarder!
end
% FACE IV: SPLINE INTERPOLATION  A L'AIDE DES ANGLES ALFA
% INVERSION OF THE LOOP ON THE ISO-ETA LINES OF NETWORK VA 
% AND THE ETA LINES OF FACE IV
for j=1:nn-1, % LOOP ON THE ETA-LINES OF FACE IV. j= order of encountering !!!
    alfaspline(1:nn)=alfa(1:nn,j);
    funspline(1:nn)=funftIV(1:nn,j);
    ppspline=spline(alfaspline,funspline);
    funaIV1(1:nn)=ppval(ppspline,alfa1(1:nn,j)); % VERIFIER SI C'EST BIEN "J" SUR LE PLAN LOGIQUE: je ne suis pas sur!!!!!
    va_fV(3*nn-3+j,1:nn)=funaIV1(1:nn);
end 
% FILTRAGE
for jline1=1:nn,
   va_fV(1:na,jline1)=ftr*va_fV(1:na,jline1);
end
% TRANSFERT FACE V+VI
for jline1=1:nn,
   funftV(1:nn,jline1)   = va_fV(1:nn,jline1);
   funftVI(nn:-1:1,jline1)    = va_fV(2*nn-1:3*nn-2,jline1);
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
for iline1=1:nn,
   vb_fV(iline1,1:nb)=vb_fV(iline1,1:nb)*ftr;
end
% TRANSFERT FACE V+VI
for iline1=1:nn,
   funftV(iline1,1:nn)  = vb_fV(iline1,1:nn);
   funftVI(nn-iline1+1,1:nn)    = vb_fV(iline1,2*nn-1:3*nn-2);
end













