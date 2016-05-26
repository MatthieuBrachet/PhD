function [funftI,funftII,funftIII,funftIV,funftV,funftVI]=...
    ftr72(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn)
global na nb;
global pts_beta ptscr_beta;
global pts_alfa ptscr_alfa;
global ftr;

betaspline=zeros(nn,1);
alfaspline=zeros(nn,1);
funspline=zeros(nn,1);
funftI(1:nn,1:nn)=funfI(1:nn,1:nn);
funftII(1:nn,1:nn)=funfII(1:nn,1:nn);
funftIII(1:nn,1:nn)=funfIII(1:nn,1:nn);
funftIV(1:nn,1:nn)=funfIV(1:nn,1:nn);
funftV(1:nn,1:nn)=funfV(1:nn,1:nn);
funftVI(1:nn,1:nn)=funfVI(1:nn,1:nn);

%% RESEAU 1 : ASSEMBLAGE DES DONNEES SUR LE RESEAU I-ALPHA

% Face I : TRANSFERT OF DATA OF FACE I 
va_fI=zeros(4*(nn-1),nn); 
for jline1=1:nn, % boucle sur les lignes iso-eta du reseaux Ia
    va_fI(1:nn-1,jline1)=funfI(1:nn-1,jline1);
end

% Face II
fun=zeros(1,nn*nn);
FUN=zeros(nn-1,nn);
compteur=1;
for j=1:nn
    for i=1:nn
        if rem(j,2) == 1
            fun(compteur)=funfII(j,i);
        elseif rem(j,2) == 0
            fun(compteur)=funfII(j,nn-i+1);
        end
            compteur=compteur+1;
    end
end
ppspline=spline(pts_beta,fun);
funspl=ppval(ppspline,ptscr_beta(1:nn*nn));

for i=1:nn-1
    if i <= nn/2
        if rem(i,2)== 1
            FUN(i,1:nn)=funspl((i-1)*nn+1:(i-1)*nn+nn);
        else
            FUN(i,1:nn)=funspl(i*nn-[1:nn]+1);
        end
    else
        if rem(i,2)== 1
            FUN(i,1:nn)=funspl((i-1)*nn+[nn:-1:1]);
        else
            FUN(i,1:nn)=funspl(i*nn-[nn:-1:1]+1);
        end
    end
end
va_fI(nn-1+[1:nn-1],1:nn)=FUN(1:nn-1,1:nn);

% Face III : transfert of data of face III
for jline1=1:nn,
    va_fI(2*nn-1:3*nn-3,jline1)=funfIII(1:nn-1,nn-jline1+1);
end

% Face IV
fun=zeros(1,nn*nn);
FUN=zeros(nn-1,nn);
compteur=1;
for j=1:nn
    for i=1:nn
        if rem(j,2) == 1
            fun(compteur)=funfIV(j,i);
        elseif rem(j,2) == 0
            fun(compteur)=funfIV(j,nn-i+1);
        end
            compteur=compteur+1;
    end
end
ppspline=spline(pts_beta,fun);
funspl=ppval(ppspline,ptscr_beta);

for i=1:nn-1
    if i > nn/2
        if rem(i,2)== 1
            FUN(i,1:nn)=funspl((i-1)*nn+[1:nn]);
        else
            FUN(i,1:nn)=funspl(i*nn-[1:nn]+1);
        end
    else
        if rem(i,2)== 1
            FUN(i,1:nn)=funspl((i-1)*nn+[nn:-1:1]);
        else
            FUN(i,1:nn)=funspl(i*nn-[nn:-1:1]+1);
        end
    end
end
va_fI(3*nn-3+[1:nn-1],1:nn)=FUN;

% FILTRAGE
for jline1=1:nn,
   va_fI(1:na,jline1)=ftr*va_fI(1:na,jline1);
end

% TRANSFERT FACE I+III
for jline1=1:nn,
   funftI(1:nn,jline1)   = va_fI(1:nn,jline1);
   funftIII(1:nn,nn-jline1+1)    = va_fI(2*nn-1:3*nn-2,jline1);
end

%% RESEAU 2: ASSEMBLAGE DES DONNEES SUR LE RESEAU I-BETA
vb_fI=zeros(nn,4*(nn-1));
funaV1=zeros(nn,1);
funaVI1=zeros(nn,1);

vb_fI=zeros(nn,4*(nn-1));
% Face I : TRANSFERT OF DATA OF FACE I 
for iline1=1:nn, % boucle sur les lignes iso-xi du reseau I-beta
    vb_fI(iline1,1:nn-1)=funfI(iline1,1:nn-1);
end

% Face V
fun=zeros(1,nn*nn);
FUN=zeros(nn,nn-1);
compteur=1;
for j=1:nn
    for i=1:nn
        if rem(j,2) == 1
            fun(compteur)=funfV(i,j);
        elseif rem(j,2) == 0
            fun(compteur)=funfV(nn-i+1,j);
        end
            compteur=compteur+1;
    end
end
ppspline=spline(pts_alfa,fun);
funspl=ppval(ppspline,ptscr_alfa(1:nn*nn));

for i=1:nn-1
    if i <= nn/2
        if rem(i,2)== 1
            FUN(1:nn,i)=funspl((i-1)*nn+[1:nn]);
        else
            FUN(1:nn,i)=funspl(i*nn-[1:nn]+1);
        end
    else
        if rem(i,2)== 1
            FUN(1:nn,i)=funspl((i-1)*nn+[nn:-1:1]);
        else
            FUN(1:nn,i)=funspl(i*nn-[nn:-1:1]+1);
        end
    end
end
vb_fI(1:nn,nn-1+[1:nn-1])=FUN(1:nn,1:nn-1);

% FACE III: TRANSFERT OF DATA
 for iline1=1:nn,
     vb_fI(iline1,2*nn-1:3*nn-3)=funfIII(iline1,nn:-1:2); % symetrie sur adresses en i + inversion en j !
 end

% Face VI
fun=zeros(1,nn*nn);
FUN=zeros(nn,nn-1);
compteur=1;
for j=1:nn
    for i=1:nn
        if rem(j,2) == 1
            fun(compteur)=funfVI(i,j);
        elseif rem(j,2) == 0
            fun(compteur)=funfVI(nn-i+1,j);
        end
            compteur=compteur+1;
    end
end
ppspline=spline(pts_beta,fun);
funspl=ppval(ppspline,ptscr_beta);

for i=1:nn-1
    if i > nn/2
        if rem(i,2)== 1
            FUN(1:nn,i)=funspl((i-1)*nn+[1:nn]);
        else
            FUN(1:nn,i)=funspl(i*nn-[1:nn]+1);
        end
    else
        if rem(i,2)== 1
            FUN(1:nn,i)=funspl((i-1)*nn+[nn:-1:1]);
        else
            FUN(1:nn,i)=funspl(i*nn-[nn:-1:1]+1);
        end
    end
end
vb_fI(1:nn,3*nn-3+[1:nn-1])=FUN(1:nn,1:nn-1);
 
 % FILTRAGE
for iline1=1:nn,
   vb_fI(iline1,1:nb)=vb_fI(iline1,1:nb)*ftr;
end

% TRANSFERT FACE I+III
for iline1=1:nn,
   funftI(iline1,1:nn)   = vb_fI(iline1,1:nn);
   funftIII(iline1,nn:-1:1)    = vb_fI(iline1,2*nn-1:3*nn-2);
end


%% RESEAU 3: ASSEMBLAGE DES DONNEES SUR LE RESEAU II-ALPHA
% Face II : TRANSFERT OF DATA OF FACE I 
va_fII=zeros(4*(nn-1),nn); 
for jline1=1:nn, % boucle sur les lignes iso-eta du reseaux Ia
    va_fII(1:nn-1,jline1)=funfII(1:nn-1,jline1);
end

% Face III
fun=zeros(1,nn*nn);
FUN=zeros(nn-1,nn);
compteur=1;
for j=1:nn
    for i=1:nn
        if rem(j,2) == 1
            fun(compteur)=funfIII(j,i);
        elseif rem(j,2) == 0
            fun(compteur)=funfIII(j,nn-i+1);
        end
            compteur=compteur+1;
    end
end
ppspline=spline(pts_beta,fun);
funspl=ppval(ppspline,ptscr_beta(1:nn*nn));

for i=1:nn-1
    if i <= nn/2
        if rem(i,2)== 1
            FUN(i,1:nn)=funspl((i-1)*nn+[1:nn]);
        else
            FUN(i,1:nn)=funspl(i*nn-[1:nn]+1);
        end
    else
        if rem(i,2)== 1
            FUN(i,1:nn)=funspl((i-1)*nn+[nn:-1:1]);
        else
            FUN(i,1:nn)=funspl(i*nn-[nn:-1:1]+1);
        end
    end
end
va_fII(nn-1+[1:nn-1],1:nn)=FUN(1:nn-1,1:nn);

% Face IV 
for jline1=1:nn,
    va_fII(2*nn-1:3*nn-3,jline1)=funfIV(1:nn-1,nn-jline1+1);
end

% Face I
fun=zeros(1,nn*nn);
FUN=zeros(nn-1,nn);
compteur=1;
for j=1:nn
    for i=1:nn
        if rem(j,2) == 1
            fun(compteur)=funfI(j,i);
        elseif rem(j,2) == 0
            fun(compteur)=funfI(j,nn-i+1);
        end
            compteur=compteur+1;
    end
end
ppspline=spline(pts_beta,fun);
funspl=ppval(ppspline,ptscr_beta);

for i=1:nn-1
    if i > nn/2
        if rem(i,2)== 1
            FUN(i,1:nn)=funspl((i-1)*nn+[1:nn]);
        else
            FUN(i,1:nn)=funspl(i*nn-[1:nn]+1);
        end
    else
        if rem(i,2)== 1
            FUN(i,1:nn)=funspl((i-1)*nn+[nn:-1:1]);
        else
            FUN(i,1:nn)=funspl(i*nn-[nn:-1:1]+1);
        end
    end
end
va_fII(3*nn-3+[1:nn-1],1:nn)=FUN(1:nn-1,1:nn);

% FILTRAGE
for jline1=1:nn,
   va_fII(1:na,jline1)=ftr*va_fII(1:na,jline1);
end

% TRANSFERT FACE II+IV
for jline1=1:nn,
   funftII(1:nn,jline1)   = va_fII(1:nn,jline1);
   funftIV(1:nn,nn-jline1+1)    = va_fII(2*nn-1:3*nn-2,jline1);
end

%% RESEAU 4: ASSEMBLAGE DES DONNEES SUR LE RESEAU II-BETA
vb_fII=zeros(nn,4*(nn-1));
% Face II 
for iline1=1:nn, 
    vb_fII(iline1,1:nn-1)=funfII(iline1,1:nn-1);
end

% Face V
fun=zeros(1,nn*nn);
FUN=zeros(nn,nn-1);
compteur=1;
for j=1:nn
    for i=1:nn
        if rem(j,2) == 1
            fun(compteur)=funfV(nn-j+1,i);
        elseif rem(j,2) == 0
            fun(compteur)=funfV(nn-j+1,nn-i+1);
        end
            compteur=compteur+1;
    end
end
ppspline=spline(pts_beta,fun);
funspl=ppval(ppspline,ptscr_beta);
for i=1:nn-1
    if i > nn/2
        if rem(i,2)== 1
            FUN(1:nn,i)=funspl((i-1)*nn+[1:nn]);
        else
            FUN(1:nn,i)=funspl(i*nn-[1:nn]+1);
        end
    else
        if rem(i,2)== 1
            FUN(1:nn,i)=funspl((i-1)*nn+[nn:-1:1]);
        else
            FUN(1:nn,i)=funspl(i*nn-[nn:-1:1]+1);
        end
    end
end
vb_fII(1:nn,nn-1+[1:nn-1])=FUN(nn:-1:1,1:nn-1);

% FACE III: TRANSFERT OF DATA OF FACE III
for iline1=1:nn,
    vb_fII(iline1,2*nn-1:3*nn-3)=funfIV(iline1,nn:-1:2); % symetrie sur adresses en i + inversion en j !
end

% Face VI
fun=zeros(1,nn*nn);
FUN=zeros(nn,nn-1);
compteur=1;
for j=1:nn
    for i=1:nn
        if rem(j,2) == 1
            fun(compteur)=funfVI(j,i);
        elseif rem(j,2) == 0
            fun(compteur)=funfVI(j,nn-i+1);
        end
            compteur=compteur+1;
    end
end
ppspline=spline(pts_beta,fun);
funspl=ppval(ppspline,ptscr_beta);

for i=1:nn-1
    if i > nn/2
        if rem(i,2)== 1
            FUN(1:nn,i)=funspl((i-1)*nn+[1:nn]);
        else
            FUN(1:nn,i)=funspl(i*nn-[1:nn]+1);
        end
    else
        if rem(i,2)== 1
            FUN(1:nn,i)=funspl((i-1)*nn+[nn:-1:1]);
        else
            FUN(1:nn,i)=funspl(i*nn-[nn:-1:1]+1);
        end
    end
end
vb_fII(1:nn,3*nn-3+[1:nn-1])=FUN(nn:-1:1,1:nn-1);

 % FILTRAGE
for iline1=1:nn,
   vb_fII(iline1,1:nb)=vb_fII(iline1,1:nb)*ftr;
end

% TRANSFERT FACE II+IV
for iline1=1:nn,
   funftII(iline1,1:nn)   = vb_fII(iline1,1:nn);
   funftIV(iline1,nn:-1:1)    = vb_fII(iline1,2*nn-1:3*nn-2);
end

%% RESEAU 5: ASSEMBLAGE DES DONNEES SUR LE RESEAU V-ALPHA
vad_fV=zeros(na,nn);
va_fV=zeros(4*(nn-1),nn); 
% Face V
for jline1=1:nn, %% Face V : TRANSFERT OF DATA OF FACE V
    va_fV(1:nn-1,jline1)=funfV(1:nn-1,jline1);
end

% Face II
fun=zeros(1,nn*nn);
FUN=zeros(nn-1,nn);
compteur=1;
for j=1:nn
    for i=1:nn
        if rem(j,2) == 0
            fun(compteur)=funfII(i,nn+1-j);
        elseif rem(j,2) == 1
            fun(compteur)=funfII(nn+1-i,nn+1-j);
        end
            compteur=compteur+1;
    end
end
ppspline=spline(pts_alfa,fun);
funspl=ppval(ppspline,ptscr_alfa(1:nn*nn));

for i=1:nn-1
    if i <= nn/2
        if rem(i,2)== 0
            FUN(i,1:nn)=funspl((i-1)*nn+[1:nn]);
        else
            FUN(i,1:nn)=funspl(i*nn-[1:nn]+1);
        end
    else
        if rem(i,2)== 0
            FUN(i,1:nn)=funspl((i-1)*nn+[nn:-1:1]);
        else
            FUN(i,1:nn)=funspl(i*nn-[nn:-1:1]+1);
        end
    end
end
va_fV(nn-1+[1:nn-1],1:nn)=FUN(1:nn-1,1:nn);

% FACE VI: TRANSFERT OF DATA
for jline1=1:nn,
    va_fV(2*nn-1:3*nn-3,jline1)=funfVI(nn:-1:2,jline1); % symetrie a bien regarder!
end

% Face IV
fun=zeros(1,nn*nn);
FUN=zeros(nn-1,nn);
compteur=1;
for j=1:nn
    for i=1:nn
        if rem(j,2) == 0
            fun(compteur)=funfIV(nn+1-i,j);
        elseif rem(j,2) == 1
            fun(compteur)=funfIV(i,j);
        end
            compteur=compteur+1;
    end
end
ppspline=spline(pts_alfa,fun);
funspl=ppval(ppspline,ptscr_alfa(1:nn*nn));

for i=1:nn-1
    if i <= nn/2
        if rem(i,2)== 1
            FUN(i,1:nn)=funspl((i-1)*nn+[1:nn]);
        else
            FUN(i,1:nn)=funspl(i*nn-[1:nn]+1);
        end
    else
        if rem(i,2)== 1
            FUN(i,1:nn)=funspl((i-1)*nn+[nn:-1:1]);
        else
            FUN(i,1:nn)=funspl(i*nn-[nn:-1:1]+1);
        end
    end
end
 va_fV(3*nn-3+[1:nn-1],1:nn)=FUN(1:nn-1,1:nn);


% FILTRAGE
for jline1=1:nn,
   va_fV(1:na,jline1)=ftr*va_fV(1:na,jline1);
end

% TRANSFERT FACE V+VI
for jline1=1:nn,
   funftV(1:nn,jline1)   = va_fV(1:nn,jline1);
   funftVI(nn:-1:1,jline1)    = va_fV(2*nn-1:3*nn-2,jline1);
end

%% RESEAU 6: ASSEMBLAGE DES DONNEES SUR LE RESEAU V-BETA
vb_fV=zeros(nn,4*(nn-1));

% Face V : TRANSFERT OF DATA selon beta
for iline1=1:nn, % boucle sur les lignes iso-xi du reseau V-beta
    vb_fV(iline1,1:nn-1)=funfV(iline1,1:nn-1);
end

% Face III
fun=zeros(1,nn*nn);
FUN=zeros(nn,nn-1);
compteur=1;
for j=1:nn
    for i=1:nn
        if rem(j,2) == 0
            fun(compteur)=funfIII(i,nn+1-j);
        elseif rem(j,2) == 1
            fun(compteur)=funfIII(nn+1-i,nn+1-j);
        end
            compteur=compteur+1;
    end
end
ppspline=spline(pts_alfa,fun);
funspl=ppval(ppspline,ptscr_alfa(1:nn*nn));

for i=1:nn-1
    if i <= nn/2
        if rem(i,2)== 0
            FUN(1:nn,i)=funspl((i-1)*nn+[1:nn]);
        else
            FUN(1:nn,i)=funspl(i*nn-[1:nn]+1);
        end
    else
        if rem(i,2)== 0
            FUN(1:nn,i)=funspl((i-1)*nn+[nn:-1:1]);
        else
            FUN(1:nn,i)=funspl(i*nn-[nn:-1:1]+1);
        end
    end
end
vb_fV(1:nn,nn:2*nn-2)=FUN(nn:-1:1,1:nn-1);

% Face V
for iline1=1:nn,
    vb_fV(iline1,2*nn-1:3*nn-3)=funfVI(nn-iline1+1,1:nn-1);
end

% Face I
fun=zeros(1,nn*nn);
FUN=zeros(nn,nn-1);
compteur=1;
for j=1:nn
    for i=1:nn
        if rem(j,2) == 0
            fun(compteur)=funfI(nn+1-i,j);
        elseif rem(j,2) == 1
            fun(compteur)=funfI(i,j);
        end
            compteur=compteur+1;
    end
end
ppspline=spline(pts_alfa,fun);
funspl=ppval(ppspline,ptscr_alfa(1:nn*nn));

for i=1:nn-1
    if i <= nn/2
        if rem(i,2)== 1
            FUN(1:nn,i)=funspl((i-1)*nn+[1:nn]);
        else
            FUN(1:nn,i)=funspl(i*nn-[1:nn]+1);
        end
    else
        if rem(i,2)== 1
            FUN(1:nn,i)=funspl((i-1)*nn+[nn:-1:1]);
        else
            FUN(1:nn,i)=funspl(i*nn-[nn:-1:1]+1);
        end
    end
end
 vb_fV(1:nn,3*nn-3+[1:nn-1])=FUN(nn:-1:1,1:nn-1);


 % FILTRAGE
for iline1=1:nn,
   vb_fV(iline1,1:nb)=vb_fV(iline1,1:nb)*ftr;
end

% TRANSFERT FACE V+VI
for iline1=1:nn,
   funftV(iline1,1:nn)  = vb_fV(iline1,1:nn);
   funftVI(nn-iline1+1,1:nn)    = vb_fV(iline1,2*nn-1:3*nn-2);
end