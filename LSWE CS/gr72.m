function [grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
    gr72(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn)
global na nb;
global pts_beta ptscr_beta;
global pts_alfa ptscr_alfa;
global p k;
global gxi_I gxi_II gxi_III gxi_IV gxi_V gxi_VI;
global geta_I geta_II geta_III geta_IV geta_V geta_VI;

%% *** Derivées alpha - beta **********************************************

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

vad_fI=zeros(na,nn);
for jline1=1:nn, 
    funa7=va_fI(:,jline1);
    funad7=p\(k*funa7); 
    funad8=funad7;
    vad_fI(:,jline1)=funad8; 
end

%% RESEAU 2: ASSEMBLAGE DES DONNEES SUR LE RESEAU I-BETA
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

vbd_fI=zeros(nn,nb);
for iline1=1:nn,
    funb1=vb_fI(iline1,:);
    funbd1=p\(k*funb1');
    funbd2=funbd1;
    vbd_fI(iline1,:)=funbd2;
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

vad_fII=zeros(na,nn);
for jline1=1:nn,
    funa9=va_fII(:,jline1);
    funad9=p\(k*funa9);
    funad10=funad9;
    vad_fII(:,jline1)=funad10;
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

vbd_fII=zeros(nn,nb);
for iline1=1:nn,
    funb2=vb_fII(iline1,:);
    funbd2=p\(k*funb2');
    funbd3=funbd2;
    vbd_fII(iline1,:)=funbd3; 
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

vad_fV=zeros(na,nn);
for jline1=1:nn,
    funa11=va_fV(:,jline1);
    funad11=p\(k*funa11);
    funad12=funad11;
    vad_fV(:,jline1)=funad12; 
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

vbd_fV=zeros(nn,nb);
for iline1=1:nn,
    funb5=vb_fV(iline1,:);
    funbd5=p\(k*funb5');
    funbd6=funbd5;
    vbd_fV(iline1,:)=funbd6; 
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

%% *** Gradient ***********************************************************

grad_I=zeros(nn,nn,3);
for i=1:nn,
    for j=1:nn,
        grad_I(i,j,1:3)=dg_alfa(i,j,1)*gxi_I(i,j,1:3) + dg_beta(i,j,1)*geta_I(i,j,1:3);
    end
end

grad_II=zeros(nn,nn,3);
for i=1:nn,
    for j=1:nn,
        grad_II(i,j,1:3)=dg_alfa(i,j,2)*gxi_II(i,j,1:3) + dg_beta(i,j,2)*geta_II(i,j,1:3);
    end
end

grad_III=zeros(nn,nn,3);
for i=1:nn,
    for j=1:nn,
        grad_III(i,j,1:3)=dg_alfa(i,j,3)*gxi_III(i,j,1:3) - dg_beta(i,j,3)*geta_III(i,j,1:3);
    end
end

grad_IV=zeros(nn,nn,3);
for i=1:nn,
    for j=1:nn,
        grad_IV(i,j,1:3)=dg_alfa(i,j,4)*gxi_IV(i,j,1:3) - dg_beta(i,j,4)*geta_IV(i,j,1:3);
    end
end

grad_V=zeros(nn,nn,3);
for i=1:nn,
    for j=1:nn,
        grad_V(i,j,1:3)=dg_alfa(i,j,5)*gxi_V(i,j,1:3) + dg_beta(i,j,5)*geta_V(i,j,1:3);
    end
end

grad_VI=zeros(nn,nn,3);

for i=1:nn,
    for j=1:nn,
        grad_VI(i,j,1:3)=dg_alfa(i,j,6)*gxi_VI(i,j,1:3) + dg_beta(i,j,6)*geta_VI(i,j,1:3);
    end
end

%% 1/2 SOMME ARRETES, 1/3 SOMME SOMMETS
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






