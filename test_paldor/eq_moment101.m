function[eq1_fI, eq1_fII, eq1_fIII, eq1_fIV, eq1_fV, eq1_fVI]=eq_moment101(ht_fI,...
    ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI)
% momentum equation (right hand side) in shallow water system.
global nn n
global gp omega
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global gr_I gr_II gr_III gr_IV gr_V gr_VI



[grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
    gr101(ht_fI,ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI,n,nn);

[lambda_fI, teta_fI, ~]=cart2sph(x_fI, y_fI, z_fI);
[lambda_fII, teta_fII, ~]=cart2sph(x_fII, y_fII, z_fII);
[lambda_fIII, teta_fIII, ~]=cart2sph(x_fIII, y_fIII, z_fIII);
[lambda_fIV, teta_fIV, ~]=cart2sph(x_fIV, y_fIV, z_fIV);
[lambda_fV, teta_fV, ~]=cart2sph(x_fV, y_fV, z_fV);
[lambda_fVI, teta_fVI, ~]=cart2sph(x_fVI, y_fVI, z_fVI);

f_I=2.*omega.*sin(teta_fI);
f_II=2.*omega.*sin(teta_fII);
f_III=2.*omega.*sin(teta_fIII);
f_IV=2.*omega.*sin(teta_fIV);
f_V=2.*omega.*sin(teta_fV);
f_VI=2.*omega.*sin(teta_fVI);

% FACE I
vect_fI(1:nn,1:nn,1)=gr_I(1:nn,1:nn,2).*vt_fI(1:nn,1:nn,3)-gr_I(1:nn,1:nn,3).*vt_fI(1:nn,1:nn,2);
vect_fI(1:nn,1:nn,2)=gr_I(1:nn,1:nn,3).*vt_fI(1:nn,1:nn,1)-gr_I(1:nn,1:nn,1).*vt_fI(1:nn,1:nn,3);
vect_fI(1:nn,1:nn,3)=gr_I(1:nn,1:nn,1).*vt_fI(1:nn,1:nn,2)-gr_I(1:nn,1:nn,2).*vt_fI(1:nn,1:nn,1);
% FACE II
vect_fII(1:nn,1:nn,1)=gr_II(1:nn,1:nn,2).*vt_fII(1:nn,1:nn,3)-gr_II(1:nn,1:nn,3).*vt_fII(1:nn,1:nn,2);
vect_fII(1:nn,1:nn,2)=gr_II(1:nn,1:nn,3).*vt_fII(1:nn,1:nn,1)-gr_II(1:nn,1:nn,1).*vt_fII(1:nn,1:nn,3);
vect_fII(1:nn,1:nn,3)=gr_II(1:nn,1:nn,1).*vt_fII(1:nn,1:nn,2)-gr_II(1:nn,1:nn,2).*vt_fII(1:nn,1:nn,1);
% FACE III
vect_fIII(1:nn,1:nn,1)=gr_III(1:nn,1:nn,2).*vt_fIII(1:nn,1:nn,3)-gr_III(1:nn,1:nn,3).*vt_fIII(1:nn,1:nn,2);
vect_fIII(1:nn,1:nn,2)=gr_III(1:nn,1:nn,3).*vt_fIII(1:nn,1:nn,1)-gr_III(1:nn,1:nn,1).*vt_fIII(1:nn,1:nn,3);
vect_fIII(1:nn,1:nn,3)=gr_III(1:nn,1:nn,1).*vt_fIII(1:nn,1:nn,2)-gr_III(1:nn,1:nn,2).*vt_fIII(1:nn,1:nn,1);
% FACE IV
vect_fIV(1:nn,1:nn,1)=gr_IV(1:nn,1:nn,2).*vt_fIV(1:nn,1:nn,3)-gr_IV(1:nn,1:nn,3).*vt_fIV(1:nn,1:nn,2);
vect_fIV(1:nn,1:nn,2)=gr_IV(1:nn,1:nn,3).*vt_fIV(1:nn,1:nn,1)-gr_IV(1:nn,1:nn,1).*vt_fIV(1:nn,1:nn,3);
vect_fIV(1:nn,1:nn,3)=gr_IV(1:nn,1:nn,1).*vt_fIV(1:nn,1:nn,2)-gr_IV(1:nn,1:nn,2).*vt_fIV(1:nn,1:nn,1);
% FACE V
vect_fV(1:nn,1:nn,1)=gr_V(1:nn,1:nn,2).*vt_fV(1:nn,1:nn,3)-gr_V(1:nn,1:nn,3).*vt_fV(1:nn,1:nn,2);
vect_fV(1:nn,1:nn,2)=gr_V(1:nn,1:nn,3).*vt_fV(1:nn,1:nn,1)-gr_V(1:nn,1:nn,1).*vt_fV(1:nn,1:nn,3);
vect_fV(1:nn,1:nn,3)=gr_V(1:nn,1:nn,1).*vt_fV(1:nn,1:nn,2)-gr_V(1:nn,1:nn,2).*vt_fV(1:nn,1:nn,1);
% FACE VI
vect_fVI(1:nn,1:nn,1)=gr_VI(1:nn,1:nn,2).*vt_fVI(1:nn,1:nn,3)-gr_VI(1:nn,1:nn,3).*vt_fVI(1:nn,1:nn,2);
vect_fVI(1:nn,1:nn,2)=gr_VI(1:nn,1:nn,3).*vt_fVI(1:nn,1:nn,1)-gr_VI(1:nn,1:nn,1).*vt_fVI(1:nn,1:nn,3);
vect_fVI(1:nn,1:nn,3)=gr_VI(1:nn,1:nn,1).*vt_fVI(1:nn,1:nn,2)-gr_VI(1:nn,1:nn,2).*vt_fVI(1:nn,1:nn,1);


%% ASSEMBLAGE
eq1_fI(1:nn,1:nn,1)=-gp*grad_I(1:nn,1:nn,1)-f_I.*vect_fI(1:nn,1:nn,1);
eq1_fI(1:nn,1:nn,2)=-gp*grad_I(1:nn,1:nn,2)-f_I.*vect_fI(1:nn,1:nn,2);
eq1_fI(1:nn,1:nn,3)=-gp*grad_I(1:nn,1:nn,3)-f_I.*vect_fI(1:nn,1:nn,3);

eq1_fII(1:nn,1:nn,1)=-gp*grad_II(1:nn,1:nn,1)-f_II.*vect_fII(1:nn,1:nn,1);
eq1_fII(1:nn,1:nn,2)=-gp*grad_II(1:nn,1:nn,2)-f_II.*vect_fII(1:nn,1:nn,2);
eq1_fII(1:nn,1:nn,3)=-gp*grad_II(1:nn,1:nn,3)-f_II.*vect_fII(1:nn,1:nn,3);

eq1_fIII(1:nn,1:nn,1)=-gp*grad_III(1:nn,1:nn,1)-f_III.*vect_fIII(1:nn,1:nn,1);
eq1_fIII(1:nn,1:nn,2)=-gp*grad_III(1:nn,1:nn,2)-f_III.*vect_fIII(1:nn,1:nn,2);
eq1_fIII(1:nn,1:nn,3)=-gp*grad_III(1:nn,1:nn,3)-f_III.*vect_fIII(1:nn,1:nn,3);

eq1_fIV(1:nn,1:nn,1)=-gp*grad_IV(1:nn,1:nn,1)-f_IV.*vect_fIV(1:nn,1:nn,1);
eq1_fIV(1:nn,1:nn,2)=-gp*grad_IV(1:nn,1:nn,2)-f_IV.*vect_fIV(1:nn,1:nn,2);
eq1_fIV(1:nn,1:nn,3)=-gp*grad_IV(1:nn,1:nn,3)-f_IV.*vect_fIV(1:nn,1:nn,3);

eq1_fV(1:nn,1:nn,1)=-gp*grad_V(1:nn,1:nn,1)-f_V.*vect_fV(1:nn,1:nn,1);
eq1_fV(1:nn,1:nn,2)=-gp*grad_V(1:nn,1:nn,2)-f_V.*vect_fV(1:nn,1:nn,2);
eq1_fV(1:nn,1:nn,3)=-gp*grad_V(1:nn,1:nn,3)-f_V.*vect_fV(1:nn,1:nn,3);

eq1_fVI(1:nn,1:nn,1)=-gp*grad_VI(1:nn,1:nn,1)-f_VI.*vect_fVI(1:nn,1:nn,1);
eq1_fVI(1:nn,1:nn,2)=-gp*grad_VI(1:nn,1:nn,2)-f_VI.*vect_fVI(1:nn,1:nn,2);
eq1_fVI(1:nn,1:nn,3)=-gp*grad_VI(1:nn,1:nn,3)-f_VI.*vect_fVI(1:nn,1:nn,3);

