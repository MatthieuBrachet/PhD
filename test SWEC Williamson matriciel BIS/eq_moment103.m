function[eq1_fI, eq1_fII, eq1_fIII, eq1_fIV, eq1_fV, eq1_fVI]=eq_moment103(ht_fI,...
    ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI)
% momentum equation (right hand side) in shallow water system.
global nn n
global gp omega alpha
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global gr_I gr_II gr_III gr_IV gr_V gr_VI
%% gradient term
% face I
norm2_I=vt_fI(:,:,1).^2+vt_fI(:,:,2).^2+vt_fI(:,:,3).^2;
fun_I=0.5*norm2_I+gp*ht_fI;

% face II
norm2_II=vt_fII(:,:,1).^2+vt_fII(:,:,2).^2+vt_fII(:,:,3).^2;
fun_II=0.5*norm2_II+gp*ht_fII;

% face III
norm2_III=vt_fIII(:,:,1).^2+vt_fIII(:,:,2).^2+vt_fIII(:,:,3).^2;
fun_III=0.5*norm2_III+gp*ht_fIII;

% face IV
norm2_IV=vt_fIV(:,:,1).^2+vt_fIV(:,:,2).^2+vt_fIV(:,:,3).^2;
fun_IV=0.5*norm2_IV+gp*ht_fIV;

% face V
norm2_V=vt_fV(:,:,1).^2+vt_fV(:,:,2).^2+vt_fV(:,:,3).^2;
fun_V=0.5*norm2_V+gp*ht_fV;

% face VI
norm2_VI=vt_fVI(:,:,1).^2+vt_fVI(:,:,2).^2+vt_fVI(:,:,3).^2;
fun_VI=0.5*norm2_VI+gp*ht_fVI;

[grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
    gr103(fun_I,fun_II,fun_III,fun_IV,fun_V,fun_VI,n,nn);

%% vectorial product
% VORTICITY
[vort_fI,vort_fII,vort_fIII,vort_fIV,vort_fV,vort_fVI]=...
    vort103(vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI,n,nn);

% PDT VECT
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
[lambda_fI, teta_fI, ~]=cart2sph(x_fI, y_fI, z_fI);
[lambda_fII, teta_fII, ~]=cart2sph(x_fII, y_fII, z_fII);
[lambda_fIII, teta_fIII, ~]=cart2sph(x_fIII, y_fIII, z_fIII);
[lambda_fIV, teta_fIV, ~]=cart2sph(x_fIV, y_fIV, z_fIV);
[lambda_fV, teta_fV, ~]=cart2sph(x_fV, y_fV, z_fV);
[lambda_fVI, teta_fVI, ~]=cart2sph(x_fVI, y_fVI, z_fVI);

f_I=2.*omega.*(-cos(lambda_fI).*cos(teta_fI).*sin(alpha)+sin(teta_fI).*cos(alpha));
f_II=2.*omega.*(-cos(lambda_fII).*cos(teta_fII).*sin(alpha)+sin(teta_fII).*cos(alpha));
f_III=2.*omega.*(-cos(lambda_fIII).*cos(teta_fIII).*sin(alpha)+sin(teta_fIII).*cos(alpha));
f_IV=2.*omega.*(-cos(lambda_fIV).*cos(teta_fIV).*sin(alpha)+sin(teta_fIV).*cos(alpha));
f_V=2.*omega.*(-cos(lambda_fV).*cos(teta_fV).*sin(alpha)+sin(teta_fV).*cos(alpha));
f_VI=2.*omega.*(-cos(lambda_fVI).*cos(teta_fVI).*sin(alpha)+sin(teta_fVI).*cos(alpha));


eq1_fI(1:nn,1:nn,1)=-grad_I(1:nn,1:nn,1)-(vort_fI+f_I).*vect_fI(1:nn,1:nn,1);
eq1_fI(1:nn,1:nn,2)=-grad_I(1:nn,1:nn,2)-(vort_fI+f_I).*vect_fI(1:nn,1:nn,2);
eq1_fI(1:nn,1:nn,3)=-grad_I(1:nn,1:nn,3)-(vort_fI+f_I).*vect_fI(1:nn,1:nn,3);

eq1_fII(1:nn,1:nn,1)=-grad_II(1:nn,1:nn,1)-(vort_fII+f_II).*vect_fII(1:nn,1:nn,1);
eq1_fII(1:nn,1:nn,2)=-grad_II(1:nn,1:nn,2)-(vort_fII+f_II).*vect_fII(1:nn,1:nn,2);
eq1_fII(1:nn,1:nn,3)=-grad_II(1:nn,1:nn,3)-(vort_fII+f_II).*vect_fII(1:nn,1:nn,3);

eq1_fIII(1:nn,1:nn,1)=-grad_III(1:nn,1:nn,1)-(vort_fIII+f_III).*vect_fIII(1:nn,1:nn,1);
eq1_fIII(1:nn,1:nn,2)=-grad_III(1:nn,1:nn,2)-(vort_fIII+f_III).*vect_fIII(1:nn,1:nn,2);
eq1_fIII(1:nn,1:nn,3)=-grad_III(1:nn,1:nn,3)-(vort_fIII+f_III).*vect_fIII(1:nn,1:nn,3);

eq1_fIV(1:nn,1:nn,1)=-grad_IV(1:nn,1:nn,1)-(vort_fIV+f_IV).*vect_fIV(1:nn,1:nn,1);
eq1_fIV(1:nn,1:nn,2)=-grad_IV(1:nn,1:nn,2)-(vort_fIV+f_IV).*vect_fIV(1:nn,1:nn,2);
eq1_fIV(1:nn,1:nn,3)=-grad_IV(1:nn,1:nn,3)-(vort_fIV+f_IV).*vect_fIV(1:nn,1:nn,3);

eq1_fV(1:nn,1:nn,1)=-grad_V(1:nn,1:nn,1)-(vort_fV+f_V).*vect_fV(1:nn,1:nn,1);
eq1_fV(1:nn,1:nn,2)=-grad_V(1:nn,1:nn,2)-(vort_fV+f_V).*vect_fV(1:nn,1:nn,2);
eq1_fV(1:nn,1:nn,3)=-grad_V(1:nn,1:nn,3)-(vort_fV+f_V).*vect_fV(1:nn,1:nn,3);

eq1_fVI(1:nn,1:nn,1)=-grad_VI(1:nn,1:nn,1)-(vort_fVI+f_VI).*vect_fVI(1:nn,1:nn,1);
eq1_fVI(1:nn,1:nn,2)=-grad_VI(1:nn,1:nn,2)-(vort_fVI+f_VI).*vect_fVI(1:nn,1:nn,2);
eq1_fVI(1:nn,1:nn,3)=-grad_VI(1:nn,1:nn,3)-(vort_fVI+f_VI).*vect_fVI(1:nn,1:nn,3);

