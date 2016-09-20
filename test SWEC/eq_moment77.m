function[u_fI, u_fII, u_fIII, u_fIV, u_fV, u_fVI]=eq_moment77(ht_fI,...
    ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI, ddt)
% solve equation of moment in Shallow Water equation with an Euler RSS
% Scheme.

global nn n
global gp omega
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global gr_I gr_II gr_III gr_IV gr_V gr_VI


%% terme en gradient
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
    gr74(fun_I,fun_II,fun_III,fun_IV,fun_V,fun_VI,n,nn);


%% terme en pdt vectoriel

% VORTICITY
[vort_fI,vort_fII,vort_fIII,vort_fIV,vort_fV,vort_fVI]=...
    vort74(vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI,n,nn);

% PDT VECT
[n1,n2]=size(x_fI);
for i=1:n1
    for j=1:n2
        vect_fI(i,j,1:3)   = cross(gr_I(i,j,1:3)   ,vt_fI(i,j,1:3)  );
        vect_fII(i,j,1:3)  = cross(gr_II(i,j,1:3)  ,vt_fII(i,j,1:3) );
        vect_fIII(i,j,1:3) = cross(gr_III(i,j,1:3) ,vt_fIII(i,j,1:3));
        vect_fIV(i,j,1:3)  = cross(gr_IV(i,j,1:3)  ,vt_fIV(i,j,1:3) );
        vect_fV(i,j,1:3)   = cross(gr_V(i,j,1:3)   ,vt_fV(i,j,1:3)  );
        vect_fVI(i,j,1:3)  = cross(gr_VI(i,j,1:3)  ,vt_fVI(i,j,1:3) );
    end
end

%% ASSEMBLAGE
[~, teta_fI, ~]   = cart2sph(x_fI   , y_fI    , z_fI);
[~, teta_fII, ~]  = cart2sph(x_fII  , y_fII   , z_fII);
[~, teta_fIII, ~] = cart2sph(x_fIII , y_fIII  , z_fIII);
[~, teta_fIV, ~]  = cart2sph(x_fIV  , y_fIV   , z_fIV);
[~, teta_fV, ~]   = cart2sph(x_fV   , y_fV    , z_fV);
[~, teta_fVI, ~]  = cart2sph(x_fVI  , y_fVI   , z_fVI);

for kk=1:3
    eq1_fI(:,:,kk)   = -ddt*  (grad_I(:,:,kk)   + (vort_fI   + 2*omega.*sin(teta_fI))   .*vect_fI(:,:,kk));
    eq1_fII(:,:,kk)  = -ddt*  (grad_II(:,:,kk)  + (vort_fII  + 2*omega.*sin(teta_fII))  .*vect_fII(:,:,kk));
    eq1_fIII(:,:,kk) = -ddt*  (grad_III(:,:,kk) + (vort_fIII + 2*omega.*sin(teta_fIII)) .*vect_fIII(:,:,kk));
    eq1_fIV(:,:,kk)  = -ddt*  (grad_IV(:,:,kk)  + (vort_fIV  + 2*omega.*sin(teta_fIV))  .*vect_fIV(:,:,kk));
    eq1_fV(:,:,kk)   = -ddt*  (grad_V(:,:,kk)   + (vort_fV   + 2*omega.*sin(teta_fV))   .*vect_fV(:,:,kk));
    eq1_fVI(:,:,kk)  = -ddt*  (grad_VI(:,:,kk)  + (vort_fVI  + 2*omega.*sin(teta_fVI))  .*vect_fVI(:,:,kk));
end

tau=0;
[zz_fI]   = solver_impli75(eq1_fI   ,vort_fI   ,ddt,tau,x_fI   ,y_fI   ,z_fI  );
[zz_fII]  = solver_impli75(eq1_fII  ,vort_fII  ,ddt,tau,x_fII  ,y_fII  ,z_fII );
[zz_fIII] = solver_impli75(eq1_fIII ,vort_fIII ,ddt,tau,x_fIII ,y_fIII ,z_fIII);
[zz_fIV]  = solver_impli75(eq1_fIV  ,vort_fIV  ,ddt,tau,x_fIV  ,y_fIV  ,z_fIV );
[zz_fV]   = solver_impli75(eq1_fV   ,vort_fV   ,ddt,tau,x_fV   ,y_fV   ,z_fV  );
[zz_fVI]  = solver_impli75(eq1_fVI  ,vort_fVI  ,ddt,tau,x_fVI  ,y_fVI  ,z_fVI );

u_fI   = vt_fI  +zz_fI;
u_fII  = vt_fII +zz_fII;
u_fIII = vt_fIII+zz_fIII;
u_fIV  = vt_fIV +zz_fIV;
u_fV   = vt_fV  +zz_fV;
u_fVI  = vt_fVI +zz_fVI;

