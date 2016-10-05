function[eq1_fI, eq1_fII, eq1_fIII, eq1_fIV, eq1_fV, eq1_fVI]=eq_moment75(ht_fI,...
    ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI)

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
        vect_fI(i,j,1:3)=cross(gr_I(i,j,1:3),vt_fI(i,j,1:3));
        vect_fII(i,j,1:3)=cross(gr_II(i,j,1:3),vt_fII(i,j,1:3));
        vect_fIII(i,j,1:3)=cross(gr_III(i,j,1:3),vt_fIII(i,j,1:3));
        vect_fIV(i,j,1:3)=cross(gr_IV(i,j,1:3),vt_fIV(i,j,1:3));
        vect_fV(i,j,1:3)=cross(gr_V(i,j,1:3),vt_fV(i,j,1:3));
        vect_fVI(i,j,1:3)=cross(gr_VI(i,j,1:3),vt_fVI(i,j,1:3));
    end
end

%% ASSEMBLAGE
[~, teta_fI, ~]=cart2sph(x_fI, y_fI, z_fI);
[~, teta_fII, ~]=cart2sph(x_fII, y_fII, z_fII);
[~, teta_fIII, ~]=cart2sph(x_fIII, y_fIII, z_fIII);
[~, teta_fIV, ~]=cart2sph(x_fIV, y_fIV, z_fIV);
[~, teta_fV, ~]=cart2sph(x_fV, y_fV, z_fV);
[~, teta_fVI, ~]=cart2sph(x_fVI, y_fVI, z_fVI);

[cor_fI] = coriolis_mat74(vt_fI, x_fI, y_fI, z_fI);
[cor_fII] = coriolis_mat74(vt_fII, x_fII, y_fII, z_fII);
[cor_fIII] = coriolis_mat74(vt_fIII, x_fIII, y_fIII, z_fIII);
[cor_fIV] = coriolis_mat74(vt_fIV, x_fIV, y_fIV, z_fIV);
[cor_fV] = coriolis_mat74(vt_fV, x_fV, y_fV, z_fV);
[cor_fVI] = coriolis_mat74(vt_fVI, x_fVI, y_fVI, z_fVI);




for kk=1:3
    eq1_fI(:,:,kk)=-grad_I(:,:,kk)-(vort_fI).*vect_fI(:,:,kk)-2*omega*sin(teta_fI).*cor_fI(:,:,kk);
    eq1_fII(:,:,kk)=-grad_II(:,:,kk)-(vort_fII).*vect_fII(:,:,kk)-2*omega*sin(teta_fII).*cor_fII(:,:,kk);
    eq1_fIII(:,:,kk)=-grad_III(:,:,kk)-(vort_fIII).*vect_fIII(:,:,kk)-2*omega*sin(teta_fIII).*cor_fIII(:,:,kk);
    eq1_fIV(:,:,kk)=-grad_IV(:,:,kk)-(vort_fIV).*vect_fIV(:,:,kk)-2*omega*sin(teta_fIV).*cor_fIV(:,:,kk);
    eq1_fV(:,:,kk)=-grad_V(:,:,kk)-(vort_fV).*vect_fV(:,:,kk)-2*omega*sin(teta_fV).*cor_fV(:,:,kk);
    eq1_fVI(:,:,kk)=-grad_VI(:,:,kk)-(vort_fVI).*vect_fVI(:,:,kk)-2*omega*sin(teta_fVI).*cor_fVI(:,:,kk);
end