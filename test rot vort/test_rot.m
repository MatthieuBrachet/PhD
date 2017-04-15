%% test rotational
clc; clear all; close all
global n nn
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr
global teta0 teta1



teta0=-3*pi/16;
teta1=3*pi/16;
opt_ftr=0;
n=127;
mod72;

[vt_fI,rotvt_fI] = sol_exacte_rot(x_fI,y_fI,z_fI);
[vt_fII,rotvt_fII] = sol_exacte_rot(x_fII,y_fII,z_fII);
[vt_fIII,rotvt_fIII] = sol_exacte_rot(x_fIII,y_fIII,z_fIII);
[vt_fIV,rotvt_fIV] = sol_exacte_rot(x_fIV,y_fIV,z_fIV);
[vt_fV,rotvt_fV] = sol_exacte_rot(x_fV,y_fV,z_fV);
[vt_fVI,rotvt_fVI] = sol_exacte_rot(x_fVI,y_fVI,z_fVI);

[rot_fI,rot_fII,rot_fIII,rot_fIV,rot_fV,rot_fVI]=...
    rot74(vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI,n,nn);

err_I=abs(rotvt_fI-rot_fI);
err_II=abs(rotvt_fII-rot_fII);
err_III=abs(rotvt_fIII-rot_fIII);
err_IV=abs(rotvt_fIV-rot_fIV);
err_V=abs(rotvt_fV-rot_fV);
err_VI=abs(rotvt_fVI-rot_fVI);

e_I=max(max(max(err_I)));
e_II=max(max(max(err_II)));
e_III=max(max(max(err_III)));
e_IV=max(max(max(err_IV)));
e_V=max(max(max(err_V)));
e_VI=max(max(max(err_VI)));

MMM=max(max(max(abs([rotvt_fI rotvt_fII rotvt_fIII rotvt_fIV rotvt_fV rotvt_fVI]))));

err=max([e_I e_II e_III e_IV e_V e_VI])./MMM

for pp=1:3
    figure(pp)
    subplot(1,3,1)
    plot_cs11(n,nn,rotvt_fI(:,:,pp),rotvt_fII(:,:,pp),rotvt_fIII(:,:,pp),rotvt_fIV(:,:,pp),rotvt_fV(:,:,pp),rotvt_fVI(:,:,pp))
    title('exact sol.')
    colorbar

    subplot(1,3,2)
    plot_cs11(n,nn,rot_fI(:,:,pp),rot_fII(:,:,pp),rot_fIII(:,:,pp),rot_fIV(:,:,pp),rot_fV(:,:,pp),rot_fVI(:,:,pp))
    title('approx sol.')
    colorbar
    
    subplot(1,3,3)
    plot_cs11(n,nn,err_I(:,:,pp),err_II(:,:,pp),err_III(:,:,pp),err_IV(:,:,pp),err_V(:,:,pp),err_VI(:,:,pp))
    title('ERROR')
    colorbar
end


