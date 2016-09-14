clc; clear all; close all;
% rot(grad(u))=0?
%% test rotationnel

clc; clear all; close all
global n nn
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr test
global teta0 teta1

n=40;
test = 0;
opt_ftr=0;
teta0=-3*pi/16;
teta1=3*pi/16;
mod72

%% *** initialisation des données
t=0;
[ ht_fI,    vt_fI] = sol_exacte(x_fI,   y_fI,   z_fI,   t);
[ ht_fII,   vt_fII] = sol_exacte(x_fII,  y_fII,  z_fII,  t);
[ ht_fIII,  vt_fIII] = sol_exacte(x_fIII, y_fIII, z_fIII, t);
[ ht_fIV,   vt_fIV] = sol_exacte(x_fIV,  y_fIV,  z_fIV,  t);
[ ht_fV,    vt_fV] = sol_exacte(x_fV,   y_fV,   z_fV,   t);
[ ht_fVI,   vt_fVI] = sol_exacte(x_fVI,  y_fVI,  z_fVI,  t);
    
% GRADIENT
[grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]= gr72( ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI  , n , nn);

[rot_fI,rot_fII,rot_fIII,rot_fIV,rot_fV,rot_fVI]=...
    rot75(grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI,n,nn);

err_I=abs(rot_fI);
err_II=abs(rot_fII);
err_III=abs(rot_fIII);
err_IV=abs(rot_fIV);
err_V=abs(rot_fV);
err_VI=abs(rot_fVI);

e_I=max(max(max(err_I)))
e_II=max(max(max(err_II)))
e_III=max(max(max(err_III)))
e_IV=max(max(max(err_IV)))
e_V=max(max(max(err_V)))
e_VI=max(max(max(err_VI)))

err=max([e_I e_II e_III e_IV e_V e_VI])


for p=1:3
    figure(p)
    plot_cs11(n,nn,err_I(:,:,p),err_II(:,:,p),err_III(:,:,p),err_IV(:,:,p),err_V(:,:,p),err_VI(:,:,p))
    title('error')
    colorbar
end
    
