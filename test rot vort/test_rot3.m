clc; clear all; close all;
% div(rot(u))=0?
%% test rotationnel

clc; clear all; close all
global n nn
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr test
global teta0 teta1

n=100;
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

[rot_fI,rot_fII,rot_fIII,rot_fIV,rot_fV,rot_fVI]=...
    rot75(vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV,vt_fVI ,n,nn);

[div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
    div72(rot_fI,rot_fII,rot_fIII,rot_fIV,rot_fV,rot_fVI,n,nn);


err_I=abs(div_fI);
err_II=abs(div_fII);
err_III=abs(div_fIII);
err_IV=abs(div_fIV);
err_V=abs(div_fV);
err_VI=abs(div_fVI);

e_I=max(max(max(err_I)))
e_II=max(max(max(err_II)))
e_III=max(max(max(err_III)))
e_IV=max(max(max(err_IV)))
e_V=max(max(max(err_V)))
e_VI=max(max(max(err_VI)))

err=max([e_I e_II e_III e_IV e_V e_VI])


for p=1:3
    figure(p)
    plot_cs11(n,nn,div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI)
    title('error')
    colorbar
end
    