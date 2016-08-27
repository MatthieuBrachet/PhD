clc; clear all; close all; format short;

global n nn
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr test

global hp gp

test=3;
opt_ftr=10;
n=500;
mod72

% time=10;
% sigma=1/5000;
% 
% %% fonction
% [ ht_fI,    vt_fI] = sol_exacte(x_fI,   y_fI,   z_fI,   time);
% [ ht_fII,   vt_fII] = sol_exacte(x_fII,  y_fII,  z_fII,  time);
% [ ht_fIII,  vt_fIII] = sol_exacte(x_fIII, y_fIII, z_fIII, time);
% [ ht_fIV,   vt_fIV] = sol_exacte(x_fIV,  y_fIV,  z_fIV,  time);
% [ ht_fV,    vt_fV] = sol_exacte(x_fV,   y_fV,   z_fV,   time);
% [ ht_fVI,   vt_fVI] = sol_exacte(x_fVI,  y_fVI,  z_fVI,  time);
% 
% 
% %% equation de conservation
% [div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
%     div72(vt_fI,vt_fII,vt_fIII,vt_fIV,vt_fV,vt_fVI,n,nn);
% 
% [forh_I] = for_h(x_fI,y_fI,z_fI,time);
% [forh_II] = for_h(x_fII,y_fII,z_fII,time);
% [forh_III] = for_h(x_fIII,y_fIII,z_fIII,time);
% [forh_IV] = for_h(x_fIV,y_fIV,z_fIV,time);
% [forh_V] = for_h(x_fV,y_fV,z_fV,time);
% [forh_VI] = for_h(x_fVI,y_fVI,z_fVI,time);
% 
% err_fI=abs(hp*div_fI-forh_I);
% err_fII=abs(hp*div_fII-forh_II);
% err_fIII=abs(hp*div_fIII-forh_III);
% err_fIV=abs(hp*div_fIV-forh_IV);
% err_fV=abs(hp*div_fV-forh_V);
% err_fVI=abs(hp*div_fVI-forh_VI);
% 
% err_h=max(max([err_fI, err_fII, err_fIII, err_fIV, err_fV, err_fVI]))
% 
% figure(1)
% plot_cs11(n,nn,err_fI, err_fII,err_fIII,err_fIV,err_fV,err_fVI);
% 
% figure(2)
% subplot(121)
% plot_cs11(n,nn,div_fI, div_fII,div_fIII,div_fIV,div_fV,div_fVI);
% title('divergence calculee')
% 
% subplot(122)
% plot_cs11(n,nn,forh_I, forh_II,forh_III,forh_IV,forh_V,forh_VI);
% title('divergence exacte')
% 
% %% equation de moment
% [grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
%     gr72(ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,n,nn);
% 
% [cor_I,cor_II,cor_III,cor_IV,cor_V,cor_VI]=coriolis(vt_fI,vt_fII,vt_fIII,vt_fIV,vt_fV,vt_fVI);
% 
% [forv_I] = for_v(x_fI,y_fI,z_fI,time);
% [forv_II] = for_v(x_fII,y_fII,z_fII,time);
% [forv_III] = for_v(x_fIII,y_fIII,z_fIII,time);
% [forv_IV] = for_v(x_fIV,y_fIV,z_fIV,time);
% [forv_V] = for_v(x_fV,y_fV,z_fV,time);
% [forv_VI] = for_v(x_fVI,y_fVI,z_fVI,time);
% 
% err_fI=abs((gp*grad_I+cor_I)-forv_I);
% err_fII=abs((gp*grad_II+cor_II)-forv_II);
% err_fIII=abs((gp*grad_III+cor_III)-forv_III);
% err_fIV=abs((gp*grad_IV+cor_IV)-forv_IV);
% err_fV=abs((gp*grad_V+cor_V)-forv_V);
% err_fVI=abs((gp*grad_VI+cor_VI)-forv_VI);
% 
% err_v=max(max([err_fI, err_fII, err_fIII, err_fIV, err_fV, err_fVI]))




for time=1:100
    clc; time
    [ ht_fI,    vt_fI] = sol_exacte(x_fI,   y_fI,   z_fI,   time);
    [ ht_fII,   vt_fII] = sol_exacte(x_fII,  y_fII,  z_fII,  time);
    [ ht_fIII,  vt_fIII] = sol_exacte(x_fIII, y_fIII, z_fIII, time);
    [ ht_fIV,   vt_fIV] = sol_exacte(x_fIV,  y_fIV,  z_fIV,  time);
    [ ht_fV,    vt_fV] = sol_exacte(x_fV,   y_fV,   z_fV,   time);
    [ ht_fVI,   vt_fVI] = sol_exacte(x_fVI,  y_fVI,  z_fVI,  time);


    figure(1)
    plot_cs11(n,nn,ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI);
    hold off
end




