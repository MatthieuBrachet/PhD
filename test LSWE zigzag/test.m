clc; clear all; close all; format short;

global n nn dxi
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr test
global hp gp u0 radius omega

test=3;
opt_ftr=10;
n=40;
mod72

time=10;

%% fonction
[ ht_fI,    vt_fI] = sol_exacte(x_fI,   y_fI,   z_fI,   time);
[ ht_fII,   vt_fII] = sol_exacte(x_fII,  y_fII,  z_fII,  time);
[ ht_fIII,  vt_fIII] = sol_exacte(x_fIII, y_fIII, z_fIII, time);
[ ht_fIV,   vt_fIV] = sol_exacte(x_fIV,  y_fIV,  z_fIV,  time);
[ ht_fV,    vt_fV] = sol_exacte(x_fV,   y_fV,   z_fV,   time);
[ ht_fVI,   vt_fVI] = sol_exacte(x_fVI,  y_fVI,  z_fVI,  time);


%% second membre

% EQUATION 1
[forv_fI] = for_v(x_fI,y_fI,z_fI,time);
[forv_fII] = for_v(x_fII,y_fII,z_fII,time);
[forv_fIII] = for_v(x_fIII,y_fIII,z_fIII,time);
[forv_fIV] = for_v(x_fIV,y_fIV,z_fIV,time);
[forv_fV] = for_v(x_fV,y_fV,z_fV,time);
[forv_fVI] = for_v(x_fVI,y_fVI,z_fVI,time);

% GRADIENT
[grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
    gr72( ht_fI  , ht_fII  , ht_fIII  , ht_fIV  , ht_fV  , ht_fVI  , n , nn);
% CALCUL DE CORIOLIS
[cor_I,cor_II,cor_III,cor_IV,cor_V,cor_VI]=coriolis( vt_fI  , vt_fII  , vt_fIII  , vt_fIV  , vt_fV  , vt_fVI );
    
% ASSEMBLAGE
K1v_fI   = -(gp*grad_I   + cor_I) + forv_fI;
K1v_fII  = -(gp*grad_II  + cor_II) + forv_fII;
K1v_fIII = -(gp*grad_III + cor_III) + forv_fIII;
K1v_fIV  = -(gp*grad_IV  + cor_IV) + forv_fIV;
K1v_fV   = -(gp*grad_V   + cor_V) + forv_fV;
K1v_fVI  = -(gp*grad_VI  + cor_VI) + forv_fVI;
    

err_fI=max(max(max(K1v_fI)));
err_fII=max(max(max(K1v_fII)));
err_fIII=max(max(max(K1v_fIII)));
err_fIV=max(max(max(K1v_fIV)));
err_fV=max(max(max(K1v_fV)));
err_fVI=max(max(max(K1v_fVI)));

err_v=max([err_fI, err_fII, err_fIII, err_fIV, err_fV, err_fVI])

% EQUATION 2
[forh_fI] = for_h(x_fI,y_fI,z_fI,t);
[forh_fII] = for_h(x_fII,y_fII,z_fII,t);
[forh_fIII] = for_h(x_fIII,y_fIII,z_fIII,t);
[forh_fIV] = for_h(x_fIV,y_fIV,z_fIV,t);
[forh_fV] = for_h(x_fV,y_fV,z_fV,t);
[forh_fVI] = for_h(x_fVI,y_fVI,z_fVI,t);
    
% divergence
[div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
    div72( vt_fI  , vt_fII  , vt_fIII  , vt_fIV  , vt_fV  , vt_fVI  , n , nn);
    
% assemblage
K1h_fI   = -hp*div_fI + forh_fI;
K1h_fII  = -hp*div_fII + forh_fII;
K1h_fIII = -hp*div_fIII + forh_fIII;
K1h_fIV  = -hp*div_fIV + forh_fIV;
K1h_fV   = -hp*div_fV + forh_fV;
K1h_fVI  = -hp*div_fVI + forh_fVI;

err_fI=max(max(max(K1h_fI)));
err_fII=max(max(max(K1h_fII)));
err_fIII=max(max(max(K1h_fIII)));
err_fIV=max(max(max(K1h_fIV)));
err_fV=max(max(max(K1h_fV)));
err_fVI=max(max(max(K1h_fVI)));

err_h=max([err_fI, err_fII, err_fIII, err_fIV, err_fV, err_fVI])


