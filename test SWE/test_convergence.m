%% test de convergence
clc; clear all; close all
global n nn
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr test
global gp
global teta0 teta1

test = 0;
opt_ftr=0;
NN=[10:10:50];
teta0=-3*pi/16;
teta1=3*pi/16;
ERRV=[];
ERRH=[];
for i=1:length(NN)

    clc; n=NN(i)
    mod72;

    %% *** initialisation des données
    t=0;
    [ ht_fI,    vt_fI] = sol_exacte(x_fI,   y_fI,   z_fI,   t);
    [ ht_fII,   vt_fII] = sol_exacte(x_fII,  y_fII,  z_fII,  t);
    [ ht_fIII,  vt_fIII] = sol_exacte(x_fIII, y_fIII, z_fIII, t);
    [ ht_fIV,   vt_fIV] = sol_exacte(x_fIV,  y_fIV,  z_fIV,  t);
    [ ht_fV,    vt_fV] = sol_exacte(x_fV,   y_fV,   z_fV,   t);
    [ ht_fVI,   vt_fVI] = sol_exacte(x_fVI,  y_fVI,  z_fVI,  t);
    
    %% K1
    
    hh_fI=ht_fI;
    hh_fII=ht_fII;
    hh_fIII=ht_fIII;
    hh_fIV=ht_fIV;
    hh_fV=ht_fV;
    hh_fVI=ht_fVI;
    
    vv_fI=vt_fI;
    vv_fII=vt_fII;
    vv_fIII=vt_fIII;
    vv_fIV=vt_fIV;
    vv_fV=vt_fV;
    vv_fVI=vt_fVI;
    
    % premiere equation
    
    % FORCAGE
    [forv_fI] = for_v(x_fI,y_fI,z_fI,t);
    [forv_fII] = for_v(x_fII,y_fII,z_fII,t);
    [forv_fIII] = for_v(x_fIII,y_fIII,z_fIII,t);
    [forv_fIV] = for_v(x_fIV,y_fIV,z_fIV,t);
    [forv_fV] = for_v(x_fV,y_fV,z_fV,t);
    [forv_fVI] = for_v(x_fVI,y_fVI,z_fVI,t);

    % GRADIENT
    [grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]= gr72( hh_fI, hh_fII, hh_fIII, hh_fIV, hh_fV, hh_fVI  , n , nn);
    % CALCUL DE CORIOLIS
    [cor_I,cor_II,cor_III,cor_IV,cor_V,cor_VI]=coriolis( vv_fI, vv_fII, vv_fIII, vv_fIV, vv_fV, vv_fVI );
    % ADVECTION
    [adv_fI,adv_fII,adv_fIII,adv_fIV,adv_fV,adv_fVI]=adv73(vv_fI, vv_fII, vv_fIII, vv_fIV, vv_fV, vv_fVI ,n ,nn );
    
    % ASSEMBLAGE
    K1v_fI   = -(gp*grad_I   + cor_I   + adv_fI  ) + forv_fI;
    K1v_fII  = -(gp*grad_II  + cor_II  + adv_fII ) + forv_fII;
    K1v_fIII = -(gp*grad_III + cor_III + adv_fIII) + forv_fIII;
    K1v_fIV  = -(gp*grad_IV  + cor_IV  + adv_fIV ) + forv_fIV;
    K1v_fV   = -(gp*grad_V   + cor_V   + adv_fV  ) + forv_fV;
    K1v_fVI  = -(gp*grad_VI  + cor_VI  + adv_fVI ) + forv_fVI;
    
    err_vt=max(max(max(max([K1v_fI, K1v_fII, K1v_fIII, K1v_fIV, K1v_fV, K1v_fVI]))));
    ERRV=[ERRV err_vt];
    
    % seconde equation
    
    [forh_fI] = for_h(x_fI,y_fI,z_fI,t);
    [forh_fII] = for_h(x_fII,y_fII,z_fII,t);
    [forh_fIII] = for_h(x_fIII,y_fIII,z_fIII,t);
    [forh_fIV] = for_h(x_fIV,y_fIV,z_fIV,t);
    [forh_fV] = for_h(x_fV,y_fV,z_fV,t);
    [forh_fVI] = for_h(x_fVI,y_fVI,z_fVI,t);
    
    % divergence
    [vec_I, vec_II, vec_III, vec_IV, vec_V, vec_VI] = fun_eq2(vv_fI, vv_fII, vv_fIII, vv_fIV, vv_fV, vv_fVI, hh_fI, hh_fII, hh_fIII, hh_fIV, hh_fV, hh_fVI);

    
    [div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
        div72( vec_I, vec_II, vec_III, vec_IV, vec_V, vec_VI,n,nn);
    
    % assemblage
     K1h_fI   = -div_fI + forh_fI;
     K1h_fII  = -div_fII + forh_fII;
     K1h_fIII = -div_fIII + forh_fIII;
     K1h_fIV  = -div_fIV + forh_fIV;
     K1h_fV   = -div_fV + forh_fV;
     K1h_fVI  = -div_fVI + forh_fVI;
     
     err_ht=max(max(max([K1h_fI, K1h_fII, K1h_fIII, K1h_fIV, K1h_fV, K1h_fVI])));
     ERRH=[ERRH err_ht];
     
end


figure(1)
semilogy(NN,ERRV,NN,ERRH)
legend('vt','ht')