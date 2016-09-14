clc; clear all; close all;

global n nn
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr test
global teta0 teta1

test=0;
%test=0 : test de J. Galewsky stationnaire;
%    = 1 : test de J. Galewsky avec perturbation.
video = 'no';
sauvegarde = 1;
opt_ftr=10;

err_1=[];
err_2=[];
HH=[];

NN=10:10:100;
for i=1:length(NN)
    clc; n=NN(i)
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

    %% *** test stationnaire


    %% EQUATION 1
    % FORCAGE
    [forv_fI] = for_v(x_fI,y_fI,z_fI,t);
    [forv_fII] = for_v(x_fII,y_fII,z_fII,t);
    [forv_fIII] = for_v(x_fIII,y_fIII,z_fIII,t);
    [forv_fIV] = for_v(x_fIV,y_fIV,z_fIV,t);
    [forv_fV] = for_v(x_fV,y_fV,z_fV,t);
    [forv_fVI] = for_v(x_fVI,y_fVI,z_fVI,t);

    % GRADIENT
    [fun_I, fun_II, fun_III, fun_IV, fun_V, fun_VI] = fun_eq1(vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI, ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI);
    [grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
         gr72( fun_I, fun_II, fun_III, fun_IV, fun_V, fun_VI  , n , nn);
    % CALCUL DE CORIOLIS
    [cor_I,cor_II,cor_III,cor_IV,cor_V,cor_VI]=coriolis( vt_fI  , vt_fII  , vt_fIII  , vt_fIV  , vt_fV  , vt_fVI );
    
    % ASSEMBLAGE
    Kv_fI   = -(grad_I   + cor_I) + forv_fI;
    Kv_fII  = -(grad_II  + cor_II) + forv_fII;
    Kv_fIII = -(grad_III + cor_III) + forv_fIII;
    Kv_fIV  = -(grad_IV  + cor_IV) + forv_fIV;
    Kv_fV   = -(grad_V   + cor_V) + forv_fV;
    Kv_fVI  = -(grad_VI  + cor_VI) + forv_fVI;
    
    err=max(max(max(max(abs([Kv_fI, Kv_fII, Kv_fIII, Kv_fIV, Kv_fV, Kv_fVI])))));
    err_1=[err_1 err];

    
    %% EQUATION 2
    [forh_fI] = for_h(x_fI,y_fI,z_fI,t);
    [forh_fII] = for_h(x_fII,y_fII,z_fII,t);
    [forh_fIII] = for_h(x_fIII,y_fIII,z_fIII,t);
    [forh_fIV] = for_h(x_fIV,y_fIV,z_fIV,t);
    [forh_fV] = for_h(x_fV,y_fV,z_fV,t);
    [forh_fVI] = for_h(x_fVI,y_fVI,z_fVI,t);
    
    % divergence
    [vec_I, vec_II, vec_III, vec_IV, vec_V, vec_VI] = fun_eq2(vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI, ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI);
    [div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
        div72( vec_I, vec_II, vec_III, vec_IV, vec_V, vec_VI  , n , nn);
    
    % assemblage
    Kh_fI   = -div_fI + forh_fI;
    Kh_fII  = -div_fII + forh_fII;
    Kh_fIII = -div_fIII + forh_fIII;
    Kh_fIV  = -div_fIV + forh_fIV;
    Kh_fV   = -div_fV + forh_fV;
    Kh_fVI  = -div_fVI + forh_fVI;

    err=max(max(max(abs([Kh_fI, Kh_fII, Kh_fIII, Kh_fIV, Kh_fV, Kh_fVI]))));
    err_2=[err_2 err];
    
    %%
    h=1/(n+1);
    HH=[HH h];
    
end


figure(1)
loglog(HH,err_1,HH, err_2)
legend('eq 1','eq 2')
grid on

figure(2)
subplot(2,2,1)
plot_cs11(n,nn,Kv_fI(:,:,1),Kv_fII(:,:,1),Kv_fIII(:,:,1),Kv_fIV(:,:,1),Kv_fV(:,:,1),Kv_fVI(:,:,1))
title('x-comp , eq1')

subplot(2,2,2)
plot_cs11(n,nn,Kv_fI(:,:,2),Kv_fII(:,:,2),Kv_fIII(:,:,2),Kv_fIV(:,:,2),Kv_fV(:,:,2),Kv_fVI(:,:,2))
title('y-comp , eq1')

subplot(2,2,3)
plot_cs11(n,nn,Kv_fI(:,:,3),Kv_fII(:,:,3),Kv_fIII(:,:,3),Kv_fIV(:,:,3),Kv_fV(:,:,3),Kv_fVI(:,:,3))
title('z-comp , eq1')

subplot(2,2,4)
plot_cs11(n,nn,Kh_fI,Kh_fII,Kh_fIII,Kh_fIV,Kh_fV,Kh_fVI)
title('eq2')

    