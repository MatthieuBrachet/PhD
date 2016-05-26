clc; clear all; close all; format short;

global n nn
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr test
global hp gp radius dxi u0

test=1;
opt_ftr=0;
disc=[10:10:200];

e_vit=[];
e_ht=[];
for i=1:length(disc)

    clc; n=disc(i)
    mod72
    
    cfl=0.5;
    ddt=radius*dxi*cfl/u0;

    %% *** initialisation des données
    [ vt_fI, ht_fI ] = fun_init(x_fI,y_fI,z_fI);
    [ vt_fII, ht_fII ] = fun_init(x_fII,y_fII,z_fII);
    [ vt_fIII, ht_fIII ] = fun_init(x_fIII,y_fIII,z_fIII);
    [ vt_fIV, ht_fIV ] = fun_init(x_fIV,y_fIV,z_fIV);
    [ vt_fV, ht_fV ] = fun_init(x_fV,y_fV,z_fV);
    [ vt_fVI, ht_fVI ] = fun_init(x_fVI,y_fVI,z_fVI);

    %% K1
    % premiere equation

    % GRADIENT
    [grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
        gr72(ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,n,nn);
    % CALCUL DE CORIOLIS
    [cor_I,cor_II,cor_III,cor_IV,cor_V,cor_VI]=coriolis(vt_fI,vt_fII,vt_fIII,vt_fIV,vt_fV,vt_fVI);
    
    % ASSEMBLAGE
    K1v_fI=-ddt*(gp*grad_I+cor_I);
    K1v_fII=-ddt*(gp*grad_II+cor_II);
    K1v_fIII=-ddt*(gp*grad_III+cor_III);
    K1v_fIV=-ddt*(gp*grad_IV+cor_IV);
    K1v_fV=-ddt*(gp*grad_V+cor_V);
    K1v_fVI=-ddt*(gp*grad_VI+cor_VI);
    
    
    
    % seconde equation
    % divergence
    [div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
        div72(vt_fI,vt_fII,vt_fIII,vt_fIV,vt_fV,vt_fVI,n,nn);
    
    % assemblage
    K1h_fI=-ddt*hp*div_fI;
    K1h_fII=-ddt*hp*div_fII;
    K1h_fIII=-ddt*hp*div_fIII;
    K1h_fIV=-ddt*hp*div_fIV;
    K1h_fV=-ddt*hp*div_fV;
    K1h_fVI=-ddt*hp*div_fVI;
    
    %% ERREUR
    % erreur eq en vitesse
    ev_I=max(max(max(abs(K1v_fI))));
    ev_II=max(max(max(abs(K1v_fII))));
    ev_III=max(max(max(abs(K1v_fIII))));
    ev_IV=max(max(max(abs(K1v_fIV))));
    ev_V=max(max(max(abs(K1v_fV))));
    ev_VI=max(max(max(abs(K1v_fVI))));
    
    ev=max([ev_I,ev_II,ev_III,ev_IV,ev_V,ev_VI]);
    
    % erreur eq en hauteur
    eh_I=max(max(max(abs(K1h_fI))));
    eh_II=max(max(max(abs(K1h_fII))));
    eh_III=max(max(max(abs(K1h_fIII))));
    eh_IV=max(max(max(abs(K1h_fIV))));
    eh_V=max(max(max(abs(K1h_fV))));
    eh_VI=max(max(max(abs(K1h_fVI))));
    
    eh=max([eh_I,eh_II,eh_III,eh_IV,eh_V,eh_VI]);
    
    e_vit=[e_vit, ev];
    e_ht=[e_ht, eh];
end

    figure(1)
    plot_cs11(n,nn,K1v_fI(:,:,1),K1v_fII(:,:,1),K1v_fIII(:,:,1),K1v_fIV(:,:,1),K1v_fV(:,:,1),K1v_fVI(:,:,1));
    
    figure(2)
    plot_cs11(n,nn,K1v_fI(:,:,2),K1v_fII(:,:,2),K1v_fIII(:,:,2),K1v_fIV(:,:,2),K1v_fV(:,:,2),K1v_fVI(:,:,2));
    
    figure(3)
    plot_cs11(n,nn,K1v_fI(:,:,3),K1v_fII(:,:,3),K1v_fIII(:,:,3),K1v_fIV(:,:,3),K1v_fV(:,:,3),K1v_fVI(:,:,3));

    figure(4)
    plot_cs11(n,nn,K1h_fI,K1h_fII,K1h_fIII,K1h_fIV,K1h_fV,K1h_fVI);
    
    figure(5)
    loglog(1./disc,e_vit,1./disc,e_ht, 1./disc, 1./disc.^4)
    legend('erreur vitesse stationnaire','erreur hauteur stationnaire','order 4',2)
    