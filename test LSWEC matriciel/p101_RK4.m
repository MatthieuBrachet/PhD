%% ************************************************************************
% Resolution de LSWEC sur la Cubed-Sphere.
%
% *** options :
% test : 0, ..., 5 choix du test à lancer.
%          test = 0 : solution stationnaire de type Galewsky
%          test = 1 : solution en exp(-sigma*t) (forcage/ammortissement).
%          test = 2 : N. Paldor test case.
% video : 'yes' ou 'no', faire une video ou non.
% sauvegarde = 0 (ne rien sauvegarder), 1 (sauvegarder toutes les valeurs
%          finales).
% opt_ftr : filtre explicite de S. Redonnet (=2, 4, 6, 8, 10, ordre du
%          filtre) (=0, pas de filtrage)
% type_ftr : classic (on the xi and eta variables) or caracteristic 
%            variables.
%
%% ************************************************************************
clc; clear all; close all; format long;
global n nn dxi
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr test scheme
global hp gp u0 radius omega
global teta0 teta1

test=2;
video = 'no';
sauvegarde = 0;
opt_ftr='redonnet10';
type_ftr='classic';
scheme='compact4';
n=31;
mod101

teta0=-pi/5;
teta1=pi/3;

cgrav=sqrt(gp*hp);
ccor=radius*omega;
c=max(cgrav,ccor);

cfl=0.9;
ddt=radius*dxi*cfl/c;
ndaymax=6;
Tmax=ndaymax*3600*24;
itermax=1;

%% *** initialisation des données
t=0;
[ ht_fI,    vt_fI] = sol_exacte(x_fI,   y_fI,   z_fI,   t);
[ ht_fII,   vt_fII] = sol_exacte(x_fII,  y_fII,  z_fII,  t);
[ ht_fIII,  vt_fIII] = sol_exacte(x_fIII, y_fIII, z_fIII, t);
[ ht_fIV,   vt_fIV] = sol_exacte(x_fIV,  y_fIV,  z_fIV,  t);
[ ht_fV,    vt_fV] = sol_exacte(x_fV,   y_fV,   z_fV,   t);
[ ht_fVI,   vt_fVI] = sol_exacte(x_fVI,  y_fVI,  z_fVI,  t);

%% quantités a conserver
[aaa,aaa,aaa,aaa,aaa,aaa,intref]=nrm101(ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,n,nn,'cor_int');
[aaa,aaa,aaa,aaa,aaa,aaa,nrmrefi]=nrm101(ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,n,nn,'infty');
[aaa,aaa,aaa,aaa,aaa,aaa,nrmref2]=nrm101(ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,n,nn,'2');
[aaa,aaa,aaa,aaa,aaa,aaa,nrmref1]=nrm101(ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,n,nn,'1');


[ Eref ] = energy(vt_fI,vt_fII,vt_fIII,vt_fIV,vt_fV,vt_fVI,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,n,nn);
%% *** video option *******************************************************
if strcmp(video,'yes')==1
    nFrames = min(itermax,floor(Tmax/ddt));
    mov(1:nFrames) = struct('cdata', [],'colormap', []);
    set(gca,'nextplot','replacechildren');
end

iter=0;
time(1)=t; erri(1)=0;
while t<Tmax && iter<itermax
    iter=iter+1;
    clc; disp([iter min(itermax,floor(Tmax/ddt)) (erri(end))]);

    %% Filtrage
    if strcmp(type_ftr,'classic')==1
        [ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI]=ftr_mixte101(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, n, nn);
        [vt_fI(:,:,1), vt_fII(:,:,1), vt_fIII(:,:,1), vt_fIV(:,:,1), vt_fV(:,:,1), vt_fVI(:,:,1)]=ftr_mixte101(vt_fI(:,:,1), vt_fII(:,:,1), vt_fIII(:,:,1), vt_fIV(:,:,1), vt_fV(:,:,1), vt_fVI(:,:,1),n,nn);
        [vt_fI(:,:,2), vt_fII(:,:,2), vt_fIII(:,:,2), vt_fIV(:,:,2), vt_fV(:,:,2), vt_fVI(:,:,2)]=ftr_mixte101(vt_fI(:,:,2), vt_fII(:,:,2), vt_fIII(:,:,2), vt_fIV(:,:,2), vt_fV(:,:,2), vt_fVI(:,:,2),n,nn);
        [vt_fI(:,:,3), vt_fII(:,:,3), vt_fIII(:,:,3), vt_fIV(:,:,3), vt_fV(:,:,3), vt_fVI(:,:,3)]=ftr_mixte101(vt_fI(:,:,3), vt_fII(:,:,3), vt_fIII(:,:,3), vt_fIV(:,:,3), vt_fV(:,:,3), vt_fVI(:,:,3),n,nn);
    else
        
    end
    %% K1
    % premiere equation
    
    % FORCAGE
    [forv_fI] = for_v(x_fI,y_fI,z_fI,t);
    [forv_fII] = for_v(x_fII,y_fII,z_fII,t);
    [forv_fIII] = for_v(x_fIII,y_fIII,z_fIII,t);
    [forv_fIV] = for_v(x_fIV,y_fIV,z_fIV,t);
    [forv_fV] = for_v(x_fV,y_fV,z_fV,t);
    [forv_fVI] = for_v(x_fVI,y_fVI,z_fVI,t);

    % GRADIENT
    [grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
        gr101( ht_fI  , ht_fII  , ht_fIII  , ht_fIV  , ht_fV  , ht_fVI  , n , nn);
    % CALCUL DE CORIOLIS
    [cor_I,cor_II,cor_III,cor_IV,cor_V,cor_VI]=coriolis( vt_fI  , vt_fII  , vt_fIII  , vt_fIV  , vt_fV  , vt_fVI );
    
    % ASSEMBLAGE
    K1v_fI   = -(gp*grad_I   + cor_I) + forv_fI;
    K1v_fII  = -(gp*grad_II  + cor_II) + forv_fII;
    K1v_fIII = -(gp*grad_III + cor_III) + forv_fIII;
    K1v_fIV  = -(gp*grad_IV  + cor_IV) + forv_fIV;
    K1v_fV   = -(gp*grad_V   + cor_V) + forv_fV;
    K1v_fVI  = -(gp*grad_VI  + cor_VI) + forv_fVI;
    
    % seconde equation
    
    [forh_fI] = for_h(x_fI,y_fI,z_fI,t);
    [forh_fII] = for_h(x_fII,y_fII,z_fII,t);
    [forh_fIII] = for_h(x_fIII,y_fIII,z_fIII,t);
    [forh_fIV] = for_h(x_fIV,y_fIV,z_fIV,t);
    [forh_fV] = for_h(x_fV,y_fV,z_fV,t);
    [forh_fVI] = for_h(x_fVI,y_fVI,z_fVI,t);
    
    % divergence
    [div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
        div101( vt_fI  , vt_fII  , vt_fIII  , vt_fIV  , vt_fV  , vt_fVI  , n , nn);
    
    err_div(iter)=max(max(abs([div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI])));
    
    % assemblage
    K1h_fI   = -hp*div_fI + forh_fI;
    K1h_fII  = -hp*div_fII + forh_fII;
    K1h_fIII = -hp*div_fIII + forh_fIII;
    K1h_fIV  = -hp*div_fIV + forh_fIV;
    K1h_fV   = -hp*div_fV + forh_fV;
    K1h_fVI  = -hp*div_fVI + forh_fVI;

    %% K2
    % premiere equation

    % FORCAGE
    [forv_fI] = for_v(x_fI,y_fI,z_fI,t+ddt/2);
    [forv_fII] = for_v(x_fII,y_fII,z_fII,t+ddt/2);
    [forv_fIII] = for_v(x_fIII,y_fIII,z_fIII,t+ddt/2);
    [forv_fIV] = for_v(x_fIV,y_fIV,z_fIV,t+ddt/2);
    [forv_fV] = for_v(x_fV,y_fV,z_fV,t+ddt/2);
    [forv_fVI] = for_v(x_fVI,y_fVI,z_fVI,t+ddt/2);
    
    % GRADIENT
    [grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
        gr101( ht_fI+ddt/2*K1h_fI  , ht_fII+ddt/2*K1h_fII  , ht_fIII+ddt/2*K1h_fIII  , ht_fIV+ddt/2*K1h_fIV  , ht_fV+ddt/2*K1h_fV  , ht_fVI+ddt/2*K1h_fVI , n , nn);
    % CALCUL DE CORIOLIS
    [cor_I,cor_II,cor_III,cor_IV,cor_V,cor_VI]=coriolis( vt_fI+ddt/2*K1v_fI , vt_fII+ddt/2*K1v_fII , vt_fIII+ddt/2*K1v_fIII , vt_fIV+ddt/2*K1v_fIV , vt_fV+ddt/2*K1v_fV , vt_fVI+ddt/2*K1v_fVI );
    
    % ASSEMBLAGE
    K2v_fI   = -(gp*grad_I   + cor_I) + forv_fI;
    K2v_fII  = -(gp*grad_II  + cor_II) + forv_fII;
    K2v_fIII = -(gp*grad_III + cor_III) + forv_fIII;
    K2v_fIV  = -(gp*grad_IV  + cor_IV) + forv_fIV;
    K2v_fV   = -(gp*grad_V   + cor_V) + forv_fV;
    K2v_fVI  = -(gp*grad_VI  + cor_VI) + forv_fVI;
    
    % seconde equation
    
    [forh_fI] = for_h(x_fI,y_fI,z_fI,t);
    [forh_fII] = for_h(x_fII,y_fII,z_fII,t);
    [forh_fIII] = for_h(x_fIII,y_fIII,z_fIII,t);
    [forh_fIV] = for_h(x_fIV,y_fIV,z_fIV,t);
    [forh_fV] = for_h(x_fV,y_fV,z_fV,t);
    [forh_fVI] = for_h(x_fVI,y_fVI,z_fVI,t);
    
    % divergence
    [div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
        div101( vt_fI+ddt/2*K1v_fI , vt_fII+ddt/2*K1v_fII , vt_fIII+ddt/2*K1v_fIII , vt_fIV+ddt/2*K1v_fIV , vt_fV+ddt/2*K1v_fV , vt_fVI+ddt/2*K1v_fVI , n , nn);
    
    % assemblage
    K2h_fI   = -hp*div_fI + forh_fI;
    K2h_fII  = -hp*div_fII + forh_fII;
    K2h_fIII = -hp*div_fIII + forh_fIII;
    K2h_fIV  = -hp*div_fIV + forh_fIV;
    K2h_fV   = -hp*div_fV + forh_fV;
    K2h_fVI  = -hp*div_fVI + forh_fVI;
    
    %% K3
    % premiere equation

    % FORCAGE
    [forv_fI] = for_v(x_fI,y_fI,z_fI,t+ddt/2);
    [forv_fII] = for_v(x_fII,y_fII,z_fII,t+ddt/2);
    [forv_fIII] = for_v(x_fIII,y_fIII,z_fIII,t+ddt/2);
    [forv_fIV] = for_v(x_fIV,y_fIV,z_fIV,t+ddt/2);
    [forv_fV] = for_v(x_fV,y_fV,z_fV,t+ddt/2);
    [forv_fVI] = for_v(x_fVI,y_fVI,z_fVI,t+ddt/2);
    
    % GRADIENT
    [grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
        gr101( ht_fI+ddt/2*K2h_fI , ht_fII+ddt/2*K2h_fII , ht_fIII+ddt/2*K2h_fIII , ht_fIV+ddt/2*K2h_fIV , ht_fV+ddt/2*K2h_fV , ht_fVI+ddt/2*K2h_fVI , n , nn);
    % CALCUL DE CORIOLIS
    [cor_I,cor_II,cor_III,cor_IV,cor_V,cor_VI]=coriolis( vt_fI+ddt/2*K2v_fI , vt_fII+ddt/2*K2v_fII , vt_fIII+ddt/2*K2v_fIII , vt_fIV+ddt/2*K2v_fIV , vt_fV+ddt/2*K2v_fV , vt_fVI+ddt/2*K2v_fVI);
    
    % ASSEMBLAGE
    K3v_fI   = -(gp*grad_I   + cor_I) + forv_fI;
    K3v_fII  = -(gp*grad_II  + cor_II) + forv_fII;
    K3v_fIII = -(gp*grad_III + cor_III) + forv_fIII;
    K3v_fIV  = -(gp*grad_IV  + cor_IV) + forv_fIV;
    K3v_fV   = -(gp*grad_V   + cor_V) + forv_fV;
    K3v_fVI  = -(gp*grad_VI  + cor_VI) + forv_fVI;
    
    % seconde equation
    
    [forh_fI] = for_h(x_fI,y_fI,z_fI,t+ddt/2);
    [forh_fII] = for_h(x_fII,y_fII,z_fII,t+ddt/2);
    [forh_fIII] = for_h(x_fIII,y_fIII,z_fIII,t+ddt/2);
    [forh_fIV] = for_h(x_fIV,y_fIV,z_fIV,t+ddt/2);
    [forh_fV] = for_h(x_fV,y_fV,z_fV,t+ddt/2);
    [forh_fVI] = for_h(x_fVI,y_fVI,z_fVI,t+ddt/2);
    
    % divergence
    [div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
        div101( vt_fI+ddt/2*K2v_fI , vt_fII+ddt/2*K2v_fII , vt_fIII+ddt/2*K2v_fIII , vt_fIV+ddt/2*K2v_fIV , vt_fV+ddt/2*K2v_fV , vt_fVI+ddt/2*K2v_fVI , n , nn);
    
    % assemblage
    K3h_fI   = -hp*div_fI + forh_fI;
    K3h_fII  = -hp*div_fII + forh_fII;
    K3h_fIII = -hp*div_fIII + forh_fIII;
    K3h_fIV  = -hp*div_fIV + forh_fIV;
    K3h_fV   = -hp*div_fV + forh_fV;
    K3h_fVI  = -hp*div_fVI + forh_fVI;
    
    %% K4
    % premiere equation

    [forv_fI] = for_v(x_fI,y_fI,z_fI,t+ddt);
    [forv_fII] = for_v(x_fII,y_fII,z_fII,t+ddt);
    [forv_fIII] = for_v(x_fIII,y_fIII,z_fIII,t+ddt);
    [forv_fIV] = for_v(x_fIV,y_fIV,z_fIV,t+ddt);
    [forv_fV] = for_v(x_fV,y_fV,z_fV,t+ddt);
    [forv_fVI] = for_v(x_fVI,y_fVI,z_fVI,t+ddt);
    
    % GRADIENT
    [grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
        gr101( ht_fI+ddt*K3h_fI , ht_fII+ddt*K3h_fII , ht_fIII+ddt*K3h_fIII , ht_fIV+ddt*K3h_fIV , ht_fV+ddt*K3h_fV , ht_fVI+ddt*K3h_fVI , n , nn);
    % CALCUL DE CORIOLIS
    [cor_I,cor_II,cor_III,cor_IV,cor_V,cor_VI]=coriolis( vt_fI+ddt*K3v_fI , vt_fII+ddt*K3v_fII , vt_fIII+ddt*K3v_fIII , vt_fIV+ddt*K3v_fIV , vt_fV+ddt*K3v_fV , vt_fVI+ddt*K3v_fVI);
    
    % ASSEMBLAGE
    K4v_fI   = -(gp*grad_I   + cor_I) + forv_fI;
    K4v_fII  = -(gp*grad_II  + cor_II) + forv_fII;
    K4v_fIII = -(gp*grad_III + cor_III) + forv_fIII;
    K4v_fIV  = -(gp*grad_IV  + cor_IV) + forv_fIV;
    K4v_fV   = -(gp*grad_V   + cor_V) + forv_fV;
    K4v_fVI  = -(gp*grad_VI  + cor_VI) + forv_fVI;
    
    % seconde equation
    
    [forh_fI] = for_h(x_fI,y_fI,z_fI,t+ddt);
    [forh_fII] = for_h(x_fII,y_fII,z_fII,t+ddt);
    [forh_fIII] = for_h(x_fIII,y_fIII,z_fIII,t+ddt);
    [forh_fIV] = for_h(x_fIV,y_fIV,z_fIV,t+ddt);
    [forh_fV] = for_h(x_fV,y_fV,z_fV,t+ddt);
    [forh_fVI] = for_h(x_fVI,y_fVI,z_fVI,t+ddt);
    
    % divergence
    [div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
        div101( vt_fI+ddt*K3v_fI , vt_fII+ddt*K3v_fII , vt_fIII+ddt*K3v_fIII , vt_fIV+ddt*K3v_fIV , vt_fV+ddt*K3v_fV , vt_fVI+ddt*K3v_fVI , n , nn );
    
    % assemblage
    K4h_fI    = -hp*div_fI +forh_fI;
    K4h_fII   = -hp*div_fII + forh_fII;
    K4h_fIII  = -hp*div_fIII + forh_fIII;
    K4h_fIV   = -hp*div_fIV + forh_fIV;
    K4h_fV    = -hp*div_fV + forh_fV;
    K4h_fVI   = -hp*div_fVI + forh_fVI;

    %% Assemblage final

    htnew_fI    = ht_fI   + ddt/6 * (K1h_fI   + 2*K2h_fI   + 2*K3h_fI   + K4h_fI);
    htnew_fII   = ht_fII  + ddt/6 * (K1h_fII  + 2*K2h_fII  + 2*K3h_fII  + K4h_fII);
    htnew_fIII  = ht_fIII + ddt/6 * (K1h_fIII + 2*K2h_fIII + 2*K3h_fIII + K4h_fIII);
    htnew_fIV   = ht_fIV  + ddt/6 * (K1h_fIV  + 2*K2h_fIV  + 2*K3h_fIV  + K4h_fIV);
    htnew_fV    = ht_fV   + ddt/6 * (K1h_fV   + 2*K2h_fV   + 2*K3h_fV   + K4h_fV);
    htnew_fVI   = ht_fVI  + ddt/6 * (K1h_fVI  + 2*K2h_fVI  + 2*K3h_fVI  + K4h_fVI);


    vtnew_fI    = vt_fI   + ddt/6 * (K1v_fI   + 2*K2v_fI   + 2*K3v_fI   + K4v_fI);
    vtnew_fII   = vt_fII  + ddt/6 * (K1v_fII  + 2*K2v_fII  + 2*K3v_fII  + K4v_fII);
    vtnew_fIII  = vt_fIII + ddt/6 * (K1v_fIII + 2*K2v_fIII + 2*K3v_fIII + K4v_fIII);
    vtnew_fIV   = vt_fIV  + ddt/6 * (K1v_fIV  + 2*K2v_fIV  + 2*K3v_fIV  + K4v_fIV);
    vtnew_fV    = vt_fV   + ddt/6 * (K1v_fV   + 2*K2v_fV   + 2*K3v_fV   + K4v_fV);
    vtnew_fVI   = vt_fVI  + ddt/6 * (K1v_fVI  + 2*K2v_fVI  + 2*K3v_fVI  + K4v_fVI);


    %% calcul de l'erreur
    t=t+ddt;
    
    [h_fI,v_fI] = sol_exacte(x_fI,y_fI,z_fI,t);
    [h_fII,v_fII] = sol_exacte(x_fII,y_fII,z_fII,t);
    [h_fIII,v_fIII] = sol_exacte(x_fIII,y_fIII,z_fIII,t);
    [h_fIV,v_fIV] = sol_exacte(x_fIV,y_fIV,z_fIV,t);
    [h_fV,v_fV] = sol_exacte(x_fV,y_fV,z_fV,t);
    [h_fVI,v_fVI] = sol_exacte(x_fVI,y_fVI,z_fVI,t);
    
    
    err_fI=htnew_fI-h_fI;
    err_fII=htnew_fII-h_fII;
    err_fIII=htnew_fIII-h_fIII;
    err_fIV=htnew_fIV-h_fIV;
    err_fV=htnew_fV-h_fV ;
    err_fVI=htnew_fVI-h_fVI;
    
    
    str='infty';
        [aaa,aaa,aaa,aaa,aaa,aaa,nrmger]=...
      nrm101(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
    erri(iter)=nrmger./nrmrefi;
    str='2';
        [aaa,aaa,aaa,aaa,aaa,aaa,nrmger]=...
      nrm101(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
    err2(iter)=nrmger./nrmref2;
    str='1';
        [aaa,aaa,aaa,aaa,aaa,aaa,nrmger]=...
      nrm101(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
    err1(iter)=nrmger./nrmref1;
    
    time(iter)=t/(24*3600);
    
    %% mise a jour
    vt_fI   = vtnew_fI;
    vt_fII  = vtnew_fII;
    vt_fIII = vtnew_fIII;
    vt_fIV  = vtnew_fIV;
    vt_fV   = vtnew_fV;
    vt_fVI  = vtnew_fVI;
    
    ht_fI   = htnew_fI;
    ht_fII  = htnew_fII;
    ht_fIII = htnew_fIII;
    ht_fIV  = htnew_fIV;
    ht_fV   = htnew_fV;
    ht_fVI  = htnew_fVI;
    
    %% conservation
    str='cor_int';
    [aaa,aaa,aaa,aaa,aaa,aaa,int]=...
    nrm101(ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,n,nn,str);
    err_int(iter)=int./intref;
    
    [ E ] = energy(vt_fI,vt_fII,vt_fIII,vt_fIV,vt_fV,vt_fVI,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,n,nn);
    err_energy(iter)=E./Eref;
    
    maxi(iter)=max(max([ht_fI; ht_fII; ht_fIII; ht_fIV; ht_fV; ht_fVI]));
    mini(iter)=min(min([ht_fI; ht_fII; ht_fIII; ht_fIV; ht_fV; ht_fVI]));

    %% film
    if strcmp(video,'yes')==1
        figure(100)
        plot_cs11(n,nn,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI)
        title(['numerical solution at ', num2str(time(end)), 'days'])
        hold off;
        mov(iter) = getframe(gcf);
    end

end
ref=floor(10000*now);
if strcmp(video,'yes') == 1
    mkdir(['./RK4_video-' date ])
    movie2avi(mov, ['./RK4_video-' date '/ref_' num2str(ref) '.avi'], 'compression', 'None');
    
    fid = fopen('AA_VIDEO_SAVE_RK4.txt','a');
    fprintf(fid,'%s\n',['date : ', date]);
    fprintf(fid,'%s\n',['ref. : ', num2str(ref)]);
    fprintf(fid,'%s\n','***********************************');
    fprintf(fid,'%s\n',['test : ', num2str(test)]);
    fprintf(fid,'%s\n','---------- numerical data ---------');
    fprintf(fid,'%s\n',['number of points  : ', num2str(n)] );
    fprintf(fid,'%s\n',['time step         : ', num2str(ddt)] );
    fprintf(fid,'%s\n',['cfl               : ', num2str(cfl)] );
    fprintf(fid,'%s\n',['ordre du filtre   : ', opt_ftr] );
    fprintf(fid,'%s\n',['schema            : ', scheme] );
    fprintf(fid,'%s\n',['type du filtre    : ' type_ftr] );
    fprintf(fid,'%s\n','---------- physical data ----------');
    fprintf(fid,'%s\n',['gravity g              : ', num2str(gp)] );
    fprintf(fid,'%s\n',['equilibrium SWE hp     : ', num2str(hp)] );
    fprintf(fid,'%s\n',['caracteristiv velocity : ', num2str(u0)] );
    fprintf(fid,'%s\n',['coriolis parameter     : ', num2str(omega)] );
    fprintf(fid,'%s\n','***********************************');
    fprintf(fid,'%s\n','  ');
    fprintf(fid,'%s\n','  ');
    fclose(fid);
end

if sauvegarde == 1
    fid = fopen('AA_RESULTS_SAVE_RK4.txt','a');
    fprintf(fid,'%s\n',['date : ', date]);
    fprintf(fid,'%s\n',['ref. : ', num2str(ref)]);
    fprintf(fid,'%s\n','***********************************');
    fprintf(fid,'%s\n',['test : ', num2str(test)]);
    fprintf(fid,'%s\n','---------- numerical data ---------');
    fprintf(fid,'%s\n',['number of points  : ', num2str(n)] );
    fprintf(fid,'%s\n',['time step         : ', num2str(ddt)] );
    fprintf(fid,'%s\n',['cfl               : ', num2str(cfl)] );
    fprintf(fid,'%s\n',['ordre du filtre   : ', opt_ftr] );
    fprintf(fid,'%s\n',['schema            : ', scheme] );
    fprintf(fid,'%s\n',['type du filtre    : ' type_ftr] );
    fprintf(fid,'%s\n','---------- physical data ----------');
    fprintf(fid,'%s\n',['gravity g              : ', num2str(gp)] );
    fprintf(fid,'%s\n',['equilibrium SWE hp     : ', num2str(hp)] );
    fprintf(fid,'%s\n',['caracteristiv velocity : ', num2str(u0)] );
    fprintf(fid,'%s\n',['coriolis parameter     : ', num2str(omega)] );
    fprintf(fid,'%s\n','***********************************');
    fprintf(fid,'%s\n','  ');
    fprintf(fid,'%s\n','  ');
    fclose(fid);
end
    
figure(1)
plot_cs11(n,nn,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI);
title(['calculated solution at time = ', num2str(time(end))])
if sauvegarde==1
    mkdir(['./RK4_results-' date ]);
    print('-dpng', ['./RK4_results-' date '/ref_' num2str(ref) '_courbe.png'])
    savefig(['./RK4_results-' date '/ref_' num2str(ref) '_courbe']);
    
    save(['./RK4_results-' date '/ref_' num2str(ref) '_erreurdata_test_' num2str(test) '.mat']);
end 

figure(2)
plot_cs11(n,nn,h_fI,h_fII,h_fIII,h_fIV,h_fV,h_fVI);
title(['exat solution at time = ', num2str(time(end))])

figure(3)
semilogy(time, erri, time, err2, time, err1)
xlabel('time')
ylabel('relative error')
legend('infty norm','norm 2','norm 1')
grid on
if sauvegarde==1
    print('-dpng', ['./RK4_results-' date '/ref_' num2str(ref) '_erreur.png'])
    savefig(['./RK4_results-' date '/ref_' num2str(ref) '_erreur']);
end 

figure(4)
semilogy(time,err_int-1,time,err_energy-1)
legend('mass','energy')
xlabel('time')
title('error on conservation')
grid on
if sauvegarde==1
    mkdir(['./RK4_results-' date ]);
    print('-dpng', ['./RK4_results-' date '/ref_' num2str(ref) '_conservation.png'])
    savefig(['./RK4_results-' date '/ref_' num2str(ref) '_conservation']);
end 

figure(5)
plot_cs11(n,nn,(h_fI-ht_fI)./nrmrefi,(h_fII-ht_fII)./nrmrefi,(h_fIII-ht_fIII)./nrmrefi,(h_fIV-ht_fIV)./nrmrefi,(h_fV-ht_fV)./nrmrefi,(h_fVI-ht_fVI)./nrmrefi);
title('relative error at final time')
if sauvegarde==1
    mkdir(['./RK4_results-' date ]);
    print('-dpng', ['./RK4_results-' date '/ref_' num2str(ref) '_space_error.png'])
    savefig(['./RK4_results-' date '/ref_' num2str(ref) '_space_error']);
end 

figure(6)
plot(time, maxi, time, mini)
title('extremums')
legend('maxi','mini')
ylabel('time')

[ vlambda_fI ,vteta_fI] = project(vt_fI,x_fI,y_fI,z_fI);
[ vlambda_fII ,vteta_fII] = project(vt_fII,x_fII,y_fII,z_fII);
[ vlambda_fIII ,vteta_fIII] = project(vt_fIII,x_fIII,y_fIII,z_fIII);
[ vlambda_fIV ,vteta_fIV] = project(vt_fIV,x_fIV,y_fIV,z_fIV);
[ vlambda_fV ,vteta_fV] = project(vt_fV,x_fV,y_fV,z_fV);
[ vlambda_fVI ,vteta_fVI] = project(vt_fVI,x_fVI,y_fVI,z_fVI);

figure(7)
plot_cs100(n,nn,vlambda_fI,vlambda_fII,vlambda_fIII,vlambda_fIV,vlambda_fV,vlambda_fVI);
title(['numerical projection of vt at time = ', num2str(time(end))])

fig_placier