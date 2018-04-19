%% ************************************************************************
% Resolution de LSWEC sur la Cubed-Sphere.
%
% *** options :
% test :   test = 0 : solution stationnaire de type Galewsky
%          test = 1 : solution en exp(-sigma*t) (forcage/ammortissement).
%          test = 2 : N. Paldor test case.
%          test = 3 : personnal test case.
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
global opt_ftr test scheme detec
global hp gp u0 radius omega
global teta0 teta1

test=1;
video = 'no';
sauvegarde = 1;
opt_ftr='redonnet10';
type_ftr='symetric';
scheme='compact4';
n=255;
mod101

teta0=-pi/3;
teta1=pi/3;

cgrav=sqrt(gp*hp);
ccor=radius*omega;
c=max(cgrav,ccor);

cfl=0.9;
ddt=radius*dxi*cfl/c;
ndaymax=5;
Tmax=3600*1.5;%ndaymax*3600*24;
itermax=100000;

%% *** initialisation des données
t=0;
[ ht_fI,    vt_fI]   = sol_exacte(x_fI,   y_fI,   z_fI,   t);
[ ht_fII,   vt_fII]  = sol_exacte(x_fII,  y_fII,  z_fII,  t);
[ ht_fIII,  vt_fIII] = sol_exacte(x_fIII, y_fIII, z_fIII, t);
[ ht_fIV,   vt_fIV]  = sol_exacte(x_fIV,  y_fIV,  z_fIV,  t);
[ ht_fV,    vt_fV]   = sol_exacte(x_fV,   y_fV,   z_fV,   t);
[ ht_fVI,   vt_fVI]  = sol_exacte(x_fVI,  y_fVI,  z_fVI,  t);
h0max=max(max(abs([ht_fI ht_fII ht_fIII ht_fIV ht_fV ht_fVI])));
v0max=max(max(max(abs([vt_fI vt_fII vt_fIII vt_fIV vt_fV vt_fVI]))));

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

    %% filtrage
    if strcmp(type_ftr,'adaptatif') == 1
        [detec]=filtre101(na,'redonnet10');
        [det_fI,det_fII,det_fIII,det_fIV,det_fV,det_fVI]=det101(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, n, nn);
        
        %[ftr]=filtre101(na,opt_ftr);
        [htff_fI, htff_fII, htff_fIII, htff_fIV, htff_fV, htff_fVI]=ftr72(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, n, nn);
        [vtff_fI(:,:,1), vtff_fII(:,:,1), vtff_fIII(:,:,1), vtff_fIV(:,:,1), vtff_fV(:,:,1), vtff_fVI(:,:,1)]=ftr72(vt_fI(:,:,1), vt_fII(:,:,1), vt_fIII(:,:,1), vt_fIV(:,:,1), vt_fV(:,:,1), vt_fVI(:,:,1),n,nn);
        [vtff_fI(:,:,2), vtff_fII(:,:,2), vtff_fIII(:,:,2), vtff_fIV(:,:,2), vtff_fV(:,:,2), vtff_fVI(:,:,2)]=ftr72(vt_fI(:,:,2), vt_fII(:,:,2), vt_fIII(:,:,2), vt_fIV(:,:,2), vt_fV(:,:,2), vt_fVI(:,:,2),n,nn);
        [vtff_fI(:,:,3), vtff_fII(:,:,3), vtff_fIII(:,:,3), vtff_fIV(:,:,3), vtff_fV(:,:,3), vtff_fVI(:,:,3)]=ftr72(vt_fI(:,:,3), vt_fII(:,:,3), vt_fIII(:,:,3), vt_fIV(:,:,3), vt_fV(:,:,3), vt_fVI(:,:,3),n,nn);
        
        ht_fI=det_fI.*htff_fI+(1-det_fI).*ht_fI;
        ht_fII=det_fII.*htff_fII+(1-det_fII).*ht_fII;
        ht_fIII=det_fIII.*htff_fIII+(1-det_fIII).*ht_fIII;
        ht_fIV=det_fIV.*htff_fIV+(1-det_fIV).*ht_fIV;
        ht_fV=det_fV.*htff_fV+(1-det_fV).*ht_fV;
        ht_fVI=det_fVI.*htff_fVI+(1-det_fVI).*ht_fVI;
        
        for comp=1:3
            vt_fI(:,:,comp)=det_fI.*vtff_fI(:,:,comp)+(1-det_fI).*vt_fI(:,:,comp);
            vt_fII(:,:,comp)=det_fII.*vtff_fII(:,:,comp)+(1-det_fII).*vt_fII(:,:,comp);
            vt_fIII(:,:,comp)=det_fIII.*vtff_fIII(:,:,comp)+(1-det_fIII).*vt_fIII(:,:,comp);
            vt_fIV(:,:,comp)=det_fIV.*vtff_fIV(:,:,comp)+(1-det_fIV).*vt_fIV(:,:,comp);
            vt_fV(:,:,comp)=det_fV.*vtff_fV(:,:,comp)+(1-det_fV).*vt_fV(:,:,comp);
            vt_fVI(:,:,comp)=det_fVI.*vtff_fVI(:,:,comp)+(1-det_fVI).*vt_fVI(:,:,comp);
        end
    
    elseif strcmp(type_ftr,'nonsymetric') == 1
        [ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI]=ftr72(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, n, nn);
        [vt_fI(:,:,1), vt_fII(:,:,1), vt_fIII(:,:,1), vt_fIV(:,:,1), vt_fV(:,:,1), vt_fVI(:,:,1)]=ftr72(vt_fI(:,:,1), vt_fII(:,:,1), vt_fIII(:,:,1), vt_fIV(:,:,1), vt_fV(:,:,1), vt_fVI(:,:,1),n,nn);
        [vt_fI(:,:,2), vt_fII(:,:,2), vt_fIII(:,:,2), vt_fIV(:,:,2), vt_fV(:,:,2), vt_fVI(:,:,2)]=ftr72(vt_fI(:,:,2), vt_fII(:,:,2), vt_fIII(:,:,2), vt_fIV(:,:,2), vt_fV(:,:,2), vt_fVI(:,:,2),n,nn);
        [vt_fI(:,:,3), vt_fII(:,:,3), vt_fIII(:,:,3), vt_fIV(:,:,3), vt_fV(:,:,3), vt_fVI(:,:,3)]=ftr72(vt_fI(:,:,3), vt_fII(:,:,3), vt_fIII(:,:,3), vt_fIV(:,:,3), vt_fV(:,:,3), vt_fVI(:,:,3),n,nn);
    
    elseif strcmp(type_ftr,'symetric') == 1
        [ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI]=ftr_mixte101(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, n, nn);
        [vt_fI(:,:,1), vt_fII(:,:,1), vt_fIII(:,:,1), vt_fIV(:,:,1), vt_fV(:,:,1), vt_fVI(:,:,1)]=ftr_mixte101(vt_fI(:,:,1), vt_fII(:,:,1), vt_fIII(:,:,1), vt_fIV(:,:,1), vt_fV(:,:,1), vt_fVI(:,:,1),n,nn);
        [vt_fI(:,:,2), vt_fII(:,:,2), vt_fIII(:,:,2), vt_fIV(:,:,2), vt_fV(:,:,2), vt_fVI(:,:,2)]=ftr_mixte101(vt_fI(:,:,2), vt_fII(:,:,2), vt_fIII(:,:,2), vt_fIV(:,:,2), vt_fV(:,:,2), vt_fVI(:,:,2),n,nn);
        [vt_fI(:,:,3), vt_fII(:,:,3), vt_fIII(:,:,3), vt_fIV(:,:,3), vt_fV(:,:,3), vt_fVI(:,:,3)]=ftr_mixte101(vt_fI(:,:,3), vt_fII(:,:,3), vt_fIII(:,:,3), vt_fIV(:,:,3), vt_fV(:,:,3), vt_fVI(:,:,3),n,nn);
    elseif strcmp(type_ftr,'inf')==1
        
    else
        error('Option ''filtre'' is uncorrect.');
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
    err_int(iter)=int./intref-1;
    
    [ E ] = energy(vt_fI,vt_fII,vt_fIII,vt_fIV,vt_fV,vt_fVI,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,n,nn);
    err_energy(iter)=E./Eref-1;
    
    maxi(iter)=max(max([ht_fI; ht_fII; ht_fIII; ht_fIV; ht_fV; ht_fVI]));
    mini(iter)=min(min([ht_fI; ht_fII; ht_fIII; ht_fIV; ht_fV; ht_fVI]));

    %% film
    if strcmp(video,'yes')==1
        hFig = figure(100);
        set(gcf,'PaperPositionMode','auto')
        set(hFig, 'Position', [50 50 1000 500])
        plot_cs100(n,nn,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI)
        caxis([-0.05*hp 0.1*hp])
        title(['numerical solution at ', num2str(time(end)), 'days'])
        hold off;
        mov(iter) = getframe(gcf);
        clf
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
colorbar
if sauvegarde==1
    mkdir(['./RK4_results-' date ]);
    print('-dpng', ['./RK4_results-' date '/ref_' num2str(ref) '_courbe.png'])
    savefig(['./RK4_results-' date '/ref_' num2str(ref) '_courbe']);
    
    save(['./RK4_results-' date '/ref_' num2str(ref) '_erreurdata_test_' num2str(test) '.mat']);
end 

figure(2)
plot_cs11(n,nn,h_fI,h_fII,h_fIII,h_fIV,h_fV,h_fVI);
title(['exat solution at time = ', num2str(time(end))])
colorbar

figure(3)
semilogy(time, erri,'k.', time, err2,'k--', time, err1,'k-')
xlabel('Temps')
ylabel('Erreur relative')
legend('norme infinie','norme 2','norme 1')
grid on
if sauvegarde==1
    print('-dpng', ['./RK4_results-' date '/ref_' num2str(ref) '_erreur.png'])
    savefig(['./RK4_results-' date '/ref_' num2str(ref) '_erreur']);
end 

figure(4)
plot(time,err_int,time,err_energy,'Linewidth',2)
legend('masse','energie')
xlabel('Temps')
grid on
if sauvegarde==1
    mkdir(['./RK4_results-' date ]);
    print('-dpng', ['./RK4_results-' date '/ref_' num2str(ref) '_conservation.png'])
    savefig(['./RK4_results-' date '/ref_' num2str(ref) '_conservation']);
end 

figure(5)
hFig = figure(5);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [50 50 1000 500])
plot_cs100(n,nn,(h_fI-ht_fI)./nrmrefi,(h_fII-ht_fII)./nrmrefi,(h_fIII-ht_fIII)./nrmrefi,(h_fIV-ht_fIV)./nrmrefi,(h_fV-ht_fV)./nrmrefi,(h_fVI-ht_fVI)./nrmrefi);
title('relative error at final time')
colorbar
if sauvegarde==1
    mkdir(['./RK4_results-' date ]);
    print('-dpng', ['./RK4_results-' date '/ref_' num2str(ref) '_space_error.png'])
    savefig(['./RK4_results-' date '/ref_' num2str(ref) '_space_error']);
end 

mm=min(min([ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI]));
MM=max(max([ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI]));

figure(6)
hFig = figure(6);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [50 50 1000 500])
plot_cs101(n,nn,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,linspace(mm,MM,10));
title('relative error at final time')
colorbar
if sauvegarde==1
    mkdir(['./RK4_results-' date ]);
    print('-dpng', ['./RK4_results-' date '/ref_' num2str(ref) '_plot.png'])
    savefig(['./RK4_results-' date '/ref_' num2str(ref) '_plot']);
end 

hFig = figure(7);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [50 50 1000 500])
plot_cs100(n,nn,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI);
title('relative error at final time')
colorbar

fig_placier

format shorte
err1(end)
err2(end)
erri(end)
M=max(max(max(abs([v_fI, v_fII,v_fIII,v_fIV,v_fV,v_fVI]))));
max(max(max(abs([vt_fI-v_fI, vt_fII-v_fII,vt_fIII-v_fIII,vt_fIV-v_fIV,vt_fV-v_fV,vt_fVI-v_fVI]))))./M
err_int(end)
err_energy(end)