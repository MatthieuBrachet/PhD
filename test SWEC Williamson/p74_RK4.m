%% ************************************************************************
% Solve SWEC on the Cubed-Sphere.
% *** options :
% test = 0 : test 2 of Williamson & al.,
%        1 : test 5 of Williamson & al..
%        2 : test 5 of Williamson with smooth mountain,
%        3 : stationnary Galewsky (exp),
%        4 : Galewsky with perturbation (exp),
%        -1 : test with Earth topography
% scheme : numerical spatial scheme used. 
% video : 'yes' ou 'no', do a video or not,
% nper  :  periodicity of frames in the video.
% sauvegarde = 0 (do not save data), 1 (save all data).
% opt_ftr : explicit (redonnet) or adaptative filtering.

%% ************************************************************************
clc; clear all; close all;
format long

global n nn na dxi
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr opt_ftr1 opt_detec test scheme
global gp h0 u0 radius omega
global alpha
global ftr detec
global teta0 teta1

comment='.';
test=1;
video = 'no';
nper=1;
sauvegarde = 0;
filtre='classic';
opt_ftr='redonnet10';
opt_detec='redonnet10';
opt_ftr1='redonnet6';
scheme='compact4';
snapshot='yes';

n=15; % for snapshot, n must be odd !
ndaymax=5;
mod74
disp('mod74 : ok')

ccor=radius*omega;
cgrav=sqrt(h0*gp);
cvit=u0;
c=max([cgrav,ccor,cvit]);

cfl=0.9;
ddt=radius*dxi*cfl/c;
Tmax=ndaymax*3600*24;
itermax=5000;

tstart=cputime;
ref=floor(10000*now);
jour=date;

%% *** test data **********************************************************

if test == -1
    alpha=0;
    u0=20;
    h0=6500;
elseif test == 0
    alpha=pi/7;
    u0=2*pi*radius/(12*24*3600);
    h0=2.94*10^4/gp;
elseif test == 1
    alpha=0;
    u0=20;
    h0=5960;
elseif test == 2
    alpha=0;
    u0=20;
    h0=5960;
elseif test == 3
    alpha=0;
    teta0=pi/7;
    teta1=pi/2-teta0;
    u0=80;
    h0=0;
elseif test == 4
    alpha=0;
    teta0=pi/7;
    teta1=pi/2-teta0;
    u0=80;
    h0=-1.581861685963503e+02+10000;
end

%% *** initial data *******************************************************
t=0;
[ ht_fI,    vt_fI]   = sol_exacte(x_fI,   y_fI,   z_fI,   t);
[ ht_fII,   vt_fII]  = sol_exacte(x_fII,  y_fII,  z_fII,  t);
[ ht_fIII,  vt_fIII] = sol_exacte(x_fIII, y_fIII, z_fIII, t);
[ ht_fIV,   vt_fIV]  = sol_exacte(x_fIV,  y_fIV,  z_fIV,  t);
[ ht_fV,    vt_fV]   = sol_exacte(x_fV,   y_fV,   z_fV,   t);
[ ht_fVI,   vt_fVI]  = sol_exacte(x_fVI,  y_fVI,  z_fVI,  t);

h_fI=ht_fI; v_fI=vt_fI;
h_fII=ht_fII; v_fII=vt_fII;
h_fIII=ht_fIII; v_fIII=vt_fIII;
h_fIV=ht_fIV; v_fIV=vt_fIV;
h_fV=ht_fV; v_fV=vt_fV;
h_fVI=ht_fVI; v_fVI=vt_fVI;

%% *** quantités a conserver **********************************************
[Mref] = mass( ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI );
[Eref] = energy(ht_fI,ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI);
[PEref] = enstrophy(ht_fI,ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI);

iter=1; FTR=0;
time(1)=t; erri(1)=0; err_int(1)=1;
%% *** video option *******************************************************
if strcmp(video,'yes')==1
    mkdir(['./RK4_video-' jour ])
    vidObj=VideoWriter(['./RK4_video-' jour '/ref_' num2str(ref) '.avi']);
    open(vidObj);
    set(gca,'nextplot','replacechildren');
end

%% *** iterations *********************************************************
while t<Tmax && iter<itermax
    %% filtrage
    if strcmp(filtre,'classic') == 1
        [ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI]=ftr_mixte74(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, n, nn);
        [vt_fI(:,:,1), vt_fII(:,:,1), vt_fIII(:,:,1), vt_fIV(:,:,1), vt_fV(:,:,1), vt_fVI(:,:,1)]=ftr_mixte74(vt_fI(:,:,1), vt_fII(:,:,1), vt_fIII(:,:,1), vt_fIV(:,:,1), vt_fV(:,:,1), vt_fVI(:,:,1),n,nn);
        [vt_fI(:,:,2), vt_fII(:,:,2), vt_fIII(:,:,2), vt_fIV(:,:,2), vt_fV(:,:,2), vt_fVI(:,:,2)]=ftr_mixte74(vt_fI(:,:,2), vt_fII(:,:,2), vt_fIII(:,:,2), vt_fIV(:,:,2), vt_fV(:,:,2), vt_fVI(:,:,2),n,nn);
        [vt_fI(:,:,3), vt_fII(:,:,3), vt_fIII(:,:,3), vt_fIV(:,:,3), vt_fV(:,:,3), vt_fVI(:,:,3)]=ftr_mixte74(vt_fI(:,:,3), vt_fII(:,:,3), vt_fIII(:,:,3), vt_fIV(:,:,3), vt_fV(:,:,3), vt_fVI(:,:,3),n,nn);

    elseif strcmp(filtre,'adaptative') == 1
        [detec]=filtre74(na,opt_detec);
        [det_fI,det_fII,det_fIII,det_fIV,det_fV,det_fVI]=det74(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, n, nn);
        
        [ftr]=filtre74(na,opt_ftr1);
        [htff_fI, htff_fII, htff_fIII, htff_fIV, htff_fV, htff_fVI]=ftr_mixte74(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, n, nn);
        [vtff_fI(:,:,1), vtff_fII(:,:,1), vtff_fIII(:,:,1), vtff_fIV(:,:,1), vtff_fV(:,:,1), vtff_fVI(:,:,1)]=ftr_mixte74(vt_fI(:,:,1), vt_fII(:,:,1), vt_fIII(:,:,1), vt_fIV(:,:,1), vt_fV(:,:,1), vt_fVI(:,:,1),n,nn);
        [vtff_fI(:,:,2), vtff_fII(:,:,2), vtff_fIII(:,:,2), vtff_fIV(:,:,2), vtff_fV(:,:,2), vtff_fVI(:,:,2)]=ftr_mixte74(vt_fI(:,:,2), vt_fII(:,:,2), vt_fIII(:,:,2), vt_fIV(:,:,2), vt_fV(:,:,2), vt_fVI(:,:,2),n,nn);
        [vtff_fI(:,:,3), vtff_fII(:,:,3), vtff_fIII(:,:,3), vtff_fIV(:,:,3), vtff_fV(:,:,3), vtff_fVI(:,:,3)]=ftr_mixte74(vt_fI(:,:,3), vt_fII(:,:,3), vt_fIII(:,:,3), vt_fIV(:,:,3), vt_fV(:,:,3), vt_fVI(:,:,3),n,nn);
        
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
    elseif strcmp(filtre,'inf')==1
        
    else
        error('Option ''filtre'' is uncorrect.');
    
    end
   
    %% K1
    hh_fI   = ht_fI;
    hh_fII  = ht_fII;
    hh_fIII = ht_fIII;
    hh_fIV  = ht_fIV;
    hh_fV   = ht_fV;
    hh_fVI  = ht_fVI;
    
    vv_fI   = vt_fI;
    vv_fII  = vt_fII;
    vv_fIII = vt_fIII;
    vv_fIV  = vt_fIV;
    vv_fV   = vt_fV;
    vv_fVI  = vt_fVI;
    
    [K1v_fI, K1v_fII, K1v_fIII, K1v_fIV, K1v_fV, K1v_fVI]=eq_moment74(hh_fI,...
                 hh_fII, hh_fIII, hh_fIV, hh_fV, hh_fVI, vv_fI, vv_fII, vv_fIII, vv_fIV, vv_fV, vv_fVI);
            
    [K1h_fI, K1h_fII, K1h_fIII, K1h_fIV, K1h_fV, K1h_fVI]=eq_cons74(hh_fI,...
                hh_fII, hh_fIII, hh_fIV, hh_fV, hh_fVI, vv_fI, vv_fII, vv_fIII, vv_fIV, vv_fV, vv_fVI);
            
    %% K2
    hh_fI   = ht_fI   + 0.5*ddt*K1h_fI;
    hh_fII  = ht_fII  + 0.5*ddt*K1h_fII;
    hh_fIII = ht_fIII + 0.5*ddt*K1h_fIII;
    hh_fIV  = ht_fIV  + 0.5*ddt*K1h_fIV;
    hh_fV   = ht_fV   + 0.5*ddt*K1h_fV;
    hh_fVI  = ht_fVI  + 0.5*ddt*K1h_fVI;
    
    vv_fI   = vt_fI   + 0.5*ddt*K1v_fI;
    vv_fII  = vt_fII  + 0.5*ddt*K1v_fII;
    vv_fIII = vt_fIII + 0.5*ddt*K1v_fIII;
    vv_fIV  = vt_fIV  + 0.5*ddt*K1v_fIV;
    vv_fV   = vt_fV   + 0.5*ddt*K1v_fV;
    vv_fVI  = vt_fVI  + 0.5*ddt*K1v_fVI;
    
    [K2v_fI, K2v_fII, K2v_fIII, K2v_fIV, K2v_fV, K2v_fVI]=eq_moment74(hh_fI,...
                hh_fII, hh_fIII, hh_fIV, hh_fV, hh_fVI, vv_fI, vv_fII, vv_fIII, vv_fIV, vv_fV, vv_fVI);
            
    [K2h_fI, K2h_fII, K2h_fIII, K2h_fIV, K2h_fV, K2h_fVI]=eq_cons74(hh_fI,...
                 hh_fII, hh_fIII, hh_fIV, hh_fV, hh_fVI, vv_fI, vv_fII, vv_fIII, vv_fIV, vv_fV, vv_fVI);
     
    %% K3
    hh_fI   = ht_fI   + 0.5*ddt*K2h_fI;
    hh_fII  = ht_fII  + 0.5*ddt*K2h_fII;
    hh_fIII = ht_fIII + 0.5*ddt*K2h_fIII;
    hh_fIV  = ht_fIV  + 0.5*ddt*K2h_fIV;
    hh_fV   = ht_fV   + 0.5*ddt*K2h_fV;
    hh_fVI  = ht_fVI  + 0.5*ddt*K2h_fVI;
    
    vv_fI   = vt_fI   + 0.5*ddt*K2v_fI;
    vv_fII  = vt_fII  + 0.5*ddt*K2v_fII;
    vv_fIII = vt_fIII + 0.5*ddt*K2v_fIII;
    vv_fIV  = vt_fIV  + 0.5*ddt*K2v_fIV;
    vv_fV   = vt_fV   + 0.5*ddt*K2v_fV;
    vv_fVI  = vt_fVI  + 0.5*ddt*K2v_fVI;
    
    [K3v_fI, K3v_fII, K3v_fIII, K3v_fIV, K3v_fV, K3v_fVI]=eq_moment74(hh_fI,...
                 hh_fII, hh_fIII, hh_fIV, hh_fV, hh_fVI, vv_fI, vv_fII, vv_fIII, vv_fIV, vv_fV, vv_fVI);
            
    [K3h_fI, K3h_fII, K3h_fIII, K3h_fIV, K3h_fV, K3h_fVI]=eq_cons74(hh_fI,...
                hh_fII, hh_fIII, hh_fIV, hh_fV, hh_fVI, vv_fI, vv_fII, vv_fIII, vv_fIV, vv_fV, vv_fVI);
            
     %% K4
    hh_fI   = ht_fI   + ddt*K3h_fI;
    hh_fII  = ht_fII  + ddt*K3h_fII;
    hh_fIII = ht_fIII + ddt*K3h_fIII;
    hh_fIV  = ht_fIV  + ddt*K3h_fIV;
    hh_fV   = ht_fV   + ddt*K3h_fV;
    hh_fVI  = ht_fVI  + ddt*K3h_fVI;
    
    vv_fI   = vt_fI   + ddt*K3v_fI;
    vv_fII  = vt_fII  + ddt*K3v_fII;
    vv_fIII = vt_fIII + ddt*K3v_fIII;
    vv_fIV  = vt_fIV  + ddt*K3v_fIV;
    vv_fV   = vt_fV   + ddt*K3v_fV;
    vv_fVI  = vt_fVI  + ddt*K3v_fVI;
    
    [K4v_fI, K4v_fII, K4v_fIII, K4v_fIV, K4v_fV, K4v_fVI]=eq_moment74(hh_fI,...
                 hh_fII, hh_fIII, hh_fIV, hh_fV, hh_fVI, vv_fI, vv_fII, vv_fIII, vv_fIV, vv_fV, vv_fVI);
            
    [K4h_fI, K4h_fII, K4h_fIII, K4h_fIV, K4h_fV, K4h_fVI]=eq_cons74(hh_fI,...
                 hh_fII, hh_fIII, hh_fIV, hh_fV, hh_fVI, vv_fI, vv_fII, vv_fIII, vv_fIV, vv_fV, vv_fVI);       
            
            
    %% End data
    Sh_fI   = (K1h_fI   + 2*K2h_fI   + 2*K3h_fI   + K4h_fI);
    Sh_fII  = (K1h_fII  + 2*K2h_fII  + 2*K3h_fII  + K4h_fII);
    Sh_fIII = (K1h_fIII + 2*K2h_fIII + 2*K3h_fIII + K4h_fIII);
    Sh_fIV  = (K1h_fIV  + 2*K2h_fIV  + 2*K3h_fIV  + K4h_fIV);
    Sh_fV   = (K1h_fV   + 2*K2h_fV   + 2*K3h_fV   + K4h_fV);
    Sh_fVI  = (K1h_fVI  + 2*K2h_fVI  + 2*K3h_fVI  + K4h_fVI);
    
    Sv_fI   = (K1v_fI   + 2*K2v_fI   + 2*K3v_fI   + K4v_fI);
    Sv_fII  = (K1v_fII  + 2*K2v_fII  + 2*K3v_fII  + K4v_fII);
    Sv_fIII = (K1v_fIII + 2*K2v_fIII + 2*K3v_fIII + K4v_fIII);
    Sv_fIV  = (K1v_fIV  + 2*K2v_fIV  + 2*K3v_fIV  + K4v_fIV);
    Sv_fV   = (K1v_fV   + 2*K2v_fV   + 2*K3v_fV   + K4v_fV);
    Sv_fVI  = (K1v_fVI  + 2*K2v_fVI  + 2*K3v_fVI  + K4v_fVI);
    
    
    htnew_fI    = ht_fI   + ddt/6 * Sh_fI;
    htnew_fII   = ht_fII  + ddt/6 * Sh_fII;
    htnew_fIII  = ht_fIII + ddt/6 * Sh_fIII;
    htnew_fIV   = ht_fIV  + ddt/6 * Sh_fIV;
    htnew_fV    = ht_fV   + ddt/6 * Sh_fV;
    htnew_fVI   = ht_fVI  + ddt/6 * Sh_fVI;

    vtnew_fI    = vt_fI   + ddt/6 * Sv_fI;
    vtnew_fII   = vt_fII  + ddt/6 * Sv_fII;
    vtnew_fIII  = vt_fIII + ddt/6 * Sv_fIII;
    vtnew_fIV   = vt_fIV  + ddt/6 * Sv_fIV;
    vtnew_fV    = vt_fV   + ddt/6 * Sv_fV;
    vtnew_fVI   = vt_fVI  + ddt/6 * Sv_fVI;


    %% error on h
    t=t+ddt;
    
    err_fI   = htnew_fI  -h_fI;
    err_fII  = htnew_fII -h_fII;
    err_fIII = htnew_fIII-h_fIII;
    err_fIV  = htnew_fIV -h_fIV;
    err_fV   = htnew_fV  -h_fV ;
    err_fVI  = htnew_fVI -h_fVI;
    
    errv_fI   = vtnew_fI   - v_fI;
    errv_fII  = vtnew_fII  - v_fII;
    errv_fIII = vtnew_fIII - v_fIII;
    errv_fIV  = vtnew_fIV  - v_fIV;
    errv_fV   = vtnew_fV   - v_fV ;
    errv_fVI  = vtnew_fVI  - v_fVI;
    
    str='infty';
        [~,~,~,~,~,~,nrmger]=...
      nrm74(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
    [~,~,~,~,~,~,nrmref]=nrm74(h_fI,h_fII,h_fIII,h_fIV,h_fV,h_fVI,n,nn,str);
    erri(iter)=nrmger./nrmref;
    
    str='2';
        [~,~,~,~,~,~,nrmger]=...
      nrm74(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
    [~,~,~,~,~,~,nrmref]=nrm74(h_fI,h_fII,h_fIII,h_fIV,h_fV,h_fVI,n,nn,str);
    err2(iter)=nrmger./nrmref;
    
    str='1';
        [~,~,~,~,~,~,nrmger]=...
      nrm74(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
    [~,~,~,~,~,~,nrmref]=nrm74(h_fI,h_fII,h_fIII,h_fIV,h_fV,h_fVI,n,nn,str);
    err1(iter)=nrmger./nrmref;
    
    %% stabilization
    stab_fI   = htnew_fI   - ht_fI;
    stab_fII  = htnew_fII  - ht_fII;
    stab_fIII = htnew_fIII - ht_fIII;
    stab_fIV  = htnew_fIV  - ht_fIV;
    stab_fV   = htnew_fV   - ht_fV ;
    stab_fVI  = htnew_fVI  - ht_fVI;

    str='infty';
    [~,~,~,~,~,~,nrmger]=...
      nrm74(stab_fI,stab_fII,stab_fIII,stab_fIV,stab_fV,stab_fVI,n,nn,str);
    stabi(iter)=nrmger;
    
    time(iter)=t/(24*3600);
    
    %% update height and velocity
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
    [M] = mass( ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI );
    err_int(iter)=M/Mref;
    [E] = energy(ht_fI,ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI);
    err_energy(iter)=E/Eref;
    [PE] = enstrophy(ht_fI,ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI);
    err_enstrophy(iter)=PE/PEref;
    
    %% video
    if strcmp(video,'yes')==1 & mod(iter,nper) == 0
        
        [vort_fI,vort_fII,vort_fIII,vort_fIV,vort_fV,vort_fVI]=...
            vort74(vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI,n,nn);
        
        clf
        figure(9)
        hFig = figure(9);
        set(gcf,'PaperPositionMode','auto')
        set(hFig, 'Position', [50 50 1000 500])
        plot_cs100(n,nn,vort_fI,vort_fII,vort_fIII,vort_fIV,vort_fV,vort_fVI);
        title(['vorticity at time : ', num2str(time(end))])
        hold off

        currFrame = getframe;
        writeVideo(vidObj,currFrame);
    end
    
    % snapshot
    if sauvegarde == 1 & mod(iter,floor(Tmax/(3*ddt))) == 0 & strcmp(snapshot,'yes')==1 
        mkdir(['./RK4_results-' jour '/' num2str(ref)])
        close all;
        
        [vort_fI,vort_fII,vort_fIII,vort_fIV,vort_fV,vort_fVI]=...
            vort74(vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI,n,nn);
        
        pas=.2*10^-4;
        mmin=min(min([vort_fI vort_fII vort_fIII vort_fIV vort_fV vort_fVI]));
        mmax=max(max([vort_fI vort_fII vort_fIII vort_fIV vort_fV vort_fVI]));
        v=mmin:pas:mmax;
 
        figure(100)
        hFig = figure(100);
        set(gcf,'PaperPositionMode','auto')
        set(hFig, 'Position', [50 50 1000 500])
        if test == 4
            plot_cs101(n,nn,vort_fI,vort_fII,vort_fIII,vort_fIV,vort_fV,vort_fVI,v);
            title(['vorticity at time : ', num2str(time(end))])
        elseif test == 1
            plot_cs101(n,nn,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,5050:50:5950);
            title(['solution at time : ', num2str(time(end))])
        else
            plot_cs100(n,nn,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI);
            title(['solution at time : ', num2str(time(end))])
        end

        print('-dpng', ['./RK4_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_snapshot_intermediaire' num2str(floor(100*time(end))) '.png'])
        savefig(['./RK4_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_snapshot_intermediaire_' num2str(floor(100*time(end))) '.fig']);
    end

    %% historique sur la divergence
    [div_fI, div_fII, div_fIII, div_fIV, div_fV, div_fVI]=...
        div74(vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI,n,nn);
    
    Mdivu(iter)=max(max([div_fI div_fII div_fIII div_fIV div_fV div_fVI]));
    
    iter=iter+1;
    clc; 
    disp(real([test iter min(itermax,floor(Tmax/ddt)) erri(end) err_enstrophy(end)]));
end
disp('calcul : OK');
disp('plot...');
tend=cputime-tstart;

if strcmp(video,'yes') == 1
    close(vidObj);
    
    fid = fopen('AA_VIDEO_SAVE_RK4.txt','a');
    fprintf(fid,'%s\n',['date : ', jour]);
    fprintf(fid,'%s\n',['ref. : ', num2str(ref)]);
    fprintf(fid,'%s\n','***********************************');
    fprintf(fid,'%s\n',['test : ', num2str(test)]);
    fprintf(fid,'%s\n','---------- numerical data ---------');
    fprintf(fid,'%s\n',['number of points  : ', num2str(n)] );
    fprintf(fid,'%s\n',['time step         : ', num2str(ddt)] );
    fprintf(fid,'%s\n',['cfl               : ', num2str(cfl)] );
    fprintf(fid,'%s\n',['type du filtre    : ',  filtre] );
    fprintf(fid,'%s\n',['ordre du filtre p.: '  opt_ftr] );
    fprintf(fid,'%s\n',['ordre du filtre a.: '  opt_ftr1] );
    fprintf(fid,'%s\n',['ordre du detect.  : '  opt_detec] );
    fprintf(fid,'%s\n',['scheme            : '  scheme] );
    fprintf(fid,'%s\n','---------- physical data ----------');
    fprintf(fid,'%s\n',['gravity g              : ', num2str(gp)] );
    fprintf(fid,'%s\n',['alpha                  : ', num2str(alpha)] );
    fprintf(fid,'%s\n',['caracteristiv velocity : ', num2str(u0)] );
    fprintf(fid,'%s\n',['coriolis parameter     : ', num2str(omega)] );
    fprintf(fid,'%s\n','------------ comment --------------');
    fprintf(fid,'%s\n', comment );
    fprintf(fid,'%s\n','***********************************');
    fprintf(fid,'%s\n','  ');
    fprintf(fid,'%s\n','  ');
    fclose(fid);
end

if sauvegarde == 1
    fid = fopen('AA_RESULTS_SAVE_RK4.txt','a');
    fprintf(fid,'%s\n',['date : ', jour]);
    fprintf(fid,'%s\n',['ref. : ', num2str(ref)]);
    fprintf(fid,'%s\n','***********************************');
    fprintf(fid,'%s\n',['test : ', num2str(test)]);
    fprintf(fid,'%s\n','---------- numerical data ---------');
    fprintf(fid,'%s\n',['number of points  : ', num2str(n)] );
    fprintf(fid,'%s\n',['time step         : ', num2str(ddt)] );
    fprintf(fid,'%s\n',['cfl               : ', num2str(cfl)] );
    fprintf(fid,'%s\n',['type du filtre    : ',  filtre] );
    fprintf(fid,'%s\n',['ordre du filtre s.: '  opt_ftr] );
    fprintf(fid,'%s\n',['ordre du filtre a.: '  opt_ftr1] );
    fprintf(fid,'%s\n',['ordre du detect.  : '  opt_detec] );
    fprintf(fid,'%s\n',['scheme            : '  scheme] );
    fprintf(fid,'%s\n','---------- physical data ----------');
    fprintf(fid,'%s\n',['gravity g              : ', num2str(gp)] );
    fprintf(fid,'%s\n',['alpha                  : ', num2str(alpha)] );
    fprintf(fid,'%s\n',['caracteristiv velocity : ', num2str(u0)] );
    fprintf(fid,'%s\n',['coriolis parameter     : ', num2str(omega)] );
    fprintf(fid,'%s\n','------------ comment --------------');
    fprintf(fid,'%s\n', comment );
    fprintf(fid,'%s\n','***********************************');
    fprintf(fid,'%s\n','  ');
    fprintf(fid,'%s\n','  ');
    fclose(fid);
end
 

figure(1)
plot_cs11(n,nn,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI);
title(['calculated solution at time = ', num2str(time(end))])
if sauvegarde==1
    mkdir(['./RK4_results-' jour '/' num2str(ref) ])
    print('-dpng', ['./RK4_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_courbe.png'])
    savefig(['./RK4_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_courbe']);
    save(['./RK4_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_erreurdata_test_' num2str(test) '.mat']);
end 

figure(2)
plot_cs11(n,nn,h_fI,h_fII,h_fIII,h_fIV,h_fV,h_fVI);
title(['exact solution at time = ', num2str(time(end))])

figure(3)
plot_cs11(n,nn,err_fI, err_fII, err_fIII, err_fIV, err_fV, err_fVI);
title(['error at time = ', num2str(time(end))])

if strcmp(snapshot,'yes')==1
    [vort_fI,vort_fII,vort_fIII,vort_fIV,vort_fV,vort_fVI]=...
            vort74(vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI,n,nn);
    
    figure(4)
    hFig = figure(4);
    set(gcf,'PaperPositionMode','auto')
    set(hFig, 'Position', [50 50 1000 500])
    plot_cs100(n,nn,vort_fI,vort_fII,vort_fIII,vort_fIV,vort_fV,vort_fVI)
    title(['vorticity at time : ', num2str(time(end))])
    colorbar
    if sauvegarde == 1
        print('-dpng', ['./RK4_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_snapshot.png'])
        savefig(['./RK4_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_snapshot']);
    end
    
    pas=.2*10^-4;
    mmin=min(min([vort_fI vort_fII vort_fIII vort_fIV vort_fV vort_fVI]));
    mmax=max(max([vort_fI vort_fII vort_fIII vort_fIV vort_fV vort_fVI]));
    v=mmin:pas:mmax;

    figure(5)
    hFig = figure(5);
    set(gcf,'PaperPositionMode','auto')
    set(hFig, 'Position', [50 50 1000 500])
    plot_cs101(n,nn,vort_fI,vort_fII,vort_fIII,vort_fIV,vort_fV,vort_fVI,v);
    title(['vorticity at time : ', num2str(time(end))])
end

figure(6)
semilogy(time, erri,time, err2, time, err1)
xlabel('time')
ylabel('relative error')
legend('infty norm','norm 2','norm 1')
grid on
if sauvegarde==1
    print('-dpng', ['./RK4_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_erreur.png'])
    savefig(['./RK4_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_erreur']);
end 

figure(7)
plot(time,err_int-1,time,err_energy-1,time,err_enstrophy-1)
legend('mass','energy','potential enstrophy')
xlabel('time')
title('relative quantity')
grid on
if sauvegarde==1
    print('-dpng', ['./RK4_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_conservation.png'])
    savefig(['./RK4_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_conservation']);
end 

figure(8)
semilogy(time,stabi)
xlabel('time')
title('stabilization of height')
grid on

figure(9)
hFig = figure(9);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [50 50 1000 500])
plot_cs100(n,nn,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI);
title(['calculated solution at time = ', num2str(time(end))])
colorbar
if sauvegarde==1
    print('-dpng', ['./RK4_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_solution.png'])
    savefig(['./RK4_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_solution']);
end


figure(10)
hFig = figure(10);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [50 50 1000 500])
plot_cs101(n,nn,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,5050:50:5950);
title(['calculated solution at time = ', num2str(time(end))])
if sauvegarde==1
    print('-dpng', ['./RK4_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_snapshot_solution.png'])
    savefig(['./RK4_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_snapshot_solution']);
end 

[div_fI, div_fII, div_fIII, div_fIV, div_fV, div_fVI]=...
        div74(vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI,n,nn);
    
figure(11)
hFig = figure(11);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [50 50 1000 500])
plot_cs100(n,nn,div_fI, div_fII, div_fIII, div_fIV, div_fV, div_fVI)
title(['divergence at time : ', num2str(time(end))])
colorbar
if sauvegarde==1
    print('-dpng', ['./RK4_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_divergence.png'])
    savefig(['./RK4_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_divergence']);
end 

figure(12)
plot(time,Mdivu)

[detec]=filtre74(na,opt_detec);
[det_fI,det_fII,det_fIII,det_fIV,det_fV,det_fVI]=det74(ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI, n, nn);

figure(13)
plot_cs11(n,nn,det_fI,det_fII,det_fIII,det_fIV,det_fV,det_fVI)
title(['detection oscillations : ', num2str(time(end))])
if sauvegarde==1
    print('-dpng', ['./RK4_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_gibbs.png'])
    savefig(['./RK4_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_gibbs']);
end 

fig_placier
disp(['temps de calcul (sans les graphiques) : ', num2str(tend)])