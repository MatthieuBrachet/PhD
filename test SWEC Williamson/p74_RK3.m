%% ************************************************************************
% Resolution de SWEC sur la Cubed-Sphere.
% *** options :
% test = 0 : test 2 of Williamson & al.,
%        1 : test 5 of Williamson & al..
%        2 : test 5 of Williamson with smooth mountain,
% scheme : numerical spatial scheme used. 
% video : 'yes' ou 'no', do a video or not,
% nper  :  periodicity of frames in the video.
% sauvegarde = 0 (do not save data), 1 (save all data).
% opt_ftr : explicit (redonnet) or implicit (visbal) filtering.
% delta_ftr : is between 0 and 1. If the filter of u is note Fu, then, the
%             filtering action is :
%                        (1-delta_ftr)*u + delta_ftr*Fu
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

comment='.';
test=1;
video = 'yes';
nper=1;
sauvegarde = 1;
filtre='classic';
opt_ftr='redonnet4';
opt_detec='redonnet10';
opt_ftr1='redonnet4';

scheme='compact4';
snapshot='yes';

n=15; % for snapshot, n must be in the form 2^m-1 !
ndaymax=5;
mod74

ccor=radius*omega;
cgrav=sqrt(h0*gp);
cvit=u0;
c=max([cgrav,ccor,cvit]);

cfl=0.5;
ddt=radius*dxi*cfl/c;
Tmax=ndaymax*3600*24;
itermax=10000;

tstart=cputime;
ref=floor(10000*now);
jour=date;

%% *** test data **********************************************************

if test == 0
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
end

%% *** initial data *******************************************************
t=0;
[ ht_fI,    vt_fI]   = sol_exacte(x_fI,   y_fI,   z_fI,   t);
[ ht_fII,   vt_fII]  = sol_exacte(x_fII,  y_fII,  z_fII,  t);
[ ht_fIII,  vt_fIII] = sol_exacte(x_fIII, y_fIII, z_fIII, t);
[ ht_fIV,   vt_fIV]  = sol_exacte(x_fIV,  y_fIV,  z_fIV,  t);
[ ht_fV,    vt_fV]   = sol_exacte(x_fV,   y_fV,   z_fV,   t);
[ ht_fVI,   vt_fVI]  = sol_exacte(x_fVI,  y_fVI,  z_fVI,  t);

%% *** quantités a conserver **********************************************
[Mref] = mass( ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI );
[Eref] = energy(ht_fI,ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI);
[PEref] = enstrophy(ht_fI,ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI);

iter=0; FTR=0;
time(1)=t; erri(1)=0; err_int(1)=1;
%% *** video option *******************************************************
if strcmp(video,'yes')==1
    mkdir(['./RK3_video-' jour ])
    vidObj=VideoWriter(['./RK3_video-' jour '/ref_' num2str(ref) '.avi']);
    open(vidObj);
    set(gca,'nextplot','replacechildren');
end

%% *** iterations *********************************************************
while t<Tmax && iter<itermax
    iter=iter+1;
    clc; 
    disp([iter min(itermax,floor(Tmax/ddt)) erri(end) err_int(end)]);
    
    %% filtrage
    if strcmp(filtre,'classic') == 1
        [ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI]=ftr74(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, n, nn);
        [vt_fI(:,:,1), vt_fII(:,:,1), vt_fIII(:,:,1), vt_fIV(:,:,1), vt_fV(:,:,1), vt_fVI(:,:,1)]=ftr74(vt_fI(:,:,1), vt_fII(:,:,1), vt_fIII(:,:,1), vt_fIV(:,:,1), vt_fV(:,:,1), vt_fVI(:,:,1),n,nn);
        [vt_fI(:,:,2), vt_fII(:,:,2), vt_fIII(:,:,2), vt_fIV(:,:,2), vt_fV(:,:,2), vt_fVI(:,:,2)]=ftr74(vt_fI(:,:,2), vt_fII(:,:,2), vt_fIII(:,:,2), vt_fIV(:,:,2), vt_fV(:,:,2), vt_fVI(:,:,2),n,nn);
        [vt_fI(:,:,3), vt_fII(:,:,3), vt_fIII(:,:,3), vt_fIV(:,:,3), vt_fV(:,:,3), vt_fVI(:,:,3)]=ftr74(vt_fI(:,:,3), vt_fII(:,:,3), vt_fIII(:,:,3), vt_fIV(:,:,3), vt_fV(:,:,3), vt_fVI(:,:,3),n,nn);

    elseif strcmp(filtre,'adaptative') == 1
        [ftr]=filtre74(na,opt_ftr);
        [htf_fI, htf_fII, htf_fIII, htf_fIV, htf_fV, htf_fVI]=ftr74(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, n, nn);
        [vtf_fI(:,:,1), vtf_fII(:,:,1), vtf_fIII(:,:,1), vtf_fIV(:,:,1), vtf_fV(:,:,1), vtf_fVI(:,:,1)]=ftr74(vt_fI(:,:,1), vt_fII(:,:,1), vt_fIII(:,:,1), vt_fIV(:,:,1), vt_fV(:,:,1), vt_fVI(:,:,1),n,nn);
        [vtf_fI(:,:,2), vtf_fII(:,:,2), vtf_fIII(:,:,2), vtf_fIV(:,:,2), vtf_fV(:,:,2), vtf_fVI(:,:,2)]=ftr74(vt_fI(:,:,2), vt_fII(:,:,2), vt_fIII(:,:,2), vt_fIV(:,:,2), vt_fV(:,:,2), vt_fVI(:,:,2),n,nn);
        [vtf_fI(:,:,3), vtf_fII(:,:,3), vtf_fIII(:,:,3), vtf_fIV(:,:,3), vtf_fV(:,:,3), vtf_fVI(:,:,3)]=ftr74(vt_fI(:,:,3), vt_fII(:,:,3), vt_fIII(:,:,3), vt_fIV(:,:,3), vt_fV(:,:,3), vt_fVI(:,:,3),n,nn);
        
        [detec]=filtre74(na,opt_detec);
        [det_fI,det_fII,det_fIII,det_fIV,det_fV,det_fVI]=det74(htf_fI, htf_fII, htf_fIII, htf_fIV, htf_fV, htf_fVI, n, nn);
        
        [ftr]=filtre74(na,opt_ftr1);
        [htff_fI, htff_fII, htff_fIII, htff_fIV, htff_fV, htff_fVI]=ftr74(htf_fI, htf_fII, htf_fIII, htf_fIV, htf_fV, htf_fVI, n, nn);
        [vtff_fI(:,:,1), vtff_fII(:,:,1), vtff_fIII(:,:,1), vtff_fIV(:,:,1), vtff_fV(:,:,1), vtff_fVI(:,:,1)]=ftr74(vtf_fI(:,:,1), vtf_fII(:,:,1), vtf_fIII(:,:,1), vtf_fIV(:,:,1), vtf_fV(:,:,1), vtf_fVI(:,:,1),n,nn);
        [vtff_fI(:,:,2), vtff_fII(:,:,2), vtff_fIII(:,:,2), vtff_fIV(:,:,2), vtff_fV(:,:,2), vtff_fVI(:,:,2)]=ftr74(vtf_fI(:,:,2), vtf_fII(:,:,2), vtf_fIII(:,:,2), vtf_fIV(:,:,2), vtf_fV(:,:,2), vtf_fVI(:,:,2),n,nn);
        [vtff_fI(:,:,3), vtff_fII(:,:,3), vtff_fIII(:,:,3), vtff_fIV(:,:,3), vtff_fV(:,:,3), vtff_fVI(:,:,3)]=ftr74(vtf_fI(:,:,3), vtf_fII(:,:,3), vtf_fIII(:,:,3), vtf_fIV(:,:,3), vtf_fV(:,:,3), vtf_fVI(:,:,3),n,nn);
        
        ht_fI=det_fI.*htff_fI+(1-det_fI).*htf_fI;
        ht_fII=det_fII.*htff_fII+(1-det_fII).*htf_fII;
        ht_fIII=det_fIII.*htff_fIII+(1-det_fIII).*htf_fIII;
        ht_fIV=det_fIV.*htff_fIV+(1-det_fIV).*htf_fIV;
        ht_fV=det_fV.*htff_fV+(1-det_fV).*htf_fV;
        ht_fVI=det_fVI.*htff_fVI+(1-det_fVI).*htf_fVI;
        
        for comp=1:3
            vt_fI(:,:,comp)=det_fI.*vtff_fI(:,:,comp)+(1-det_fI).*vtf_fI(:,:,comp);
            vt_fII(:,:,comp)=det_fII.*vtff_fII(:,:,comp)+(1-det_fII).*vtf_fII(:,:,comp);
            vt_fIII(:,:,comp)=det_fIII.*vtff_fIII(:,:,comp)+(1-det_fIII).*vtf_fIII(:,:,comp);
            vt_fIV(:,:,comp)=det_fIV.*vtff_fIV(:,:,comp)+(1-det_fIV).*vtf_fIV(:,:,comp);
            vt_fV(:,:,comp)=det_fV.*vtff_fV(:,:,comp)+(1-det_fV).*vtf_fV(:,:,comp);
            vt_fVI(:,:,comp)=det_fVI.*vtff_fVI(:,:,comp)+(1-det_fVI).*vtf_fVI(:,:,comp);
        end
        
        
    else
        error('Option ''filtre'' is uncorrect.');
    
    end
   
    %% step 1
    hI=ht_fI;
    hII=ht_fII;
    hIII=ht_fIII;
    hIV=ht_fIV;
    hV=ht_fV;
    hVI=ht_fVI;
    
    vI=vt_fI;
    vII=vt_fII;
    vIII=vt_fIII;
    vIV=vt_fIV;
    vV=vt_fV;
    vVI=vt_fVI;
    
    [K1v_fI, K1v_fII, K1v_fIII, K1v_fIV, K1v_fV, K1v_fVI]=eq_moment74(...
        hI,hII,hIII,hIV,hV,hVI,vI,vII,vIII,vIV,vV,vVI);
            
    [K1h_fI, K1h_fII, K1h_fIII, K1h_fIV, K1h_fV, K1h_fVI]=eq_cons74(...
        hI,hII,hIII,hIV,hV,hVI,vI,vII,vIII,vIV,vV,vVI);
    
    hh1_fI=ht_fI+ddt*K1h_fI;
    hh1_fII=ht_fII+ddt*K1h_fII;
    hh1_fIII=ht_fIII+ddt*K1h_fIII;
    hh1_fIV=ht_fIV+ddt*K1h_fIV;
    hh1_fV=ht_fV+ddt*K1h_fV;
    hh1_fVI=ht_fVI+ddt*K1h_fVI;
    
    vv1_fI=vt_fI+ddt*K1v_fI;
    vv1_fII=vt_fII+ddt*K1v_fII;
    vv1_fIII=vt_fIII+ddt*K1v_fIII;
    vv1_fIV=vt_fIV+ddt*K1v_fIV;
    vv1_fV=vt_fV+ddt*K1v_fV;
    vv1_fVI=vt_fVI+ddt*K1v_fVI;
    
    %% step 2
    hI=hh1_fI;
    hII=hh1_fII;
    hIII=hh1_fIII;
    hIV=hh1_fIV;
    hV=hh1_fV;
    hVI=hh1_fVI;
    
    vI=vv1_fI;
    vII=vv1_fII;
    vIII=vv1_fIII;
    vIV=vv1_fIV;
    vV=vv1_fV;
    vVI=vv1_fVI;
    
    
    [K1v_fI, K1v_fII, K1v_fIII, K1v_fIV, K1v_fV, K1v_fVI]=eq_moment74(...
        hI,hII,hIII,hIV,hV,hVI,vI,vII,vIII,vIV,vV,vVI);
            
    [K1h_fI, K1h_fII, K1h_fIII, K1h_fIV, K1h_fV, K1h_fVI]=eq_cons74(...
        hI,hII,hIII,hIV,hV,hVI,vI,vII,vIII,vIV,vV,vVI);
    
    hh2_fI=(3/4)*ht_fI+(1/4)*hh1_fI+(1/4)*ddt*K1h_fI;
    hh2_fII=(3/4)*ht_fII+(1/4)*hh1_fII+(1/4)*ddt*K1h_fII;
    hh2_fIII=(3/4)*ht_fIII+(1/4)*hh1_fIII+(1/4)*ddt*K1h_fIII;
    hh2_fIV=(3/4)*ht_fIV+(1/4)*hh1_fIV+(1/4)*ddt*K1h_fIV;
    hh2_fV=(3/4)*ht_fV+(1/4)*hh1_fV+(1/4)*ddt*K1h_fV;
    hh2_fVI=(3/4)*ht_fVI+(1/4)*hh1_fVI+(1/4)*ddt*K1h_fVI;

    vv2_fI=(3/4)*vt_fI+(1/4)*vv1_fI+(1/4)*ddt*K1v_fI;
    vv2_fII=(3/4)*vt_fII+(1/4)*vv1_fII+(1/4)*ddt*K1v_fII;
    vv2_fIII=(3/4)*vt_fIII+(1/4)*vv1_fIII+(1/4)*ddt*K1v_fIII;
    vv2_fIV=(3/4)*vt_fIV+(1/4)*vv1_fIV+(1/4)*ddt*K1v_fIV;
    vv2_fV=(3/4)*vt_fV+(1/4)*vv1_fV+(1/4)*ddt*K1v_fV;
    vv2_fVI=(3/4)*vt_fVI+(1/4)*vv1_fVI+(1/4)*ddt*K1v_fVI;
    
    %% step 3
    hI=hh2_fI;
    hII=hh2_fII;
    hIII=hh2_fIII;
    hIV=hh2_fIV;
    hV=hh2_fV;
    hVI=hh2_fVI;
    
    vI=vv2_fI;
    vII=vv2_fII;
    vIII=vv2_fIII;
    vIV=vv2_fIV;
    vV=vv2_fV;
    vVI=vv2_fVI;
    
    
    [K1v_fI, K1v_fII, K1v_fIII, K1v_fIV, K1v_fV, K1v_fVI]=eq_moment74(...
        hI,hII,hIII,hIV,hV,hVI,vI,vII,vIII,vIV,vV,vVI);
            
    [K1h_fI, K1h_fII, K1h_fIII, K1h_fIV, K1h_fV, K1h_fVI]=eq_cons74(...
        hI,hII,hIII,hIV,hV,hVI,vI,vII,vIII,vIV,vV,vVI);
    
    htnew_fI=(1/3)*ht_fI+(2/3)*hh1_fI+(2/3)*ddt*K1h_fI;
    htnew_fII=(1/3)*ht_fII+(2/3)*hh1_fII+(2/3)*ddt*K1h_fII;
    htnew_fIII=(1/3)*ht_fIII+(2/3)*hh1_fIII+(2/3)*ddt*K1h_fIII;
    htnew_fIV=(1/3)*ht_fIV+(2/3)*hh1_fIV+(2/3)*ddt*K1h_fIV;
    htnew_fV=(1/3)*ht_fV+(2/3)*hh1_fV+(2/3)*ddt*K1h_fV;
    htnew_fVI=(1/3)*ht_fVI+(2/3)*hh1_fVI+(2/3)*ddt*K1h_fVI;
    
    vtnew_fI=(1/3)*vt_fI+(2/3)*vv1_fI+(2/3)*ddt*K1v_fI;
    vtnew_fII=(1/3)*vt_fII+(2/3)*vv1_fII+(2/3)*ddt*K1v_fII;
    vtnew_fIII=(1/3)*vt_fIII+(2/3)*vv1_fIII+(2/3)*ddt*K1v_fIII;
    vtnew_fIV=(1/3)*vt_fIV+(2/3)*vv1_fIV+(2/3)*ddt*K1v_fIV;
    vtnew_fV=(1/3)*vt_fV+(2/3)*vv1_fV+(2/3)*ddt*K1v_fV;
    vtnew_fVI=(1/3)*vt_fVI+(2/3)*vv1_fVI+(2/3)*ddt*K1v_fVI;
    
            
    %% error on h
    t=t+ddt;
    
    [h_fI,v_fI]     = sol_exacte(x_fI  ,y_fI   ,z_fI  ,t);
    [h_fII,v_fII]   = sol_exacte(x_fII ,y_fII  ,z_fII ,t);
    [h_fIII,v_fIII] = sol_exacte(x_fIII,y_fIII ,z_fIII,t);
    [h_fIV,v_fIV]   = sol_exacte(x_fIV ,y_fIV  ,z_fIV ,t);
    [h_fV,v_fV]     = sol_exacte(x_fV  ,y_fV   ,z_fV  ,t);
    [h_fVI,v_fVI]   = sol_exacte(x_fVI ,y_fVI  ,z_fVI ,t);
    
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
        figure(9)
        plot_cs17(n,nn,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,5050,5950,25);
        title(['calculated solution at time = ', num2str(time(end))])
        hold off

        currFrame = getframe;
        writeVideo(vidObj,currFrame);
    end
    
    % snapshot
    if sauvegarde == 1 & mod(iter,floor(Tmax/(3*ddt))) == 0 & strcmp(snapshot,'yes')==1 
        mkdir(['./RK3_results-' jour '/' num2str(ref)])
        close all;
        
        figure(100)
        hFig = figure(100);
        set(gcf,'PaperPositionMode','auto')
        set(hFig, 'Position', [50 50 1000 500])
        plot_cs17(n,nn,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,5050,5950,50);
        title(['calculated solution at time = ', num2str(time(end))])

        print('-dpng', ['./RK3_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_snapshot_intermediaire' num2str(floor(time(end))) '.png'])
        savefig(['./RK3_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_snapshot_intermediaire_' num2str(floor(time(end))) '.fig']);
    end

    %% historique sur la divergence
    [div_fI, div_fII, div_fIII, div_fIV, div_fV, div_fVI]=...
        div74(vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI,n,nn);
    
    Mdivu(iter)=max(max([div_fI div_fII div_fIII div_fIV div_fV div_fVI]));
    
    
end
tend=cputime-tstart;

if strcmp(video,'yes') == 1
    close(vidObj);
    
    fid = fopen('AA_VIDEO_SAVE_RK3.txt','a');
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
    fprintf(fid,'%s\n',['ordre du filtre   : '  opt_ftr1] );
    fprintf(fid,'%s\n',['ordre du filtre   : '  opt_detec] );
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
    fid = fopen('AA_RESULTS_SAVE_RK3.txt','a');
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
    fprintf(fid,'%s\n',['ordre du filtre   : '  opt_ftr1] );
    fprintf(fid,'%s\n',['ordre du filtre   : '  opt_detec] );
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
    mkdir(['./RK3_results-' jour '/' num2str(ref) ])
    print('-dpng', ['./RK3_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_courbe.png'])
    savefig(['./RK3_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_courbe']);
    save(['./RK3_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_erreurdata_test_' num2str(test) '.mat']);
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
    plot_cs7(n,nn,vort_fI,vort_fII,vort_fIII,vort_fIV,vort_fV,vort_fVI)
    title(['vorticity at time : ', num2str(time(end))])
    if sauvegarde == 1
        print('-dpng', ['./RK3_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_snapshot.png'])
        savefig(['./RK3_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_snapshot']);
    end
end

figure(5)
semilogy(time, erri,time, err2, time, err1)
xlabel('time')
ylabel('relative error')
legend('infty norm','norm 2','norm 1')
grid on
if sauvegarde==1
    print('-dpng', ['./RK3_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_erreur.png'])
    savefig(['./RK3_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_erreur']);
end 

figure(6)
plot(time,err_int,time,err_energy,time,err_enstrophy)
legend('mass','energy','potential enstrophy')
xlabel('time')
title('relative quantity')
grid on
if sauvegarde==1
    print('-dpng', ['./RK3_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_conservation.png'])
    savefig(['./RK3_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_conservation']);
end 

figure(7)
semilogy(time,stabi)
xlabel('time')
title('stabilization of height')
grid on

figure(8)
hFig = figure(8);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [50 50 1000 500])
plot_cs7(n,nn,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI);
title(['calculated solution at time = ', num2str(time(end))])
if sauvegarde==1
    print('-dpng', ['./RK3_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_solution.png'])
    savefig(['./RK3_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_solution']);
end


figure(9)
hFig = figure(9);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [50 50 1000 500])
plot_cs17(n,nn,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,5050,5950,50);
title(['calculated solution at time = ', num2str(time(end))])
if sauvegarde==1
    print('-dpng', ['./RK3_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_snapshot_solution.png'])
    savefig(['./RK3_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_snapshot_solution']);
end 

[div_fI, div_fII, div_fIII, div_fIV, div_fV, div_fVI]=...
        div74(vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI,n,nn);
    
figure(10)
hFig = figure(10);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [50 50 1000 500])
plot_cs7(n,nn,div_fI, div_fII, div_fIII, div_fIV, div_fV, div_fVI)
title(['divergence at time : ', num2str(time(end))])

figure(11)
plot(time,Mdivu)

[detec]=filtre74(na,opt_detec);
[det_fI,det_fII,det_fIII,det_fIV,det_fV,det_fVI]=det74(ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI, n, nn);

figure(12)
hFig = figure(12);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [50 50 1000 500])
plot_cs7(n,nn,det_fI,det_fII,det_fIII,det_fIV,det_fV,det_fVI)
title(['detection oscillations : ', num2str(time(end))])
if sauvegarde==1
    print('-dpng', ['./RK3_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_gibbs.png'])
    savefig(['./RK3_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_gibbs']);
end 


fig_placier
disp(['temps de calcul (sans les graphiques) : ', num2str(tend)])