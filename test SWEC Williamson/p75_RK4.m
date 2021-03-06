%% ************************************************************************
% Resolution de SWEC sur la Cubed-Sphere.
syms 
% *** options :
% test = 0 : test 2 of Williamson & al.,
%        1 : test 5 of Williamson & al..
% scheme : numerical spatial scheme used. 
% video : 'yes' ou 'no', do a video or not,
% nper  :  periodicity of frames in the video.
% sauvegarde = 0 (do not save data), 1 (save all data).
% opt_ftr : explicit filtering of S. Redonnet (=2, 4, 6, 8, 10, filter 
%           order = 0, no filtering)
% alfa_ftr : parameter for implicit fliter (type visbal only, 
%            if alpha_ftr=0, the filter is equivalent to Redonnet filter,
%            if alpha_ftr=+-0.5, the filter is inexistant).
%
%% ************************************************************************
clc; clear all; close all;
format short
global n nn dxi
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr alfa_ftr test scheme
global gp h0 u0 radius omega
global alpha

test=1;
video = 'yes';
nper=10;
sauvegarde = 0;
opt_ftr='redonnet0';
alfa_ftr=0.45;
scheme='compact4';
snapshot='no';

n=31; % for snapshot, n must be in the form 2^m-1 !
mod74

ccor=radius*omega;
cgrav=sqrt(h0*gp);
cvit=u0;
c=max([cgrav,ccor,cvit]);

cfl=0.9;
ddt=radius*dxi*cfl/c;
ndaymax=5;
Tmax=ndaymax*3600*24;
itermax=5000;
comment='.';

tstart=cputime;
%% *** test data **********************************************************

if test == 0
    alpha=pi/7;
    u0=2*pi*radius/(12*24*3600);
    h0=2.94*10^4/gp;
elseif test == 1
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

if strcmp(snapshot,'yes')==1
    disp('Please wait')
    [vort_I, vort_II, vort_III, vort_IV, vort_V, vort_VI]=vort74(vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI, n, nn);

    figure(1)
    plot_cs7(n,nn, vort_I, vort_II, vort_III, vort_IV, vort_V, vort_VI)

    pause
end

%% *** quantités a conserver **********************************************
[Mref] = mass( ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI );
[Eref] = energy(ht_fI,ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI);
[PEref] = enstrophy(ht_fI,ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI);

iter=0; FTR=0;
time(1)=t; erri(1)=0; err_int(1)=1;
%% *** video option *******************************************************
if strcmp(video,'yes')==1
    nFrames = min(itermax,floor(Tmax/(nper*ddt)));
    mov(1:nFrames) = struct('cdata', [],'colormap', []);
    set(gca,'nextplot','replacechildren');
end

%% *** iterations *********************************************************
ref=floor(10000*now);
jour=date;
while t<Tmax && iter<itermax
    iter=iter+1;
    clc; 
    disp([iter min(itermax,floor(Tmax/ddt)) erri(end) err_int(end)]);
    %% Filtering
    e=1;
    iterf=0;
    if strcmp(opt_ftr,'redonnet0')==0
        [hs_fI] = relief(x_fI,y_fI,z_fI);
        hstar_fI=ht_fI-hs_fI;
        [hs_fII] = relief(x_fII,y_fII,z_fII);
        hstar_fII=ht_fII-hs_fII;
        [hs_fIII] = relief(x_fIII,y_fIII,z_fIII);
        hstar_fIII=ht_fIII-hs_fIII;
        [hs_fIV] = relief(x_fIV,y_fIV,z_fIV);
        hstar_fIV=ht_fIV-hs_fIV;
        [hs_fV] = relief(x_fV,y_fV,z_fV);
        hstar_fV=ht_fV-hs_fV;
        [hs_fVI] = relief(x_fVI,y_fVI,z_fVI);
        hstar_fVI=ht_fVI-hs_fVI;
        
        for p=1:3
            [vt_fI(:,:,p),vt_fII(:,:,p),vt_fIII(:,:,p),vt_fIV(:,:,p),vt_fV(:,:,p),vt_fVI(:,:,p)]=...
                ftr74(vt_fI(:,:,p),vt_fII(:,:,p),vt_fIII(:,:,p),vt_fIV(:,:,p),vt_fV(:,:,p),vt_fVI(:,:,p),n,nn);
        end
        [hstar_fI,hstar_fII,hstar_fIII,hstar_fIV,hstar_fV,hstar_fVI]=...
            ftr74(hstar_fI,hstar_fII,hstar_fIII,hstar_fIV,hstar_fV,hstar_fVI,n,nn);
    
        ht_fI=hstar_fI+hs_fI;
        ht_fII=hstar_fII+hs_fII;
        ht_fIII=hstar_fIII+hs_fIII;
        ht_fIV=hstar_fIV+hs_fIV;
        ht_fV=hstar_fV+hs_fV;
        ht_fVI=hstar_fVI+hs_fVI;
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
        close all
        
        figure(9)
        plot_cs17(n,nn,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,5050,5950,25);
        title(['calculated solution at time = ', num2str(time(end))])
        hold off;
        mov(iter) = getframe(gcf, [0 0 560 420]);
    end
    
    % snapshot
    if sauvegarde == 1 & mod(iter,floor(Tmax/(3*ddt))+1) == 0 & strcmp(snapshot,'yes')==1 
        mkdir(['./RK4_results2-' jour '/' num2str(ref) ])
        close all;
        
        figure(100)
        plot_cs17(n,nn,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,5050,5950,10);
        title(['calculated solution at time = ', num2str(time(end))])

        print('-dpng', ['./RK4_results2-' jour '/ref_' num2str(ref) '_snapshot_intermediaire' num2str(floor(time(end))) '.png'])
        savefig(['./RK4_results2-' jour '/ref_' num2str(ref) '_snapshot_intermediaire_' num2str(floor(time(end))) '.fig']);
    end

    
end
tend=cputime-tstart;

if strcmp(video,'yes') == 1
    mkdir(['./RK4_video2-' jour ])
    movie2avi(mov, ['./RK4_video2-' jour '/ref_' num2str(ref) '.avi'], 'compression', 'None');
    
    fid = fopen('AA_VIDEO_SAVE_RK4.txt','a');
    fprintf(fid,'%s\n',['date : ', jour]);
    fprintf(fid,'%s\n',['ref. : ', num2str(ref)]);
    fprintf(fid,'%s\n','***********************************');
    fprintf(fid,'%s\n',['test : ', num2str(test)]);
    fprintf(fid,'%s\n','---------- numerical data ---------');
    fprintf(fid,'%s\n',['number of points  : ', num2str(n)] );
    fprintf(fid,'%s\n',['time step         : ', num2str(ddt)] );
    fprintf(fid,'%s\n',['cfl               : ', num2str(cfl)] );
    fprintf(fid,'%s\n',['ordre du filtre   : '  opt_ftr] );
    fprintf(fid,'%s\n',['alpha_ftr         : ', num2str(alfa_ftr)] );
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
    fprintf(fid,'%s\n',['ordre du filtre   : '  opt_ftr] );
    fprintf(fid,'%s\n',['alpha_ftr         : ', num2str(alfa_ftr)] );
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
    mkdir(['./RK4_results2-' jour '/' num2str(ref) ])
    print('-dpng', ['./RK4_results2-' jour '/ref_' num2str(ref) '_courbe.png'])
    savefig(['./RK4_results2-' jour '/' num2str(ref) '/ref_' num2str(ref) '_courbe']);
    save(['./RK4_results2-' jour '/' num2str(ref) '/ref_' num2str(ref) '_erreurdata_test_' num2str(test) '.mat']);
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
    plot_cs7(n,nn,vort_fI,vort_fII,vort_fIII,vort_fIV,vort_fV,vort_fVI)
    title(['vorticity at time : ', num2str(time(end))])
    print('-dpng', ['./RK4_results2-' jour '/' num2str(ref) '/ref_' num2str(ref) '_snapshot.png'])
    savefig(['./RK4_results2-' jour '/' num2str(ref) '/ref_' num2str(ref) '_snapshot']);
end

figure(5)
semilogy(time, erri,time, err2, time, err1)
xlabel('time')
ylabel('relative error')
legend('infty norm','norm 2','norm 1')
grid on
if sauvegarde==1
    print('-dpng', ['./RK4_results2-' jour '/' num2str(ref) '/ref_' num2str(ref) '_erreur.png'])
    savefig(['./RK4_results2-' jour '/' num2str(ref) '/ref_' num2str(ref) '_erreur']);
end 

figure(6)
plot(time,err_int,time,err_energy,time,err_enstrophy)
legend('mass','energy','potential enstrophy')
xlabel('time')
title('relative quantity')
grid on
if sauvegarde==1
    print('-dpng', ['./RK4_results2-' jour '/' num2str(ref) '/ref_' num2str(ref) '_conservation.png'])
    savefig(['./RK4_results2-' jour '/' num2str(ref) '/ref_' num2str(ref) '_conservation']);
end 

figure(7)
semilogy(time,stabi)
xlabel('time')
title('stabilization of height')
grid on

figure(8)
plot_cs7(n,nn,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI);
title(['calculated solution at time = ', num2str(time(end))])


figure(9)
plot_cs17(n,nn,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,5050,5950,25);
title(['calculated solution at time = ', num2str(time(end))])
if sauvegarde==1
    print('-dpng', ['./RK4_results2-' jour '/' num2str(ref) '/ref_' num2str(ref) '_snapshot_solution.png'])
    savefig(['./RK4_results2-' jour '/' num2str(ref) '/ref_' num2str(ref) '_snapshot_solution']);
end 

if test == 1
    [hs_fI] = relief(x_fI,y_fI,z_fI);
    [hs_fII] = relief(x_fII,y_fII,z_fII);
    [hs_fIII] = relief(x_fIII,y_fIII,z_fIII);
    [hs_fIV] = relief(x_fIV,y_fIV,z_fIV);
    [hs_fV] = relief(x_fV,y_fV,z_fV);
    [hs_fVI] = relief(x_fVI,y_fVI,z_fVI);
    figure(10)
    plot_cs7(n,nn,hs_fI,hs_fII,hs_fIII,hs_fIV,hs_fV,hs_fVI);
    title('relief on the sphere')
end

fig_placier
disp(['temps de calcul (sans les graphiques) : ', num2str(tend)])