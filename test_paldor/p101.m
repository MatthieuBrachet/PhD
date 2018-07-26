%% ************************************************************************
% Solve LSWEC on the Cubed-Sphere.
% *** options :
% scheme : numerical spatial scheme used. 
% video  : 'yes' ou 'no', do a video or not,
% nper   :  periodicity of frames in the video.
% sauvegarde = 0 (do not save data), 1 (save all data).
% opt_ftr : explicit (redonnet10, redonnet8, ...).
%% ************************************************************************
clc; clear all; close all;timest=clock;
format long

global n nn
global radius dxi
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr scheme nrm
global gp hp omega
global waveFlag

comment='.';
video = 'no';
sauvegarde = 0;
opt_ftr='redonnet10';
scheme='compact4';
nrm='int';

n=63;
mod101
disp('mod101 : ok')

%% ************************************************************************
waveFlag='EIG';
if strcmp(waveFlag,'EIG')==1 || strcmp(waveFlag,'WIG')==1
    lat_deg=36;
    per=3.16/24;
elseif strcmp(waveFlag,'Rossby') == 1
    lat_deg=44;
    per=12.03*24*3600;
end
per_s=per*3600*24;
nday=100*per_s;

%% ************************************************************************
nper=5;
tstart=cputime;
ref=floor(10000*now);
jour=date;

%% ************************************************************************
cfl=.9;
ddt=cfl*radius*dxi/sqrt(gp*hp);
tmax=nday;
itermax=10e15;

tstart=cputime;
ref=floor(10000*now);
jour=date;

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

iter=1;
time(1)=t; erri(1)=0; err_int(1)=1;
%% *** video option *******************************************************
if strcmp(video,'yes')==1
    mkdir(['./RK4_video-' jour ])
    vidObj=VideoWriter(['./RK4_video-' jour '/ref_' num2str(ref)],'Motion JPEG AVI');
    vidObj.Quality = 100;
    open(vidObj);
    set(gca,'nextplot','replacechildren');
end

%% *** iterations *********************************************************
equateur=[];
hovmoller=[];
time_hov=[];
point=[];
while iter<itermax && t+ddt<tmax
    %% filtrage
    [ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI]=ftr_mixte101(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, n, nn);
    [vt_fI(:,:,1), vt_fII(:,:,1), vt_fIII(:,:,1), vt_fIV(:,:,1), vt_fV(:,:,1), vt_fVI(:,:,1)]=ftr_mixte101(vt_fI(:,:,1), vt_fII(:,:,1), vt_fIII(:,:,1), vt_fIV(:,:,1), vt_fV(:,:,1), vt_fVI(:,:,1),n,nn);
    [vt_fI(:,:,2), vt_fII(:,:,2), vt_fIII(:,:,2), vt_fIV(:,:,2), vt_fV(:,:,2), vt_fVI(:,:,2)]=ftr_mixte101(vt_fI(:,:,2), vt_fII(:,:,2), vt_fIII(:,:,2), vt_fIV(:,:,2), vt_fV(:,:,2), vt_fVI(:,:,2),n,nn);
    [vt_fI(:,:,3), vt_fII(:,:,3), vt_fIII(:,:,3), vt_fIV(:,:,3), vt_fV(:,:,3), vt_fVI(:,:,3)]=ftr_mixte101(vt_fI(:,:,3), vt_fII(:,:,3), vt_fIII(:,:,3), vt_fIV(:,:,3), vt_fV(:,:,3), vt_fVI(:,:,3),n,nn);
    
    %% conservation
    [M] = mass( ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI );
    err_int(iter)=M/Mref;
   
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
    
    [K1v_fI, K1v_fII, K1v_fIII, K1v_fIV, K1v_fV, K1v_fVI]=eq_moment101(hh_fI,...
                 hh_fII, hh_fIII, hh_fIV, hh_fV, hh_fVI, vv_fI, vv_fII, vv_fIII, vv_fIV, vv_fV, vv_fVI);
            
    [K1h_fI, K1h_fII, K1h_fIII, K1h_fIV, K1h_fV, K1h_fVI]=eq_cons101(hh_fI,...
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
    
    [K2v_fI, K2v_fII, K2v_fIII, K2v_fIV, K2v_fV, K2v_fVI]=eq_moment101(hh_fI,...
                hh_fII, hh_fIII, hh_fIV, hh_fV, hh_fVI, vv_fI, vv_fII, vv_fIII, vv_fIV, vv_fV, vv_fVI);
            
    [K2h_fI, K2h_fII, K2h_fIII, K2h_fIV, K2h_fV, K2h_fVI]=eq_cons101(hh_fI,...
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
    
    [K3v_fI, K3v_fII, K3v_fIII, K3v_fIV, K3v_fV, K3v_fVI]=eq_moment101(hh_fI,...
                 hh_fII, hh_fIII, hh_fIV, hh_fV, hh_fVI, vv_fI, vv_fII, vv_fIII, vv_fIV, vv_fV, vv_fVI);
            
    [K3h_fI, K3h_fII, K3h_fIII, K3h_fIV, K3h_fV, K3h_fVI]=eq_cons101(hh_fI,...
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
    
    [K4v_fI, K4v_fII, K4v_fIII, K4v_fIV, K4v_fV, K4v_fVI]=eq_moment101(hh_fI,...
                 hh_fII, hh_fIII, hh_fIV, hh_fV, hh_fVI, vv_fI, vv_fII, vv_fIII, vv_fIV, vv_fV, vv_fVI);
            
    [K4h_fI, K4h_fII, K4h_fIII, K4h_fIV, K4h_fV, K4h_fVI]=eq_cons101(hh_fI,...
                 hh_fII, hh_fIII, hh_fIV, hh_fV, hh_fVI, vv_fI, vv_fII, vv_fIII, vv_fIV, vv_fV, vv_fVI);       
            
            
    %% Assemblage data
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
    
    %% coupe equatoriale
    if t>=95*per_s-ddt 
        %% plot u
        [ft_fI,ft_fII,ft_fIII,ft_fIV,ft_fV,ft_fVI] = plot_u(vtnew_fI,vtnew_fII,vtnew_fIII,vtnew_fIV,vtnew_fV,vtnew_fVI);
        funI=ft_fI;
        funII=ft_fII;
        funIII=ft_fIII;
        funIV=ft_fIV;
        funV=ft_fV;
        funVI=ft_fVI;
        
        %% coupes
        [equa_x,coupe] = coupe_eq(funI,funII,funIII,funIV);
        equateur=[equateur; coupe];

        lat=lat_deg/180*pi;
        [lambdai,coupe2,cc] = coupe_lat(funI,funII,funIII,funIV,funV,funVI,lat);
        point=[point cc];
        
        hovmoller=[hovmoller; coupe2];
        time_hov=[time_hov t];
    end

    %% error on h
    t=t+ddt;
    
    [ h_fI,    v_fI]   = sol_exacte(x_fI,   y_fI,   z_fI,   t);
    [ h_fII,   v_fII]  = sol_exacte(x_fII,  y_fII,  z_fII,  t);
    [ h_fIII,  v_fIII] = sol_exacte(x_fIII, y_fIII, z_fIII, t);
    [ h_fIV,   v_fIV]  = sol_exacte(x_fIV,  y_fIV,  z_fIV,  t);
    [ h_fV,    v_fV]   = sol_exacte(x_fV,   y_fV,   z_fV,   t);
    [ h_fVI,   v_fVI]  = sol_exacte(x_fVI,  y_fVI,  z_fVI,  t);
    
    err_fI   = htnew_fI  - h_fI;
    err_fII  = htnew_fII - h_fII;
    err_fIII = htnew_fIII- h_fIII;
    err_fIV  = htnew_fIV - h_fIV;
    err_fV   = htnew_fV  - h_fV ;
    err_fVI  = htnew_fVI - h_fVI;
    
    str='infty';
        [~,~,~,~,~,~,nrmger]=...
      nrm101(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
    [~,~,~,~,~,~,nrmref]=nrm101(h_fI,h_fII,h_fIII,h_fIV,h_fV,h_fVI,n,nn,str);
    erri(iter)=nrmger./nrmref;
    
    str='2';
        [~,~,~,~,~,~,nrmger]=...
      nrm101(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
    [~,~,~,~,~,~,nrmref]=nrm101(h_fI,h_fII,h_fIII,h_fIV,h_fV,h_fVI,n,nn,str);
    err2(iter)=nrmger./nrmref;
    
    str='1';
        [~,~,~,~,~,~,nrmger]=...
      nrm101(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
    [~,~,~,~,~,~,nrmref]=nrm101(h_fI,h_fII,h_fIII,h_fIV,h_fV,h_fVI,n,nn,str);
    err1(iter)=nrmger./nrmref;
       
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

    %% video
    if strcmp(video,'yes')==1 && mod(iter,nper) == 0
        [ft_fI,ft_fII,ft_fIII,ft_fIV,ft_fV,ft_fVI] = plot_u(vt_fI,vt_fII,vt_fIII,vt_fIV,vt_fV,vt_fVI);
        %ft_fI=ht_fI; ft_fII=ht_fII;
        %ft_fIII=ht_fIII; ft_fIV=ht_fIV;
        %ft_fV=ht_fV; ft_fVI=ht_fVI;
        
        
        [nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,fmax]=...
            nrm101(ft_fI,ft_fII,ft_fIII,ft_fIV,ft_fV,ft_fVI,n,nn,'infty');
        
        clf
        hFig = figure(100);
        set(gcf,'PaperPositionMode','auto')
        set(hFig, 'Position', [500 500 1600 800])
        plot_cs100(n,nn,ft_fI./fmax,ft_fII./fmax,ft_fIII./fmax,ft_fIV./fmax,ft_fV./fmax,ft_fVI./fmax);
        title(['solution at time : ', num2str(time(end))])
        colorbar

        hold off
        currFrame = getframe;
        writeVideo(vidObj,currFrame);
        
    end
    
    %% time update
    time(iter)=iter*ddt/3600/24;
    iter=iter+1;
    
    %% affichage intermédiaire
    clc; 
    disp(real([iter min(itermax,floor(tmax/ddt))]));
end

%% résultats
close all;
disp('calcul : OK');
disp('plot...');
timeend=clock;
tend=cputime-tstart;

if strcmp(video,'yes') == 1
    close(vidObj);
    
    fid = fopen('AA_VIDEO_SAVE_RK4.txt','a');
    fprintf(fid,'%s\n',['date : ', jour]);
    fprintf(fid,'%s\n',['ref. : ', num2str(ref)]);
    fprintf(fid,'%s\n','***********************************');
    fprintf(fid,'%s\n','---------- numerical data ---------');
    fprintf(fid,'%s\n',['number of points  : ', num2str(n)] );
    fprintf(fid,'%s\n',['time step         : ', num2str(ddt)] );
    fprintf(fid,'%s\n',['filter order      : '  opt_ftr] );
    fprintf(fid,'%s\n',['scheme            : '  scheme] );
    fprintf(fid,'%s\n',['quadrature        : '  nrm] );
    fprintf(fid,'%s\n','---------- physical data ----------');
    fprintf(fid,'%s\n',['gravity g              : ', num2str(gp)] );
    fprintf(fid,'%s\n',['caracteristic height   : ', num2str(hp)] );
    fprintf(fid,'%s\n',['coriolis parammeter    : ', num2str(omega)] );
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
    fprintf(fid,'%s\n','---------- numerical data ---------');
    fprintf(fid,'%s\n',['number of points  : ', num2str(n)] );
    fprintf(fid,'%s\n',['time step         : ', num2str(ddt)] );
    fprintf(fid,'%s\n',['filter order      : '  opt_ftr] );
    fprintf(fid,'%s\n',['scheme            : '  scheme] );
    fprintf(fid,'%s\n',['quadrature        : '  nrm] );
    fprintf(fid,'%s\n','---------- physical data ----------');
    fprintf(fid,'%s\n',['gravity g              : ', num2str(gp)] );
    fprintf(fid,'%s\n',['caracteristic height   : ', num2str(hp)] );
    fprintf(fid,'%s\n',['coriolis parammter     : ', num2str(omega)] );
    fprintf(fid,'%s\n','------------ comment --------------');
    fprintf(fid,'%s\n', comment );
    fprintf(fid,'%s\n','***********************************');
    fprintf(fid,'%s\n','  ');
    fprintf(fid,'%s\n','  ');
    fclose(fid);
end
 
%% *** PLOT ***************************************************************
figure(1)
plot_cs11(n,nn,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI);
title(['Time : ', num2str(time(end))])
colorbar

[up_fI,up_fII,up_fIII,up_fIV,up_fV,up_fVI] = plot_u(vt_fI,vt_fII,vt_fIII,vt_fIV,vt_fV,vt_fVI);
figure(2)
plot_cs100(n,nn,up_fI,up_fII,up_fIII,up_fIV,up_fV,up_fVI);
title(['Time : ', num2str(time(end))])
colorbar
%axis([-pi pi -pi/4 pi/4])

figure(3)
plot(time,err1,time,err2,time,erri,'Linewidth',2)
legend('norm 1','norm 2','norm \infty')
grid on

[X,Y]=meshgrid(equa_x,time_hov);
figure(4)
surf(X,Y,equateur)
shading interp
title('equatorial section')
ylabel('time')
xlabel('longitude')
colorbar
view(2)

% f = fit(time_hov',point','fourier2');
% Cn=f.w/(2*pi);
% 
% figure(5)
% plot(f,time_hov,point)
% title(['point at : (0,' num2str(lat_deg) ')'])
% grid on



hov_max=max(max(hovmoller));
[X,Y]=meshgrid(lambdai,time_hov);

Cn= fit_hov(X/pi*180,Y,hovmoller./hov_max);

C=getPhaseSpeed2(waveFlag);
ad=1./C;
bd=97.5*per_s;
yd=ad.*lambdai+bd;
adn=Cn;
ydn=adn.*lambdai+bd;



figure(6)
hold on
surf(X,Y,hovmoller./hov_max)
plot3(lambdai,yd,ones(size(yd)),'w-',lambdai,ydn,ones(size(ydn)),'w--','Linewidth',2)
hold off
shading interp
title(['Hovmöller diagram of ' waveFlag ' at \lambda = ' num2str(lat_deg)])
ylabel('time')
xlabel('longitude')
axis([-pi/2 pi/2 95*per_s 100*per_s-ddt])
colorbar
view(2)

abs(C-Cn)*100


fig_placier