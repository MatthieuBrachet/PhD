%% ************************************************************************
% Solve a Conservation Law on the Cubed-Sphere grid.
% *** options :
% test = 0  : test 3 in M. Ben-Artzi, J. Falcovitz and P. G. Lefloch
%             article ('Steady state solution in spherical cap.').
%        1  : test 1 in M. Ben-Artzi, J. Falcovitz and P. G. Lefloch
%             article ('Equatorial periodic solution').
%        2  : test 4 in M. Ben-Artzi, J. Falcovitz and P. G. Lefloch
%             article ('Confined solutions').
%        3  : test 2 in M. Ben-Artzi, J. Falcovitz and P. G. Lefloch
%             article ('Steady state solution').
% scheme : numerical spatial scheme used. 
% video  : 'yes' ou 'no', do a video or not,
% nper   :  periodicity of frames in the video.
% sauvegarde = 0 (do not save data), 1 (save all data).
% opt_ftr : explicit (redonnet10, redonnet8, ...).
%% ************************************************************************
clc; clear all; close all;format shorte; timest=clock;
format long

global n nn
global radius dxi
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr test scheme nrm detec

comment='.';
test=1;
video = 'yes';
sauvegarde =1;
filtre='symetric';
opt_ftr='redonnet2';
scheme='compact4';
snapshot='no';
nrm='int';

n=31 ; % for snapshot and better spherical integration (B. Portenelle works), n must be odd !
ndaymax=10/(2*pi);
mod101
ddt=.005;%.96/pi*dxi;
disp('mod101 : ok')

%% ************************************************************************
nper=1;
tstart=cputime;
ref=floor(10000*now);
jour=date;
%% ************************************************************************
Tmax=ndaymax;
itermax=200000;

tstart=cputime;
ref=floor(10000*now);
jour=date;

%% *** initial data *******************************************************
t=0;
[ ht_fI   ] = sol_exacte(x_fI,   y_fI,   z_fI,   t);
[ ht_fII  ] = sol_exacte(x_fII,  y_fII,  z_fII,  t);
[ ht_fIII ] = sol_exacte(x_fIII, y_fIII, z_fIII, t);
[ ht_fIV  ] = sol_exacte(x_fIV,  y_fIV,  z_fIV,  t);
[ ht_fV   ] = sol_exacte(x_fV,   y_fV,   z_fV,   t);
[ ht_fVI  ] = sol_exacte(x_fVI,  y_fVI,  z_fVI,  t);

h_fI=ht_fI;
h_fII=ht_fII;
h_fIII=ht_fIII;
h_fIV=ht_fIV;
h_fV=ht_fV; 
h_fVI=ht_fVI;

%% *** quantités a conserver **********************************************
[Mref] = mass( ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI );

iter=1;
time(1)=t; erri(1)=0; err_int(1)=1;
%% *** video option *******************************************************
if strcmp(video,'yes')==1
    mkdir(['./RK4_video-' jour ])
    vidObj=VideoWriter(['./RK4_video-' jour '/ref_' num2str(ref)],'Motion JPEG AVI');
    
    open(vidObj);
    set(gca,'nextplot','replacechildren');
end

%% *** iterations *********************************************************
while t<Tmax && iter<itermax
    %% filtrage
    if strcmp(filtre,'adaptatif') == 1
        [detec]=filtre101(na,'redonnet10');
        [det_fI,det_fII,det_fIII,det_fIV,det_fV,det_fVI]=det101(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, n, nn);
        
        %[ftr]=filtre101(na,opt_ftr);
        [htff_fI, htff_fII, htff_fIII, htff_fIV, htff_fV, htff_fVI]=ftr72(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, n, nn);
        ht_fI=det_fI.*htff_fI+(1-det_fI).*ht_fI;
        ht_fII=det_fII.*htff_fII+(1-det_fII).*ht_fII;
        ht_fIII=det_fIII.*htff_fIII+(1-det_fIII).*ht_fIII;
        ht_fIV=det_fIV.*htff_fIV+(1-det_fIV).*ht_fIV;
        ht_fV=det_fV.*htff_fV+(1-det_fV).*ht_fV;
        ht_fVI=det_fVI.*htff_fVI+(1-det_fVI).*ht_fVI;
        
    elseif strcmp(filtre,'nonsymetric') == 1
        [ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI]=ftr72(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, n, nn);
        
    elseif strcmp(filtre,'symetric') == 1
        [ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI]=ftr_mixte101(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, n, nn);
        
    elseif strcmp(filtre,'inf')==1
        
    else
        error('Option ''filtre'' is uncorrect.');
    end
    
    %% conservation
    [M] = mass( ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI );
    err_int(iter)=M;

    %% K1
    hh_fI   = ht_fI;
    hh_fII  = ht_fII;
    hh_fIII = ht_fIII;
    hh_fIV  = ht_fIV;
    hh_fV   = ht_fV;
    hh_fVI  = ht_fVI;

    [K1h_fI, K1h_fII, K1h_fIII, K1h_fIV, K1h_fV, K1h_fVI]=eq_cons101(hh_fI,...
                hh_fII, hh_fIII, hh_fIV, hh_fV, hh_fVI);
            
    %% K2
    hh_fI   = ht_fI   + 0.5*ddt*K1h_fI;
    hh_fII  = ht_fII  + 0.5*ddt*K1h_fII;
    hh_fIII = ht_fIII + 0.5*ddt*K1h_fIII;
    hh_fIV  = ht_fIV  + 0.5*ddt*K1h_fIV;
    hh_fV   = ht_fV   + 0.5*ddt*K1h_fV;
    hh_fVI  = ht_fVI  + 0.5*ddt*K1h_fVI;
     
    [K2h_fI, K2h_fII, K2h_fIII, K2h_fIV, K2h_fV, K2h_fVI]=eq_cons101(hh_fI,...
                 hh_fII, hh_fIII, hh_fIV, hh_fV, hh_fVI);
     
    %% K3
    hh_fI   = ht_fI   + 0.5*ddt*K2h_fI;
    hh_fII  = ht_fII  + 0.5*ddt*K2h_fII;
    hh_fIII = ht_fIII + 0.5*ddt*K2h_fIII;
    hh_fIV  = ht_fIV  + 0.5*ddt*K2h_fIV;
    hh_fV   = ht_fV   + 0.5*ddt*K2h_fV;
    hh_fVI  = ht_fVI  + 0.5*ddt*K2h_fVI;
  
    [K3h_fI, K3h_fII, K3h_fIII, K3h_fIV, K3h_fV, K3h_fVI]=eq_cons101(hh_fI,...
                hh_fII, hh_fIII, hh_fIV, hh_fV, hh_fVI);
            
     %% K4
    hh_fI   = ht_fI   + ddt*K3h_fI;
    hh_fII  = ht_fII  + ddt*K3h_fII;
    hh_fIII = ht_fIII + ddt*K3h_fIII;
    hh_fIV  = ht_fIV  + ddt*K3h_fIV;
    hh_fV   = ht_fV   + ddt*K3h_fV;
    hh_fVI  = ht_fVI  + ddt*K3h_fVI;
     
    [K4h_fI, K4h_fII, K4h_fIII, K4h_fIV, K4h_fV, K4h_fVI]=eq_cons101(hh_fI,...
                 hh_fII, hh_fIII, hh_fIV, hh_fV, hh_fVI);       
            
            
    %% Assemblage data
    Sh_fI   = (K1h_fI   + 2*K2h_fI   + 2*K3h_fI   + K4h_fI);
    Sh_fII  = (K1h_fII  + 2*K2h_fII  + 2*K3h_fII  + K4h_fII);
    Sh_fIII = (K1h_fIII + 2*K2h_fIII + 2*K3h_fIII + K4h_fIII);
    Sh_fIV  = (K1h_fIV  + 2*K2h_fIV  + 2*K3h_fIV  + K4h_fIV);
    Sh_fV   = (K1h_fV   + 2*K2h_fV   + 2*K3h_fV   + K4h_fV);
    Sh_fVI  = (K1h_fVI  + 2*K2h_fVI  + 2*K3h_fVI  + K4h_fVI);
    
    htnew_fI    = ht_fI   + ddt/6 * Sh_fI;
    htnew_fII   = ht_fII  + ddt/6 * Sh_fII;
    htnew_fIII  = ht_fIII + ddt/6 * Sh_fIII;
    htnew_fIV   = ht_fIV  + ddt/6 * Sh_fIV;
    htnew_fV    = ht_fV   + ddt/6 * Sh_fV;
    htnew_fVI   = ht_fVI  + ddt/6 * Sh_fVI;

    %% error on h
    t=t+ddt;
    
    err_fI   = htnew_fI  - h_fI;
    err_fII  = htnew_fII - h_fII;
    err_fIII = htnew_fIII- h_fIII;
    err_fIV  = htnew_fIV - h_fIV;
    err_fV   = htnew_fV  - h_fV ;
    err_fVI  = htnew_fVI - h_fVI;
    
    [Er1,Er2,Eri] = err101(err_fI, err_fII, err_fIII,err_fIV,err_fV,err_fVI);
    
    str='infty';
    [nrmref,~,~,~,~,~,~]=nrm101(h_fI,h_fII,h_fIII,h_fIV,h_fV,h_fVI,n,nn,str);
    erri(iter)=Eri./nrmref;
    
    str='2';
    [nrmref,~,~,~,~,~,~]=nrm101(h_fI,h_fII,h_fIII,h_fIV,h_fV,h_fVI,n,nn,str);
    err2(iter)=Er2./nrmref;
    
    str='1';
    [nrmref,~,~,~,~,~,~]=nrm101(h_fI,h_fII,h_fIII,h_fIV,h_fV,h_fVI,n,nn,str);
    err1(iter)=Er1./nrmref;
  
    
    maxi(iter)=max(max([ht_fI ht_fII ht_fIII ht_fIV ht_fV ht_fVI]));
    mini(iter)=min(min([ht_fI ht_fII ht_fIII ht_fIV ht_fV ht_fVI]));
       
    %% update height
    ht_fI   = htnew_fI;
    ht_fII  = htnew_fII;
    ht_fIII = htnew_fIII;
    ht_fIV  = htnew_fIV;
    ht_fV   = htnew_fV;
    ht_fVI  = htnew_fVI;

    %% video
    if strcmp(video,'yes')==1 && mod(iter,nper) == 0
        
        clf
        hFig=figure(100);
        set(gcf,'PaperPositionMode','auto')
        set(hFig, 'Position', [50 50 1000 500])
        
        subplot(121)
        plot_cs11(n,nn,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI);
        title(['Temps = ' num2str(t)])
        caxis([-1 1])
        view([-1 0 0])
        
        [ lambdae, hte ] = equateur(ht_fI, ht_fII, ht_fIII, ht_fIV,ht_fV, ht_fVI);
        [x,y] = burgers( 4*(n+1)-1, ddt, t );
        
        subplot(122)
        plot(lambdae,hte,'-o',x,y,'Linewidth',2)
        xticks([0 pi/2 pi 3*pi/2 2*pi])
        xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
        axis([0 2*pi -1.2 1.2])
        grid minor
        legend('Coupe équatoriale','Burgers 1D')
        

        hold off
        currFrame = getframe;
        writeVideo(vidObj,currFrame);
        
        
    end
    
    %% snapshot
    if strcmp(snapshot,'yes') == 1 && mod(iter,floor(Tmax/(2*ddt))) == 0
        mkdir(['./RK4_results-' jour '/' num2str(ref)])
        clf;
 
        figure(100)
        hFig = figure(100);
        set(gcf,'PaperPositionMode','auto')
        set(hFig, 'Position', [50 50 1000 500])
        if test == 0
            plot_cs100(n,nn,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI);
            title(['solution at time : ', num2str(time(end))])
            colorbar
        end
            print('-dpng', ['./RK4_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_snapshot_intermediaire' num2str(floor(100*time(end))) '.png'])
            savefig(['./RK4_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_snapshot_intermediaire_' num2str(floor(100*time(end))) '.fig']);
    end


    
    %% time update
    time(iter)=iter*ddt;
    iter=iter+1;
    
    %% affichage intermédiaire
    clc; 
    disp(real([test iter min(itermax,floor(Tmax/ddt)) erri(end)]));
end


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
    fprintf(fid,'%s\n',['test : ', num2str(test)]);
    fprintf(fid,'%s\n','---------- numerical data ---------');
    fprintf(fid,'%s\n',['number of points  : ', num2str(n)] );
    fprintf(fid,'%s\n',['time step         : ', num2str(ddt)] );
    fprintf(fid,'%s\n',['filtering proc.   : ',  filtre] );
    fprintf(fid,'%s\n',['filter order      : '  opt_ftr] );
    fprintf(fid,'%s\n',['scheme            : '  scheme] );
    fprintf(fid,'%s\n',['quadrature        : '  nrm] );
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
    fprintf(fid,'%s\n',['filtering proc.   : ',  filtre] );
    fprintf(fid,'%s\n',['filter order      : '  opt_ftr] );
    fprintf(fid,'%s\n',['scheme            : '  scheme] );
    fprintf(fid,'%s\n',['quadrature        : '  nrm] );
    fprintf(fid,'%s\n','------------ comment --------------');
    fprintf(fid,'%s\n', comment );
    fprintf(fid,'%s\n','***********************************');
    fprintf(fid,'%s\n','  ');
    fprintf(fid,'%s\n','  ');
    fclose(fid);
end
 
%% *** PLOT ***************************************************************
hmax=max(max(abs([h_fI,h_fII,h_fIII,h_fIV,h_fV,h_fVI])));
hFig = figure(1);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [50 50 1000 500])
plot_cs102(n,nn,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI);
title(['calculated solution at time = ', num2str(time(end))])
colorbar
if sauvegarde==1
    mkdir(['./RK4_results-' jour '/' num2str(ref) ])
    print('-dpng', ['./RK4_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_courbe.png'])
    savefig(['./RK4_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_courbe']);
    save(['./RK4_results-' jour '/' num2str(ref) '/ref' num2str(ref) '_test' num2str(test) '.mat']);
end 

hmax=max(max(abs([h_fI,h_fII,h_fIII,h_fIV,h_fV,h_fVI])));hmax=1;
hFig = figure(2);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [50 50 1000 500])
plot_cs102(n,nn,err_fI./hmax, err_fII./hmax, err_fIII./hmax, err_fIV./hmax, err_fV./hmax, err_fVI./hmax);
title(['Relative error at time = ', num2str(time(end))])
colorbar
if sauvegarde == 1
    print('-dpng', ['./RK4_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_snapshot_err.png'])
    savefig(['./RK4_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_snapshot_err']);
end

figure(3)
plot(time,err1,'-k',time,err2,'--k',time,erri,'.k')
legend('norme 1','norme 2','norme \infty')
grid on
if sauvegarde == 1
    print('-dpng', ['./RK4_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_err.png'])
    savefig(['./RK4_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_err']);
end

[ lambdae, hte ] = equateur(ht_fI, ht_fII, ht_fIII, ht_fIV,ht_fV, ht_fVI);
[x,y] = burgers( 4*(n+1)-1, ddt, t );
figure(4)
plot(lambdae,hte,'-o',x,y,'Linewidth',2)
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
axis([0 2*pi -.4 .4])
grid minor
legend('Coupe équatoriale','Burgers 1D')
if sauvegarde == 1
    print('-dpng', ['./RK4_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_equateur.png'])
    savefig(['./RK4_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_equateur']);
end

figure(5)
plot(time,err_int,'b-','Linewidth',2)
grid minor
xlabel('Temps')
ylabel('Erreur de conservation')
if sauvegarde == 1
    print('-dpng', ['./RK4_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_conservation.png'])
    savefig(['./RK4_results-' jour '/' num2str(ref) '/ref_' num2str(ref) '_conservation']);
end

figure(6)
plot(time,maxi,time,mini)
legend('maximum','minimum')
grid on

fig_placier
clc
disp(['temps de calcul (sans les graphiques) : ', num2str(tend)])
err1(end)
err2(end)
erri(end)

timest
timeend

%% save(['sol_ref_test_' num2str(test) '_N_' num2str(n) '.mat'])