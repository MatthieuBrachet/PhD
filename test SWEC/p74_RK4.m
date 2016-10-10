clc; clear all; close all;
format short

global n nn dxi
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr test scheme
global gp h0 u0 radius omega
global teta0 teta1

%% ************************************************************************
% Resolution de SWEC sur la Cubed-Sphere.
%
% *** options :
% test = 0 : test de J. Galewsky stationnaire;
%      = 1 : test de J. Galewsky avec perturbation.
% scheme : numerical spatial scheme used. 
% video : 'yes' ou 'no', faire une video ou non.
% sauvegarde = 0 (ne rien sauvegarder), 1 (sauvegarder toutes les valeurs
%          finales).
% opt_ftr : filtre explicite de S. Redonnet (=2, 4, 6, 8, 10, ordre du
%          filtre; =0, pas de filtrage)
%
%% ************************************************************************

test=1;
video = 'no';
sauvegarde = 1;
opt_ftr=10;
scheme='compact4';
snapshot='no';

n=31;
teta0=pi/7;
teta1=pi/2-teta0;
mod74

ccor=radius*omega;
cgrav=sqrt(h0*gp);
cvit=u0;
c=max([cgrav,ccor,cvit]);

cfl=0.5;
ddt=radius*dxi*cfl/c;
ndaymax=1;
Tmax=ndaymax*3600*24;
itermax=5000;

comment='Start Galewsky benchmark.';

tstart=cputime;
%% *** initialisation des données *****************************************
t=0;
[ ht_fI,    vt_fI]   = sol_exacte(x_fI,   y_fI,   z_fI,   t);
[ ht_fII,   vt_fII]  = sol_exacte(x_fII,  y_fII,  z_fII,  t);
[ ht_fIII,  vt_fIII] = sol_exacte(x_fIII, y_fIII, z_fIII, t);
[ ht_fIV,   vt_fIV]  = sol_exacte(x_fIV,  y_fIV,  z_fIV,  t);
[ ht_fV,    vt_fV]   = sol_exacte(x_fV,   y_fV,   z_fV,   t);
[ ht_fVI,   vt_fVI]  = sol_exacte(x_fVI,  y_fVI,  z_fVI,  t);

if strcmp(snapshot,'yes')==1
%    disp('vort... please wait')
%     [vort_I, vort_II, vort_III, vort_IV, vort_V, vort_VI]=vort74(vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI, n, nn);
% 
%     figure(1)
%     plot_cs7(n,nn, vort_I, vort_II, vort_III, vort_IV, vort_V, vort_VI)
    [div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=div74(vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI,n,nn);
          
    figure(100)
    plot_cs7(n,nn,div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI)
    pause
end

%% *** quantités a conserver **********************************************
[~,~,~,~,~,~,intref]=nrm74(ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,n,nn,'int');
[~,~,~,~,~,~,nrmref]=nrm74(ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,n,nn,'infty');

iter=0; FTR=0;
time(1)=t; erri(1)=0;
%% *** video option *******************************************************
if strcmp(video,'yes')==1
    nFrames = min(itermax,floor(Tmax/ddt));
    mov(1:nFrames) = struct('cdata', [],'colormap', []);
    set(gca,'nextplot','replacechildren');
end

%% *** iterations *********************************************************
while t<Tmax & iter<itermax
    iter=iter+1;
    clc; 
    disp([iter min(itermax,floor(Tmax/ddt)) erri(end)]);
    %% Filtrage
    e=1;
    iterf=0;
    if opt_ftr==0
        iterfmax=-1;
    else
        iterfmax=1;
    end
    while e>1.e-4 & iterf<iterfmax
        iterf=iterf+1;
        for ppp=1:3
            [vt_fI(:,:,ppp),vt_fII(:,:,ppp),vt_fIII(:,:,ppp),vt_fIV(:,:,ppp),vt_fV(:,:,ppp),vt_fVI(:,:,ppp)]=...
                ftr74(vt_fI(:,:,ppp),vt_fII(:,:,ppp),vt_fIII(:,:,ppp),vt_fIV(:,:,ppp),vt_fV(:,:,ppp),vt_fVI(:,:,ppp),n,nn);
        end
        [htf_fI,htf_fII,htf_fIII,htf_fIV,htf_fV,htf_fVI]=...
            ftr74(ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,n,nn);
        
        eh_fI   = max(max(abs(htf_fI   - ht_fI   )./abs(ht_fI   )));
        eh_fII  = max(max(abs(htf_fII  - ht_fII  )./abs(ht_fII  )));
        eh_fIII = max(max(abs(htf_fIII - ht_fIII )./abs(ht_fIII )));
        eh_fIV  = max(max(abs(htf_fIV  - ht_fIV  )./abs(ht_fIV  )));
        eh_fV   = max(max(abs(htf_fV   - ht_fV   )./abs(ht_fV   )));
        eh_fVI  = max(max(abs(htf_fVI  - ht_fVI  )./abs(ht_fVI  )));
        e=max([eh_fI,eh_fII,eh_fIII,eh_fIV,eh_fV,eh_fVI]);
        
        ht_fI=htf_fI;     ht_fII=htf_fII;    ht_fIII=htf_fIII;
        ht_fIV=htf_fIV;   ht_fV=htf_fV;      ht_fVI=htf_fVI;
    end
    FTR=max(iterf,FTR);

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


    %% calcul de l'erreur sur h
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
    
    errv_fI=vtnew_fI-v_fI;
    errv_fII=vtnew_fII-v_fII;
    errv_fIII=vtnew_fIII-v_fIII;
    errv_fIV=vtnew_fIV-v_fIV;
    errv_fV=vtnew_fV-v_fV ;
    errv_fVI=vtnew_fVI-v_fVI;
    
    str='infty';
        [~,~,~,~,~,~,nrmger]=...
      nrm74(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
    erri(iter)=nrmger./nrmref;
    
    %% stabilisation
    stab_fI=htnew_fI-ht_fI;
    stab_fII=htnew_fII-ht_fII;
    stab_fIII=htnew_fIII-ht_fIII;
    stab_fIV=htnew_fIV-ht_fIV;
    stab_fV=htnew_fV-ht_fV ;
    stab_fVI=htnew_fVI-ht_fVI;

    [~,~,~,~,~,~,nrmger]=...
      nrm74(stab_fI,stab_fII,stab_fIII,stab_fIV,stab_fV,stab_fVI,n,nn,str);
    stabi(iter)=nrmger;
    
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
    str='int';
    [~,~,~,~,~,~,int]=...
    nrm74(ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,n,nn,str);
    err_int(iter)=int./intref;
    
    %% film
    if strcmp(video,'yes')==1
%         [vort_fI,vort_fII,vort_fIII,vort_fIV,vort_fV,vort_fVI]=...
%             vort74(vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI,n,nn);
%         
%         
%         figure(100)
%         plot_cs7(n,nn,vort_fI,vort_fII,vort_fIII,vort_fIV,vort_fV,vort_fVI)
%         title(['numerical solution at ', num2str(time(end)), 'days'])
%         hold off;
%         mov(iter) = getframe(gcf);

         [div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
            div74(vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI,n,nn);
        
        
        figure(100)
        plot_cs7(n,nn,div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI)
        title(['numerical solution at ', num2str(time(end)*24), 'hours.'])
        hold off;
        mov(iter) = getframe(gcf);


    end

end
tend=cputime-tstart;

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
    fprintf(fid,'%s\n',['ordre du filtre   : ', num2str(opt_ftr)] );
    fprintf(fid,'%s\n',['iterfmax          : ', num2str(iterfmax)] );
    fprintf(fid,'%s\n','---------- physical data ----------');
    fprintf(fid,'%s\n',['gravity g              : ', num2str(gp)] );
    fprintf(fid,'%s\n',['caracteristiv velocity : ', num2str(u0)] );
    fprintf(fid,'%s\n',['coriolis parameter     : ', num2str(omega)] );
    fprintf(fid,'%s\n','------------ comment --------------');
    fprintf(fid,'%s\n',[comment] );
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
    fprintf(fid,'%s\n',['ordre du filtre   : ', num2str(opt_ftr)] );
    fprintf(fid,'%s\n',['iterfmax          : ', num2str(iterfmax)] );
    fprintf(fid,'%s\n','---------- physical data ----------');
    fprintf(fid,'%s\n',['gravity g              : ', num2str(gp)] );
    fprintf(fid,'%s\n',['caracteristiv velocity : ', num2str(u0)] );
    fprintf(fid,'%s\n',['coriolis parameter     : ', num2str(omega)] );
    fprintf(fid,'%s\n','------------ comment --------------');
    fprintf(fid,'%s\n',[comment] );
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
[vort_fI,vort_fII,vort_fIII,vort_fIV,vort_fV,vort_fVI]=...
     vort74(vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI,n,nn);

plot_cs7(n,nn,vort_fI,vort_fII,vort_fIII,vort_fIV,vort_fV,vort_fVI)
title(['vorticity at time : ', num2str(time(end))])

if strcmp(snapshot,'yes')==1
    figure(3)
    plot_cs13(n,nn,vort_fI,vort_fII,vort_fIII,vort_fIV,vort_fV,vort_fVI)
    title(['vorticity at time : ', num2str(time(end))])
    print('-dpng', ['./RK4_results-' date '/ref_' num2str(ref) '_snapshot.png'])
    savefig(['./RK4_results-' date '/ref_' num2str(ref) '_snapshot']);
end

figure(4)
semilogy(time, erri)
xlabel('time')
ylabel('relative error')
legend('infty norm')
grid on
if sauvegarde==1
    print('-dpng', ['./RK4_results-' date '/ref_' num2str(ref) '_erreur.png'])
    savefig(['./RK4_results-' date '/ref_' num2str(ref) '_erreur']);
end 

figure(5)
semilogy(time,err_int)
legend('mass')
xlabel('time')
title('relative quantity')
grid on
if sauvegarde==1
    mkdir(['./RK4_results-' date ]);
    print('-dpng', ['./RK4_results-' date '/ref_' num2str(ref) '_conservation.png'])
    savefig(['./RK4_results-' date '/ref_' num2str(ref) '_conservation']);
end 

disp(['temps de calcul (sans les graphiques) : ', num2str(tend)])
