clc; clear all; close all; format short;

global n nn dxi
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr test
global hp gp u0

test=2;
video = 'no';
opt_ftr=10;
n=20;
mod72

cfl=0.5;
ddt=radius*dxi*cfl/u0;
ndaymax=12;
Tmax=ndaymax*3600*24;
itermax=20;

%% *** initialisation des données
[ vt_fI, ht_fI ] = fun_init(x_fI,y_fI,z_fI);
[ vt_fII, ht_fII ] = fun_init(x_fII,y_fII,z_fII);
[ vt_fIII, ht_fIII ] = fun_init(x_fIII,y_fIII,z_fIII);
[ vt_fIV, ht_fIV ] = fun_init(x_fIV,y_fIV,z_fIV);
[ vt_fV, ht_fV ] = fun_init(x_fV,y_fV,z_fV);
[ vt_fVI, ht_fVI ] = fun_init(x_fVI,y_fVI,z_fVI);

hte_fI=ht_fI;
hte_fII=ht_fII;
hte_fIII=ht_fIII;
hte_fIV=ht_fIV;
hte_fV=ht_fV;
hte_fVI=ht_fVI;
%% *** video option *******************************************************
if strcmp(video,'yes')==1
    nFrames = min(itermax,floor(Tmax/ddt));
    mov(1:nFrames) = struct('cdata', [],'colormap', []);
    set(gca,'nextplot','replacechildren');
end

t=0;iter=0;
time(1)=0; err(1)=0;
while t<Tmax & iter<itermax
    t=t+ddt;
    iter=iter+1;
    
    clc; [iter min(itermax,floor(Tmax/ddt)) err(end)]
    %% Filtrage
    for p=1:3
        [vt_fI(:,:,p),vt_fII(:,:,p),vt_fIII(:,:,p),vt_fIV(:,:,p),vt_fV(:,:,p),vt_fVI(:,:,p)]=...
            ftr72(vt_fI(:,:,p),vt_fII(:,:,p),vt_fIII(:,:,p),vt_fIV(:,:,p),vt_fV(:,:,p),vt_fVI(:,:,p),n,nn);
    end
    [ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI]=...
        ftr72(ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,n,nn);

    %% K1
    % premiere equation
    
    % GRADIENT
    [grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
        gr72(ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,n,nn);
    % CALCUL DE CORIOLIS
    [cor_I,cor_II,cor_III,cor_IV,cor_V,cor_VI]=coriolis(vt_fI,vt_fII,vt_fIII,vt_fIV,vt_fV,vt_fVI);
    
    % ASSEMBLAGE
    vt1_fI=vt_fI-ddt*(gp*grad_I+cor_I);
    vt1_fII=vt_fII-ddt*(gp*grad_II+cor_II);
    vt1_fIII=vt_fIII-ddt*(gp*grad_III+cor_III);
    vt1_fIV=vt_fIV-ddt*(gp*grad_IV+cor_IV);
    vt1_fV=vt_fV-ddt*(gp*grad_V+cor_V);
    vt1_fVI=vt_fVI-ddt*(gp*grad_VI+cor_VI);
    
    % seconde equation
    % divergence
    [div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
        div72(vt_fI,vt_fII,vt_fIII,vt_fIV,vt_fV,vt_fVI,n,nn);
    
    % assemblage
    ht1_fI=ht_fI-ddt*hp*div_fI;
    ht1_fII=ht_fII-ddt*hp*div_fII;
    ht1_fIII=ht_fIII-ddt*hp*div_fIII;
    ht1_fIV=ht_fIV-ddt*hp*div_fIV;
    ht1_fV=ht_fV-ddt*hp*div_fV;
    ht1_fVI=ht_fVI-ddt*hp*div_fVI;
    
    %% K2
    % premiere equation
    
    % GRADIENT
    [grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
        gr72(ht1_fI,ht1_fII,ht1_fIII,ht1_fIV,ht1_fV,ht1_fVI,n,nn);
    % CALCUL DE CORIOLIS
    [cor_I,cor_II,cor_III,cor_IV,cor_V,cor_VI]=coriolis(vt1_fI,vt1_fII,vt1_fIII,vt1_fIV,vt1_fV,vt1_fVI);
    
    % ASSEMBLAGE
    vt2_fI=3/4*vt_fI+1/4*vt1_fI-1/4*ddt*(gp*grad_I+cor_I);
    vt2_fII=3/4*vt_fII+1/4*vt1_fII-1/4*ddt*(gp*grad_II+cor_II);
    vt2_fIII=3/4*vt_fIII+1/4*vt1_fIII-1/4*ddt*(gp*grad_III+cor_III);
    vt2_fIV=3/4*vt_fIV+1/4*vt1_fIV-1/4*ddt*(gp*grad_IV+cor_IV);
    vt2_fV=3/4*vt_fV+1/4*vt1_fV-1/4*ddt*(gp*grad_V+cor_V);
    vt2_fVI=3/4*vt_fVI+1/4*vt1_fVI-1/4*ddt*(gp*grad_VI+cor_VI);

    
    % seconde equation
    % divergence
    [div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
        div72(vt1_fI,vt1_fII,vt1_fIII,vt1_fIV,vt1_fV,vt1_fVI,n,nn);
    
    % assemblage
    ht2_fI=3/4*ht_fI+1/4*ht1_fI-ddt/4*hp*div_fI;
    ht2_fII=3/4*ht_fII+1/4*ht1_fII-ddt/4*hp*div_fII;
    ht2_fIII=3/4*ht_fIII+1/4*ht1_fIII-ddt/4*hp*div_fIII;
    ht2_fIV=3/4*ht_fIV+1/4*ht1_fIV-ddt/4*hp*div_fIV;
    ht2_fV=3/4*ht_fV+1/4*ht1_fV-ddt/4*hp*div_fV;
    ht2_fVI=3/4*ht_fVI+1/4*ht1_fVI-ddt/4*hp*div_fVI;

    %% K2
    % premiere equation
    
    % GRADIENT
    [grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
        gr72(ht2_fI,ht2_fII,ht2_fIII,ht2_fIV,ht2_fV,ht2_fVI,n,nn);
    % CALCUL DE CORIOLIS
    [cor_I,cor_II,cor_III,cor_IV,cor_V,cor_VI]=coriolis(vt2_fI,vt2_fII,vt2_fIII,vt2_fIV,vt2_fV,vt2_fVI);
    
    % ASSEMBLAGE
    vtnew_fI=1/3*vt_fI+2/3*vt2_fI-2/3*ddt*(gp*grad_I+cor_I);
    vtnew_fII=1/3*vt_fII+2/3*vt2_fII-2/3*ddt*(gp*grad_II+cor_II);
    vtnew_fIII=1/3*vt_fIII+2/3*vt2_fIII-2/3*ddt*(gp*grad_III+cor_III);
    vtnew_fIV=1/3*vt_fIV+2/3*vt2_fIV-2/3*ddt*(gp*grad_IV+cor_IV);
    vtnew_fV=1/3*vt_fV+2/3*vt2_fV-2/3*ddt*(gp*grad_V+cor_V);
    vtnew_fVI=1/3*vt_fVI+2/3*vt2_fVI-2/3*ddt*(gp*grad_VI+cor_VI);


    % seconde equation
    % divergence
    [div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
        div72(vt2_fI,vt2_fII,vt2_fIII,vt2_fIV,vt2_fV,vt2_fVI,n,nn);
    
    % assemblage
    htnew_fI=1/3*ht_fI+2/3*ht2_fI-ddt*hp*div_fI;
    htnew_fII=1/3*ht_fII+2/3*ht2_fII-ddt*hp*div_fII;
    htnew_fIII=1/3*ht_fIII+2/3*ht2_fIII-ddt*hp*div_fIII;
    htnew_fIV=1/3*ht_fIV+2/3*ht2_fIV-ddt*hp*div_fIV;
    htnew_fV=1/3*ht_fV+2/3*ht2_fV-ddt*hp*div_fV;
    htnew_fVI=1/3*ht_fVI+2/3*ht2_fVI-ddt*hp*div_fVI;

    %% calcul de l'erreur
    
    err_fI(iter)=max(max(abs(htnew_fI-hte_fI)./abs(hte_fI)));
    err_fII(iter)=max(max(abs(htnew_fII-hte_fII)./abs(hte_fII)));
    err_fIII(iter)=max(max(abs(htnew_fIII-hte_fIII)./abs(hte_fIII)));
    err_fIV(iter)=max(max(abs(htnew_fIV-hte_fIV)./abs(hte_fIV)));
    err_fV(iter)=max(max(abs(htnew_fV-hte_fV)./abs(hte_fV)));
    err_fVI(iter)=max(max(abs(htnew_fVI-hte_fVI)./abs(hte_fVI)));
    
    err(iter)=max([err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI]);
    time(iter)=t/(3600*24);
    
    
    %% mise a jour
    vt_fI=vtnew_fI;
    vt_fII=vtnew_fII;
    vt_fIII=vtnew_fIII;
    vt_fIV=vtnew_fIV;
    vt_fV=vtnew_fV;
    vt_fVI=vtnew_fVI;
    
    ht_fI=htnew_fI;
    ht_fII=htnew_fII;
    ht_fIII=htnew_fIII;
    ht_fIV=htnew_fIV;
    ht_fV=htnew_fV;
    ht_fVI=htnew_fVI;

    %% film
    if strcmp(video,'yes')==1
        figure(100)
        plot_cs14(n,nn,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,6000,11000)
        hold off;
        mov(iter) = getframe(gcf);
    end
    

    
end
ref=floor(10000*now);
if strcmp(video,'yes') == 1
    mkdir(['./video-' date ])
    movie2avi(mov, ['./video-' date '/ref_' num2str(ref) '.avi'], 'compression', 'None');
    
    fid = fopen('AA_VIDEO_SAVE.txt','a');
    fprintf(fid,'%s\n',['date : ', date]);
    fprintf(fid,'%s\n',['ref. : ', num2str(ref)]);
    fprintf(fid,'%s\n','***********************************');
    fprintf(fid,'%s\n',['test : ', num2str(test)]);
    fprintf(fid,'%s\n','---------- numerical data ---------');
    fprintf(fid,'%s\n',['number of points  : ', num2str(n)] );
    fprintf(fid,'%s\n',['time step         : ', num2str(ddt)] );
    fprintf(fid,'%s\n',['cfl               : ', num2str(cfl)] );
    fprintf(fid,'%s\n',['ordre du filtre   : ', num2str(opt_ftr)] );
    fprintf(fid,'%s\n','-------- mathematical data --------');
    fprintf(fid,'%s\n',['ndaymax           : ', num2str(ndaymax)] );
    fprintf(fid,'%s\n','***********************************');
    fprintf(fid,'%s\n','  ');
    fprintf(fid,'%s\n','  ');
    fclose(fid);
end
figure(1)
plot_cs11(n,nn,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI);

figure(2)
plot(time,err)
xlabel('time')
ylabel('erreur relative')
