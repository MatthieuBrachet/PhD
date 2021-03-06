clc; clear all; close all; format short;

global n nn dxi
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr test
global hp gp u0 radius omega

test=1;
video = 'no';
save = 1;
opt_ftr=10;
n=40;
mod72

cfl=0.5;
ddt=radius*dxi*cfl/u0;
ndaymax=12;
Tmax=ndaymax*3600*24;
itermax=1000;

%% *** initialisation des données
[ vt_fI,    ht_fI   ] = fun_init(x_fI,   y_fI,   z_fI);
[ vt_fII,   ht_fII  ] = fun_init(x_fII,  y_fII,  z_fII);
[ vt_fIII,  ht_fIII ] = fun_init(x_fIII, y_fIII, z_fIII);
[ vt_fIV,   ht_fIV  ] = fun_init(x_fIV,  y_fIV,  z_fIV);
[ vt_fV,    ht_fV   ] = fun_init(x_fV,   y_fV,   z_fV);
[ vt_fVI,   ht_fVI  ] = fun_init(x_fVI,  y_fVI,  z_fVI);

hte_fI   = ht_fI;
hte_fII  = ht_fII;
hte_fIII = ht_fIII;
hte_fIV  = ht_fIV;
hte_fV   = ht_fV;
hte_fVI  = ht_fVI;

vte_fI   = vt_fI;
vte_fII  = vt_fII;
vte_fIII = vt_fIII;
vte_fIV  = vt_fIV;
vte_fV   = vt_fV;
vte_fVI  = vt_fVI;

%% quantités a conserver
[~,~,~,~,~,~,intref]=...
    nrm72(hte_fI,hte_fII,hte_fIII,hte_fIV,hte_fV,hte_fVI,n,nn,'int');

[ Eref ] = energy(vte_fI,vte_fII,vte_fIII,vte_fIV,vte_fV,vte_fVI,hte_fI,hte_fII,hte_fIII,hte_fIV,hte_fV,hte_fVI,n,nn);
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
    K1v_fI=-(gp*grad_I+cor_I);
    K1v_fII=-(gp*grad_II+cor_II);
    K1v_fIII=-(gp*grad_III+cor_III);
    K1v_fIV=-(gp*grad_IV+cor_IV);
    K1v_fV=-(gp*grad_V+cor_V);
    K1v_fVI=-(gp*grad_VI+cor_VI);
    
    % seconde equation
    % divergence
    [div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
        div72(vt_fI,vt_fII,vt_fIII,vt_fIV,vt_fV,vt_fVI,n,nn);
    
    % assemblage
    K1h_fI=-hp*div_fI;
    K1h_fII=-hp*div_fII;
    K1h_fIII=-hp*div_fIII;
    K1h_fIV=-hp*div_fIV;
    K1h_fV=-hp*div_fV;
    K1h_fVI=-hp*div_fVI;

    %% K2
    % premiere equation

    % GRADIENT
    [grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
        gr72(ht_fI+ddt*K1h_fI,ht_fII+ddt*K1h_fII,ht_fIII+ddt*K1h_fIII,ht_fIV+ddt*K1h_fIV,ht_fV+ddt*K1h_fV,ht_fVI+ddt*K1h_fVI,n,nn);
    % CALCUL DE CORIOLIS
    [cor_I,cor_II,cor_III,cor_IV,cor_V,cor_VI]=coriolis(vt_fI+ddt*K1v_fI,vt_fII+ddt*K1v_fII,vt_fIII+ddt*K1v_fIII,vt_fIV+ddt*K1v_fIV,vt_fV+ddt*K1v_fV,vt_fVI+ddt*K1v_fVI);
    
    % ASSEMBLAGE
    K2v_fI=-(gp*grad_I+cor_I);
    K2v_fII=-(gp*grad_II+cor_II);
    K2v_fIII=-(gp*grad_III+cor_III);
    K2v_fIV=-(gp*grad_IV+cor_IV);
    K2v_fV=-(gp*grad_V+cor_V);
    K2v_fVI=-(gp*grad_VI+cor_VI);
    
    % seconde equation
    % divergence
    [div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
        div72(vt_fI+ddt*K1v_fI,vt_fII+ddt*K1v_fII,vt_fIII+ddt*K1v_fIII,vt_fIV+ddt*K1v_fIV,vt_fV+ddt*K1v_fV,vt_fVI+ddt*K1v_fVI,n,nn);
    
    % assemblage
    K2h_fI=-hp*div_fI;
    K2h_fII=-hp*div_fII;
    K2h_fIII=-hp*div_fIII;
    K2h_fIV=-hp*div_fIV;
    K2h_fV=-hp*div_fV;
    K2h_fVI=-hp*div_fVI;

    %% Assemblage final

    htnew_fI=ht_fI+ddt/2*(K1h_fI+K2h_fI);
    htnew_fII=ht_fII+ddt/2*(K1h_fII+K2h_fII);
    htnew_fIII=ht_fIII+ddt/2*(K1h_fIII+K2h_fIII);
    htnew_fIV=ht_fIV+ddt/2*(K1h_fIV+K2h_fIV);
    htnew_fV=ht_fV+ddt/2*(K1h_fV+K2h_fV);
    htnew_fVI=ht_fVI+ddt/2*(K1h_fVI+K2h_fVI);

    vtnew_fI=vt_fI+ddt/2*(K1v_fI+K2v_fI);
    vtnew_fII=vt_fII+ddt/2*(K1v_fII+K2v_fII);
    vtnew_fIII=vt_fIII+ddt/2*(K1v_fIII+K2v_fIII);
    vtnew_fIV=vt_fIV+ddt/2*(K1v_fIV+K2v_fIV);
    vtnew_fV=vt_fV+ddt/2*(K1v_fV+K2v_fV);
    vtnew_fVI=vt_fVI+ddt/2*(K1v_fVI+K2v_fVI);

    %% calcul de l'erreur
    
    err_fI(iter)=max(max(abs(htnew_fI-hte_fI)./abs(hte_fI)));
    err_fII(iter)=max(max(abs(htnew_fII-hte_fII)./abs(hte_fII)));
    err_fIII(iter)=max(max(abs(htnew_fIII-hte_fIII)./abs(hte_fIII)));
    err_fIV(iter)=max(max(abs(htnew_fIV-hte_fIV)./abs(hte_fIV)));
    err_fV(iter)=max(max(abs(htnew_fV-hte_fV)./abs(hte_fV)));
    err_fVI(iter)=max(max(abs(htnew_fVI-hte_fVI)./abs(hte_fVI)));
    
    err(iter)=max([err_fI(iter),err_fII(iter),err_fIII(iter),err_fIV(iter),err_fV(iter),err_fVI(iter)]);
    time(iter)=t/(24*3600);
    
    
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

    %% conservation
    str='int';
    [~,~,~,~,~,~,int]=...
    nrm72(ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,n,nn,str);
    err_int(iter)=abs(int-intref)./abs(intref);
    
    E = energy(vt_fI,vt_fII,vt_fIII,vt_fIV,vt_fV,vt_fVI,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,n,nn);
    err_energy(iter)=abs(E-Eref)./abs(Eref);

    
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
    mkdir(['./Heun_video-' date ])
    movie2avi(mov, ['./Heun_video-' date '/ref_' num2str(ref) '.avi'], 'compression', 'None');
    
    fid = fopen('AA_VIDEO_SAVE_Heun.txt','a');
    fprintf(fid,'%s\n',['date : ', date]);
    fprintf(fid,'%s\n',['ref. : ', num2str(ref)]);
    fprintf(fid,'%s\n','***********************************');
    fprintf(fid,'%s\n',['test : ', num2str(test)]);
    fprintf(fid,'%s\n','---------- numerical data ---------');
    fprintf(fid,'%s\n',['number of points  : ', num2str(n)] );
    fprintf(fid,'%s\n',['time step         : ', num2str(ddt)] );
    fprintf(fid,'%s\n',['cfl               : ', num2str(cfl)] );
    fprintf(fid,'%s\n',['ordre du filtre   : ', num2str(opt_ftr)] );
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
if save==1
    mkdir(['./Heun_results-' date ]);
    print('-dpng', ['./Heun_results-' date '/ref_' num2str(ref) '_courbe.png'])
    savefig(['./Heun_results-' date '/ref_' num2str(ref) '_courbe']);
end 

figure(2)
plot(time,err)
xlabel('time')
ylabel('erreur relative')
grid on
if save==1
    mkdir(['./Heun_results-' date ]);
    print('-dpng', ['./Heun_results-' date '/ref_' num2str(ref) '_erreur.png'])
    savefig(['./Heun_results-' date '/ref_' num2str(ref) '_erreur']);
end 

figure(3)
semilogy(time,err_int,time,err_energy)
legend('matiere','energie')
xlabel('time')
title('erreur relative sur les conservations')
grid on
if save==1
    mkdir(['./Heun_results-' date ]);
    print('-dpng', ['./RK4_results-' date '/ref_' num2str(ref) '_conservation.png'])
    savefig(['./RK4_results-' date '/ref_' num2str(ref) '_conservation']);
end 

figure(4)
plot_cs11(n,nn,abs(ht_fI-hte_fI)./abs(hte_fI),abs(ht_fII-hte_fII)./abs(hte_fII),abs(ht_fIII-hte_fIII)./abs(hte_fIII),abs(ht_fIV-hte_fIV)./abs(hte_fIV),abs(ht_fV-hte_fV)./abs(hte_fV),abs(ht_fVI-hte_fVI)./abs(hte_fVI));
title('relative error')

figure(5)
plot(time,err_energy)
title('energy')
