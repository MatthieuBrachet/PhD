clc; clear all; close all; format shorte

global n X Y
global gp hp coriolis
global opt_ftr
global FTR

opt_ftr='redonnet100';

n=64;
mod100
gg=6*n^2*gp*hp+coriolis^2;
ddt=2*sqrt(2)/sqrt(gg);
ddt=1e-4;


[h,u,v] = fun100(X,Y);
cons=[];
ener=[];
time=[];
consref=sum(h)/(n^2);
energyref=energy(h,u,v);
t=0;tmax=2;
while t<tmax
    t=t+ddt;
    time=[time t];
    
    %% K1
    hh=h;
    uu=u;
    vv=v;
    [kh1,ku1,kv1] = semi_disc(hh,uu,vv);
    
    %% K2
    hh=h+ddt/2*kh1;
    uu=u+ddt/2*ku1;
    vv=v+ddt/2*kv1;
    [kh2,ku2,kv2] = semi_disc(hh,uu,vv);
    
    %% K3
    hh=h+ddt/2*kh2;
    uu=u+ddt/2*ku2;
    vv=v+ddt/2*kv2;
    [kh3,ku3,kv3] = semi_disc(hh,uu,vv);
    
    %% K4
    hh=h+ddt*kh3;
    uu=u+ddt*ku3;
    vv=v+ddt*kv3;
    [kh4,ku4,kv4] = semi_disc(hh,uu,vv);
    
    %% Assemblage
    h=h+ddt/6*(kh1+2*kh2+2*kh3+kh4);
    u=u+ddt/6*(ku1+2*ku2+2*ku3+ku4);
    v=v+ddt/6*(kv1+2*kv2+2*kv3+kv4);
    
    %% Filtrage
    h=FTR*h;
    u=FTR*u;
    v=FTR*v;
    
    %% Conservation
    cons=[cons sum(h)/(n^2)];
    ener=[ener energy(h,u,v)];
    
%     %% plot
%     hplot=reshape(h,n,n);
%     xplot=reshape(X,n,n);
%     yplot=reshape(Y,n,n);
%     
%     pause(1e-10)
%     clf
%     figure(1)
%     surf(xplot,yplot,hplot)
%     shading interp;
%     axis([0 1 0 1 -1/3 1])
%     alpha 0.7
%     colormap winter

    clc; [t ener(end)/energyref-1]
end

figure(2)
plot(time,cons/consref-1,'Linewidth',2)
title('Erreur sur la conservation de la masse')


figure(3)
plot(time,ener./energyref-1,'Linewidth',2)
title('Erreur sur la conservation de l''energie')
grid on

%% plot
hplot=reshape(h,n,n);
xplot=reshape(X,n,n);
yplot=reshape(Y,n,n);

figure(4)
surf(xplot,yplot,hplot)
shading interp;
axis([0 1 0 1 -1/3 1])
%alpha 0.7
colormap jet
colorbar
caxis([-0.5 1])

hFig=figure(5);
contourf(xplot,yplot,hplot)
set(hFig, 'Position', [50 50 500 500])
colorbar


max(abs(ener./energyref-1))

