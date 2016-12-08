clc; clear all; close all;

n=1000;
dxi=(1./n);
x=[0:dxi:1-dxi]';

% filtre adaptatif
opt_det='redonnet6';
opt_ftra='redonnet4';
[ detec ] = filtre74( n , opt_det );
[ LAP_adap, B_adap, ftra ] = adaptative74( n, opt_ftra );
epsilon=10^-16;
rth=100;

% filtre conventionnel
opt_ftr4='redonnet4';
[ ftr4 ] = filtre74( n , opt_ftr4 );

opt_ftr10='redonnet10';
[ ftr10 ] = filtre74( n , opt_ftr10 );

Ea=[];
E4=[];
E10=[];
N=[];
for i=1:floor(n/2)-1
    u=cos(2*pi*i*x);
    
    %% filtre adaptatif
    % detection
    uhf=u-detec*u;
    r=0.5*((B_adap*uhf).^2+((B_adap')*uhf).^2)./dxi.^2+epsilon;
    sigma=0.5*(1-rth./r+abs(1-rth./r));
    
    S(i)=max(sigma);
    
    % filtrage
    uf=ftra*u;
    ufa=(1-sigma).*(u)+sigma.*uf;
    
    e=norm(ufa,2)./norm(u,2);
    Ea=[Ea e];
    
    %% filtre conventionnel
    uf4=ftr4*u;
    e=norm(uf4,2)./norm(u,2);
    E4=[E4 e];

    %% filtre conventionnel
    uf10=ftr10*u;
    e=norm(uf10,2)./norm(u,2);
    E10=[E10 e];

    N=[N i*2*pi/n];
end

figure(1)
plot(N,Ea,'xk-',N,E4,'b--',N,E10,'m--',N,S,'r.-')
legend('filtrage adaptatif',opt_ftr4,opt_ftr10,'detecteur')
grid on
