clc; clear all; close all; format shorte

c=1/5;

n=100;
h=1./n;
x=[h:h:1]';

cfl=1.5;
ddt=cfl*h/c;
tmax=10;

k=diag(ones(n-1,1),1)-diag(ones(n-1,1),-1);
k(1,end)=-1; k(end,1)=1;
k=sparse(k./(2*h));
p=1/6*diag(ones(n-1,1),1)+1/6*diag(ones(n-1,1),-1)+4/6*eye(n,n);
p(end,1)=1/6; p(1,end)=1/6;
p=sparse(p);


opt_ftr='redonnet6';
if strcmp(opt_ftr,'redonnet10')==1
    ftr0=772/1024;
    ftr1=420/1024;
    ftr2=-240/1024;
    ftr3=90/1024;
    ftr4=-20/1024;
    ftr5=2/1024;
    J=diag(ones(n-1,1),1); J(end,1)=1;
    ftr=ftr0.*eye(n,n)+ftr1/2*(J+J^(n-1))+ftr2/2*(J^2+J^(n-2))+ftr3/2*(J^3+J^(n-3))+ftr4/2*(J^4+J^(n-4))+ftr5/2*(J^5+J^(n-5));
    
elseif strcmp(opt_ftr,'redonnet8')==1
    ftr0=186/256; 
    ftr1=112/256;
    ftr2=-56/256;
    ftr3=16/256;
    ftr4=-2/256;
    ftr5=0;
    J=diag(ones(n-1,1),1); J(end,1)=1;
    ftr=ftr0.*eye(n,n)+ftr1/2*(J+J^(n-1))+ftr2/2*(J^2+J^(n-2))+ftr3/2*(J^3+J^(n-3))+ftr4/2*(J^4+J^(n-4))+ftr5/2*(J^5+J^(n-5));
    
elseif strcmp(opt_ftr,'redonnet6')==1
    ftr0=44/64;
    ftr1=30/64;
    ftr2=-12/64;
    ftr3=2/64;
    ftr4=0;
    ftr5=0;
    J=diag(ones(n-1,1),1); J(end,1)=1;
    ftr=ftr0.*eye(n,n)+ftr1/2*(J+J^(n-1))+ftr2/2*(J^2+J^(n-2))+ftr3/2*(J^3+J^(n-3))+ftr4/2*(J^4+J^(n-4))+ftr5/2*(J^5+J^(n-5));
    
elseif strcmp(opt_ftr,'redonnet4')==1
    ftr0=10/16;
    ftr1=8/16;
    ftr2=-2/16;
    J=diag(ones(n-1,1),1); J(end,1)=1;
    ftr=ftr0.*eye(n,n)+ftr1/2*(J+J^(n-1))+ftr2/2*(J^2+J^(n-2));
elseif strcmp(opt_ftr,'redonnet2')==1
    ftr0=1/2;
    ftr1=1/2;
    J=diag(ones(n-1,1),1); J(end,1)=1;
    ftr=ftr0.*eye(n,n)+ftr1/2*(J+J^(n-1));
elseif strcmp(opt_ftr,'bogey6')==1
    d0=0.234810479761700;
    d1=-.199250131285813;
    d2=0.120198310245186;
    d3=-.049303775636020;
    d4=0.012396449873964;
    d5=-.001446093078167;
    
    ftr0=1-d0;
    ftr1=-d1;
    ftr2=-d2;
    ftr3=-d3;
    ftr4=-d4;
    ftr5=-d5;
    J=diag(ones(n-1,1),1); J(end,1)=1;
    ftr=ftr0.*eye(n,n)+ftr1*(J+J^(n-1))+ftr2*(J^2+J^(n-2))+ftr3*(J^3+J^(n-3))+ftr4*(J^4+J^(n-4))+ftr5*(J^5+J^(n-5));
else
    ftr=speye(n,n);
end
ftr=sparse(ftr);


u=fun(x);
int=sum(u)*h;

e1=[]; e2=[]; ei=[];
cons=[];
t=0;
while t+ddt<tmax
     clc; t=t+ddt
    
    du=k*u;
    w=p\du;
    k1=-c*w;
    
    du=k*(u+ddt/2*k1);
    w=p\du;
    k2=-c*w;
    
    du=k*(u+ddt/2*k2);
    w=p\du;
    k3=-c*w;
    
    du=k*(u+ddt*k3);
    w=p\du;
    k4=-c*w;
    
    fu=u+ddt/6*(k1+2*k2+2*k3+k4);
    u=ftr*fu;
    
    uex=fun(x-c*t);
    e1=[e1 norm(u-uex,1)./norm(uex,1)];
    e2=[e2 norm(u-uex,2)./norm(uex,2)];
    ei=[ei norm(u-uex,inf)./norm(uex,inf)];
    cons=[cons sum(u)*h-int];

    
    pause(0.001)
    clf
    figure(1)
    plot(x,u,'Linewidth',2)
    axis([0 1 -1.2 1.2])
end
time=[1:length(e1)]*ddt;

figure(2)
semilogy(time, 1.4*10^-4*time,time,e2,time,ei,'Linewidth',2)
legend('1.4 \times 10^{-4} \times temps','norme 2','norme \infty','Location','SouthEast')
xlabel('Temps')
ylabel('Erreur')
title(['\lambda = ', num2str(cfl)])
grid on

% figure(3)
% plot([1:length(e1)]*ddt,cons,'Linewidth',2)

    figure(4)
    plot(x,u,'Linewidth',2)
    axis([0 1 -1.2 1.2])

    uex=fun(x);
    figure(5)
    plot(x,uex,x,u,'Linewidth',2)
    axis([0 1 -1.2 1.2])
    xlabel('x')
    title('Filtrage d''ordre 6')

[e1(end) e2(end) ei(end)]

