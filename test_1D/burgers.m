clc; clear all; close all;

n=200;
h=1./n;
x=[h:h:1]';

cfl=.8;
ddt=cfl*h/.5;
tmax=1;

k=diag(ones(n-1,1),1)-diag(ones(n-1,1),-1);
k(1,end)=-1; k(end,1)=1;
k=sparse(k./(2*h));
p=1/6*diag(ones(n-1,1),1)+1/6*diag(ones(n-1,1),-1)+4/6*eye(n,n);
p(end,1)=1/6; p(1,end)=1/6;
p=sparse(p);

opt_ftr='redonnet10';
if strcmp(opt_ftr,'redonnet10')==1
    ftr0=772/1024;
    ftr1=420/1024;
    ftr2=-240/1024;
    ftr3=90/1024;
    ftr4=-20/1024;
    ftr5=2/1024;
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
else
    ftr=speye(n,n);
end

u=.5+.4*sin(2*pi*x);
%u=exp(-500*(x-.4).^2);
int=sum(u)*h;

e=[];cons=[];
t=0;
while t+ddt<tmax
     clc; t=t+ddt
    
    fu=u;
    du=k*(fu.^2);
    w=p\du;
    k1=-.5*w;
    
    fu=u+ddt/2.*k1;
    du=k*(fu.^2);
    w=p\du;
    k2=-.5*w;
    
    fu=u+ddt/2.*k2;
    du=k*(fu.^2);
    w=p\du;
    k3=-.5*w;
    
    fu=u+ddt.*k3;
    du=k*(fu.^2);
    w=p\du;
    k4=-.5*w;
    
    fu=u+ddt/6*(k1+2*k2+2*k3+k4);
    u=ftr*fu;
    
    cons=[cons sum(u)*h-int];
    
%     w=k*u;
%     dux=p\w;
%     w=k*dux;
%     duxx=p\w;
    
    pause(0.00001)
    clf
    figure(1)
    plot(x,u)
    axis([0 1 0 1])
end

figure(2)
plot(cons)
