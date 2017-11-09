clc; clear all; close all;

c=1/5;

n=200;
h=1./n;
x=[h:h:1]';

cfl=2*sqrt(2/3);
ddt=cfl*h/c;
tmax=8;

k=diag(ones(n-1,1),1)-diag(ones(n-1,1),-1);
k(1,end)=-1; k(end,1)=1;
k=sparse(k./(2*h));
p=1/6*diag(ones(n-1,1),1)+1/6*diag(ones(n-1,1),-1)+4/6*eye(n,n);
p(end,1)=1/6; p(1,end)=1/6;
p=sparse(p);


opt_ftr='redonnet';
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
ftr=sparse(ftr);


u=(abs(x-.3)<0.2);
%u=exp(-500*(x-.3).^2);
%u=.4*cos(2*pi*x)+.5;
int=sum(u)*h;

e=[]; cons=[];
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
    
    %uex=exp(-500*(x-c*t-.3).^2);
    uex=(abs(x-c*t-.3)<0.2);
    %uex=.4*cos(2*pi*(x-c*t))+.5;
    e=[e norm(u-uex,inf)];
    cons=[cons sum(u)*h-int];
    
    w=k*u;
    dux=p\w;
    w=k*dux;
    duxx=p\w;
    
    pause(0.001)
    clf
    figure(1)
    plot(x,log(abs(duxx)))
   % axis([0 1 -.2 1.2])
end


figure(2)
semilogy(e)

figure(3)
plot(cons)

