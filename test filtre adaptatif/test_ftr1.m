clc; clear all; close all;
global scheme

N=100;
h=1./N;
x=[h:h:1]';

%% initial function
u=10*(x<0.5)+0.*rand(size(x));

%% classic filter
opt_ftr1='redonnet4';
[ ftr ] = filtre74( N , opt_ftr1 );
uf1=ftr*u;

%% adaptative filter
opt_ftr2='redonnet2';
scheme='compact4';
[ LAP, B, ftra, p, k ] = adaptative74( N, opt_ftr2 );

% detection
u1=ftr*u;
u1=LAP*u1;
u1=p\(k*u1);
dmag=0.5*((B*u1).^2+((B')*u1).^2);
epsilon=10^-16;
c=u1+(u1==0)*epsilon;
r=dmag./(c.^2./(h.^2))+epsilon;
rth=10^-5;
sigma=0.5*(1-rth./r+abs(1-rth./r));
uf2=sigma.*(ftra*u1)+(1-sigma).*u1;


figure(1)
plot(x,u,'^-',x,uf1,'x-',x,uf2,'o-')
legend('initial function','classic filter','adaptative filter')
grid on

figure(2)
semilogy(x,abs(u-uf1)+10^-16,x,abs(u-uf2)+10^-16)
legend('classic filter','adaptative filter')
grid on

figure(3)
plot(x,sigma,'x-')
legend('detecteur')
grid on


