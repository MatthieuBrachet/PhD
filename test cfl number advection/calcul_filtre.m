clc; clear all; close all

theta = 0:0.01:pi;
f=10;
filtre = 'visbal';
delta=0.4;
[ ftr1 ] = ftr( theta,f,filtre,delta );
n=20; ftr2 = 1-sin(theta/2).^n;
[ ftr3 ] = ftr( theta,2,'redonnet',0 );
[ ftr4 ] = ftr( theta,4,'redonnet',0 );
[ ftr5 ] = ftr( theta,6,'redonnet',0 );
[ ftr6 ] = ftr( theta,8,'redonnet',0 );
[ ftr7 ] = ftr( theta,10,'redonnet',0 );


figure(1)
plot(theta,ftr1,'x-',theta, ftr2, 'x-',theta,ftr3,theta,ftr4,theta,ftr5,theta,ftr6,theta,ftr7)
legend(['imp. \delta=', num2str(delta)],['Shapiro n=', num2str(n)],'expl. 2','expl. 4','expl. 6','expl. 8','expl. 10') 
grid on