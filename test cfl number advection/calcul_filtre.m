clc; clear all; close all

theta = 0:0.01:pi;
f=10;
filtre = 'visbal';
delta=0.48;

[ ftr3 ] = ftr( theta,2,'redonnet',0 );
[ ftr4 ] = ftr( theta,4,'redonnet',0 );
[ ftr5 ] = ftr( theta,6,'redonnet',0 );
[ ftr6 ] = ftr( theta,8,'redonnet',0 );
[ ftr7 ] = ftr( theta,10,'redonnet',0 );


figure(1)
plot(theta,ftr3,theta,ftr4,theta,ftr5,theta,ftr6,theta,ftr7)
legend('expl. 2','expl. 4','expl. 6','expl. 8','expl. 10') 
xlabel('\theta')
ylabel('F(\theta)')
grid on