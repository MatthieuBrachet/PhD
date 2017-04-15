clc; clear all; close all

theta = linspace(0,pi,50);
f=10;
filtre = 'redonnet';
delta=0.48;

[ ftr3 ] = ftr( theta,2,filtre,0 );
[ ftr4 ] = ftr( theta,4,filtre,0 );
[ ftr5 ] = ftr( theta,6,filtre,0 );
[ ftr6 ] = ftr( theta,8,filtre,0 );
[ ftr7 ] = ftr( theta,10,filtre,0 );


figure(1)
plot(theta,ftr3,theta,ftr4,theta,ftr5,theta,ftr6,theta,ftr7,theta,ones(size(theta)),'k-','linewidth',2)
legend('expl. 2','expl. 4','expl. 6','expl. 8','expl. 10','theoric.') 
xlabel('\theta')
ylabel('F(\theta)')
grid on