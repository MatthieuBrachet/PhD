clc; clear all; close all

theta = linspace(0,pi,50);
f=10;
filtre = 'visbal';
delta=0.0;

[ ftr3 ] = ftr( theta,2,filtre,delta );
[ ftr4 ] = ftr( theta,4,filtre,delta );
[ ftr5 ] = ftr( theta,6,filtre,delta );
[ ftr6 ] = ftr( theta,8,filtre,delta );
[ ftr7 ] = ftr( theta,10,filtre,delta );


figure(1)
plot(theta,ftr3,theta,ftr4,theta,ftr5,theta,ftr6,theta,ftr7,theta,ones(size(theta)),'k-','linewidth',2)
legend('expl. 2','expl. 4','expl. 6','expl. 8','expl. 10','theorique') 
xlabel('\theta')
ylabel('F(\theta)')
xticks([0 pi/4 pi/2 3*pi/4 pi])
xticklabels({'0','\pi/4','\pi/2','3\pi /4','\pi'})
grid on