clc; clear all; close all;

teta=linspace(0,pi,100);
se2=sin(teta);
se4=(4/3)*sin(teta)-(1/3)*sin(2*teta)/2;
se6=3/2*sin(teta)-3/5*sin(2*teta)/2+1/10*sin(3*teta)/3;
sc4=sin(teta)./(2/3+2*(1/6*cos(teta)));
sc6=(14/9*sin(teta)+1/18*sin(2*teta))./(1+(2/3)*cos(teta));

figure(1)
plot(teta,se2,'b-',teta,se4,'r-',teta,se6,'g-',teta,sc4,'m-',teta,sc6,'-',teta,teta,'k-','Linewidth',2)
legend('expl. 2','expl. 4','expl. 6','compact 4','compact 6','theoric.','Location','northwest')
xlabel('\theta= \Deltax \xi')
ylabel('modified \theta : \theta_m')
grid on