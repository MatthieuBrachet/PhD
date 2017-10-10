clc; clear all; close all;

teta=linspace(0,pi,100);

se2=sin(teta);
se4=(4/3)*sin(teta)-(1/3)*sin(2*teta)/2;
se6=3/2*sin(teta)-3/5*sin(2*teta)/2+1/10*sin(3*teta)/3;
se8=8/5*sin(teta)-4/5*sin(2*teta)/2+8/35*sin(3*teta)/3-1/35*sin(4*teta)/4;

sc4=sin(teta)./(2/3+2*(1/6*cos(teta)));
sc6=(14/9*sin(teta)+1/18*sin(2*teta))./(1+(2/3)*cos(teta));
sc8=((25/16)*sin(teta)+(1/5)*sin(2*teta)/2-(1/80)*sin(3*teta)/3)./(1+2*(3/8)*cos(teta)); !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! pbm here!!!!

alpha=.435181352;
a=1.551941906;
b=.361328195;
c=-.042907397;
sk4=((a)*sin(teta)+(b)*sin(2*teta)/2+(c)*sin(3*teta)/3)./(1+2*(alpha)*cos(teta));

figure(1)
plot(teta,se2,'--',teta,se4,'--',teta,se6,'--',teta,se8,'--',teta,sc4,'-',teta,sc6,'-',teta,sc8,'-',teta,teta,'k-','Linewidth',2)
legend('expl. 2','expl. 4','expl. 6','expl. 8','comp. 4','comp. 6','comp. 8','theoric.','Location','northwest')
xlabel('\theta')
ylabel('modified \theta : \theta_m')
grid on

figure(2)
plot(teta,se2,'-',teta,se4,'-',teta,se6,'-',teta,se8,'-',teta,teta,'k-','Linewidth',2)
legend('expl. 2','expl. 4','expl. 6','expl. 8','theoric.','Location','northwest')
xlabel('\theta')
ylabel('modified \theta : \theta_m')
grid on

figure(3)
plot(teta,se4,teta,sc4,teta,sk4,teta,teta,'k-','Linewidth',2)
legend('expli. 4','comp. 4','opt. 4','theoric.','Location','northwest')
xlabel('\theta')
ylabel('modified \theta : \theta_m')
grid on
