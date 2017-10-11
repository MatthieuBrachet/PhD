clc; clear all; close all;

teta=linspace(0,pi,100);

se2=sin(teta);
se4=(4/3)*sin(teta)-(1/3)*sin(2*teta)/2;
se6=3/2*sin(teta)-3/5*sin(2*teta)/2+1/10*sin(3*teta)/3;
se8=8/5*sin(teta)-4/5*sin(2*teta)/2+8/35*sin(3*teta)/3-1/35*sin(4*teta)/4;

sc4=sin(teta)./(2/3+2*(1/6*cos(teta)));
sc6=(14/9*sin(teta)+1/18*sin(2*teta))./(1+2*(1/3)*cos(teta));
sc8=((25/16)*sin(teta)+(1/5)*sin(2*teta)/2-(1/80)*sin(3*teta)/3)./(1+2*(3/8)*cos(teta));

a=-103/72; b=91/36; c=-7/4; d=29/36; e=-11/72;
sdc4=(a.*exp(1i*teta)+b.*exp(1i*2*teta)+c.*exp(1i*3*teta)+d.*exp(1i*4*teta)+e.*exp(1i*5*teta))./(2/3+1/6*exp(1i*teta));

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

figure(4)
plot(teta,se4,teta,sc4,teta,imag(sdc4),teta,teta,'k-','Linewidth',2)
legend('expli. 4','comp. 4','decentré. 4','theoric.','Location','northwest')
xlabel('\theta')
ylabel('modified \theta : \theta_m')
grid on

figure(5)
plot(real(sdc4),'Linewidth',2)
grid on