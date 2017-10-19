clc; clear all; close all;

teta=linspace(0,pi,10000);

J1=1;
ampli1=1-(1/(-2)^J1).*(cos(teta)-1).^J1;

J2=2;
ampli2=1-(1/(-2)^J2).*(cos(teta)-1).^J2;

J3=3;
ampli3=1-(1/(-2)^J3).*(cos(teta)-1).^J3;

J4=4;
ampli4=1-(1/(-2)^J4).*(cos(teta)-1).^J4;

J5=5;
ampli5=1-(1/(-2)^J5).*(cos(teta)-1).^J5;

figure(1)
plot(teta,ampli1,teta,ampli2,teta,ampli3,teta,ampli4,teta,ampli5,'Linewidth',1.5)
legend(['order ' num2str(2*J1)],['order ' num2str(2*J2)],['order ' num2str(2*J3)],['order ' num2str(2*J4)],['order ' num2str(2*J5)],'Location','southwest');
%title('Filter amplification')
grid on
xlabel('\theta')
ylabel('\beta(\theta)')

for J=1:5
    teta_m=acos(1-2*(.05)^(1/J))
end

J=linspace(1,100,50);
teta_m=acos(1-2.*(.05).^(1./J));

figure(2)
plot(J,teta_m)