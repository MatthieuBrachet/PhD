clc; clear all; close all;

teta=linspace(0,pi,10000);

J1=5;
[ ampli1 ] = ampli_ftr( J1, teta );
max(ampli1)
J2=10;
[ ampli2 ] = ampli_ftr( J2, teta );
max(ampli2)
J3=15;
[ ampli3 ] = ampli_ftr( J3, teta );
max(ampli3)
J4=20;
[ ampli4 ] = ampli_ftr( J4, teta );
max(ampli4)
J5=30;
[ ampli5 ] = ampli_ftr( J5, teta );
max(ampli5)


figure(1)
plot(teta,ampli1,teta,ampli2,teta,ampli3,teta,ampli4,teta,ampli5,'Linewidth',2)
legend(['order ' num2str(2*J1)],['order ' num2str(2*J2)],['order ' num2str(2*J3)],['order ' num2str(2*J4)],['order ' num2str(2*J5)],'Location','southwest');
title('Filter amplification')