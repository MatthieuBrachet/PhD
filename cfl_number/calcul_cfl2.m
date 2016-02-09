clc; clear all; close all

tt=-pi:0.001:pi;
N=floor(length(tt)/2)+1;
cfl=1;
c=4;
space='implicit';
f=10;
filtre='redonnet';
delta=0.1;

%% schemas explicites

% *** RK1
time1='dirk12a';
c1 = atsf(cfl, tt,time1,c,space,f,filtre,delta);
g1=log(c1);
h1=abs(c1);
h1=h1(N:end);
a1=angle(c1)./(-cfl*tt);
a1=a1(N:end);

% *** RK2
time2='dirk23a';
c2 = atsf(cfl, tt,time2,c,space,f,filtre,delta);
g2=log(c2);
h2=abs(c2);
h2=h2(N:end);
a2=angle(c2)./(-cfl*tt);
a2=a2(N:end);

% *** RK4
time3='dirk34a';
c3 = atsf(cfl, tt,time3,c,space,f,filtre,delta);
g3=log(c3);
h3=abs(c3);
h3=h3(N:end);
a3=angle(c3)./(-cfl*tt);
a3=a3(N:end);

%% schemas implicites

% *** DIRK12
time4='dirk22b';
c4 = atsf(cfl, tt,time4,c,space,f,filtre,delta);
g4=log(c4);
h4=abs(c4);
h4=h4(N:end);
a4=angle(c4)./(-cfl*tt);
a4=a4(N:end);

% *** DIRK23
time5='dirk33b';
c5 = atsf(cfl, tt,time5,c,space,f,filtre,delta);
g5=log(c5);
h5=abs(c5);
h5=h5(N:end);
a5=angle(c5)./(-cfl*tt);
a5=a5(N:end);

% *** DIRK34
time6='rk4';
c6 = atsf(cfl, tt,time6,c,space,f,filtre,delta);
g6=log(c6);
h6=abs(c6);
h6=h6(N:end);
a6=angle(c6)./(-cfl*tt);
a6=a6(N:end);

%% figures

figure(1)
hold on
grid on
plot(g1,'b-x')
plot(g2,'r-x')
plot(g3,'g-x')
plot(g4,'k-o')
plot(g5,'c-o')
plot(g6,'m-o')
legend(time1,time2,time3,time4,time5,time6,2)
hold off

figure(2)
subplot(121)
hold on
grid on
plot(tt(N:end),h1,'b-')
plot(tt(N:end),h2,'r-')
plot(tt(N:end),h3,'g-')
plot(tt(N:end),h4,'k--')
plot(tt(N:end),h5,'c--')
plot(tt(N:end),h6,'m--')
legend(time1,time2,time3,time4,time5,time6,3)
axis([0 pi 0 1.1])
title('Dissipation')
hold off

subplot(122)
hold on
grid on
plot(tt(N:end),a1,'b-')
plot(tt(N:end),a2,'r-')
plot(tt(N:end),a3,'g-')
plot(tt(N:end),a4,'k--')
plot(tt(N:end),a5,'c--')
plot(tt(N:end),a6,'m--')
legend(time1,time2,time3,time4,time5,time6,3)
title('Dispersion')
axis([0 pi 0 1.1])
hold off
