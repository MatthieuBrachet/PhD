clc; clear all; close all

%% paramétres numériques
epsi=10^-10;
itermax=10000;
iter=0;
qual=0.95;

%% schema testé
order=2;
genre='redonnet';
delta=0;

%% algorithme
a=0;
b=pi;
while abs(b-a)>epsi & iter<itermax
    theta=0.5*(a+b);
    m=ftr(theta,order,genre,delta);
    iter=iter+1;
    if m>qual
        a=theta;
        b=b;
    elseif m<qual
        a=a;
        b=theta;
    end
end

theta
theta/pi*100
iter

%% courbe
teta=0:0.001:pi;
ff=ftr(teta,order,genre,delta);

figure(1)
plot(teta,ff,'b-',theta,qual,'ro')
grid on
xlabel('freq.')
ylabel('quality')
title(['filter of ', genre, ' with ', num2str(delta) ,' and order ', num2str(order)])
