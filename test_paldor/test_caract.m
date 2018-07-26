clc; clear all; close all

Nt=150;
t=linspace(0,2,Nt);
Nx=100;
x=linspace(0,2*pi,Nx);
[X,T]=meshgrid(x,t);

c=2.5
Y=X-c*T;
U=cos(Y).*cos(5*Y)+sin(7*Y)+.5*rand(size(Y));
Cn= fit_hov( X,T,U )

figure(1)
hold on
surf(X,T,U)
shading interp
view(2)
axis([0 2*pi 0 1])
colorbar
xlabel('space')
ylabel('time')
plot3(x,x/c,2.5*ones(size(x)),'r-','Linewidth',2)
plot3(x,x/Cn,2.5*ones(size(x)),'r--','Linewidth',2)
hold off