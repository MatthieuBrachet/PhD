clc; clear all; close all;
% Specify x range and number of points
x0 = -3;
x1 = 3;
Nx = 5000;
% Specify y range and number of points
y0 = -3;
y1 = 3;
Ny = 5000;
% Construct mesh
xv = linspace(x0,x1,Nx);
yv = linspace(y0,y1,Ny);
[x,y] = meshgrid(xv,yv);
% Calculate z
z = x + 1i*y;
% 2nd order Runge-Kutta growth factor
g = 1 + z + 0.5*z.^2 + (1/6)*z.^3 + (1/24)*z.^4;
% Calculate magnitude of g
gmag = abs(g);
% Plot contours of gmag
hFig=figure(2)
contourf(x,y,(gmag<1),[1 1],'k-','LineWidth',2.0);
grid minor
colormap winter
axis([-3 0.5 -3 3])
yticks([-2*sqrt(2) -sqrt(2) 0 sqrt(2) 2*sqrt(2)])
yticklabels({'-2\surd{2}','-\surd{2}','0','\surd{2}','2\surd{2}'})
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [50 50 350 600])