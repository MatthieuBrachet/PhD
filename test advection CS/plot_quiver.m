function []=plot_quiver(nn,vitx_I,vitx_II, vitx_III, vitx_IV, vitx_V, vitx_VI, ...
    vity_I,vity_II, vity_III, vity_IV, vity_V, vity_VI, ...
    vitz_I,vitz_II, vitz_III, vitz_IV, vitz_V, vitz_VI, time)
% PLOT VECTORS THE CUBED SPHERE GRID
% M. Brachet - SEPT 9, 2015

global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;
global radius;

global coef;
global u0 lambda0 teta0;

%% PLOT OF A QUIVER ON THE SPHERE

% FACE F - I;
quiver3(x_fI, y_fI, z_fI,vitx_I, vity_I, vitz_I,'color',[1,0,0]);hold on;
axis([-radius radius -radius radius -radius radius]);
text(radius,0,0,'F','Color','g','FontSize',15);

% FACE E - II;
quiver3(x_fII, y_fII, z_fII,vitx_II, vity_II, vitz_II,'color',[1,0,0]);hold on;
axis([-radius radius -radius radius -radius radius]);
text(0,radius,0,'E','Color','g','FontSize',15);

% FACE B - III;
quiver3(x_fIII, y_fIII, z_fIII,vitx_III, vity_III, vitz_III,'color',[1,0,0]);hold on;
axis([-radius radius -radius radius -radius radius]);
text(-radius,0,0,'B','Color','g','FontSize',15);

% FACE W - IV;
quiver3(x_fIV, y_fIV, z_fIV,vitx_IV, vity_IV, vitz_IV,'color',[1,0,0]);hold on;
axis([-radius radius -radius radius -radius radius]);
text(0,-radius,0,'W','Color','g','FontSize',15);

% FACE N - V;
quiver3(x_fV, y_fV, z_fV,vitx_V, vity_V, vitz_V,'color',[1,0,0]);hold on;
axis([-radius radius -radius radius -radius radius]);
text(0,0,radius,'N','Color','g','FontSize',15);

% FACE S - VI;
quiver3(x_fVI, y_fVI, z_fVI,vitx_VI, vity_VI, vitz_VI,'color',[1,0,0]);hold on;
axis([-radius radius -radius radius -radius radius]);
text(0,0,-radius,'S','Color','g','FontSize',15);

axis equal;
axis([-radius radius -radius radius -radius radius]);

ws=u0/radius;
if coef == 2
%% position du centre du vortex
[lambda0_prime,teta0_prime ] = rotated_coord( lambda0, teta0 );
lambdac_prime=lambda0_prime+ws*time;
tetac_prime=teta0_prime;
[lambdac,tetac ] = unrotated_coord( lambdac_prime, tetac_prime );
[xc,yc,zc]=sph2cart(lambdac,tetac,radius);
plot3(xc,yc,zc,'ko'); hold on
end
funfIe=fun4_b(x_fI,y_fI,z_fI,time);funfIIe=fun4_b(x_fII,y_fII,z_fII,time);
funfIIIe=fun4_b(x_fIII,y_fIII,z_fIII,time); funfIVe=fun4_b(x_fIV,y_fIV,z_fIV,time);
funfVe=fun4_b(x_fV,y_fV,z_fV,time); funfVIe=fun4_b(x_fVI,y_fVI,z_fVI,time);
plot_cs5(nn-2,nn,funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe);colorbar; hold on;
hold off

