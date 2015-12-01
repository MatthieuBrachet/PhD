clc; clear all; close all

% test 1
% observons si la solution est bien stationnaire

global n nn;
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI;
global opt_ftr;
global g radius OMEGA 
global gamma u0 phi0 h0 phi1 umax
global n_simpson
global test



%% *** choix du test ******************************************************
test = 0;
% test = 0 test court de Galewsky, Scott et Polvani.
%        1 test long de Galewsky, Scott et Polvani.



% *************************************************************************

%% discretisation data
n=10;
n_simpson=1000;

%% physical data
g=9.80616;
radius=6.37122*10^6;
opt_ftr=10;
OMEGA=7.292*10^-5;

%% functional data
% test court de Galewski, Scott et Polvani
gamma=pi/18;
u0=80;
phi0=pi/4;
h0=1;

% test long de Galewski, Scott et Polvani
umax=80;
phi0=pi/10;
phi1=pi/2-phi0;

%% ************************************************************************
data_mod;

%% calcul de la vitesse
disp('1- calcul du champ de vitesse...')
V_fI=vitesse(x_fI,y_fI,z_fI);
V_fII=vitesse(x_fII,y_fII,z_fII);
V_fIII=vitesse(x_fIII,y_fIII,z_fIII);
V_fIV=vitesse(x_fIV,y_fIV,z_fIV);
V_fV=vitesse(x_fV,y_fV,z_fV);
V_fVI=vitesse(x_fVI,y_fVI,z_fVI);

%% calcul du geopotentiel gp
disp('2- calcul du geopotentiel...')
[ gp_fI ] = geopotential( x_fI,y_fI,z_fI );disp('   face I   : ok');
[ gp_fII ] = gp_fI;disp('   face II  : ok');
[ gp_fIII ] = gp_fI;disp('   face III : ok');
[ gp_fIV ] = gp_fI;disp('   face IV  : ok');
[ gp_fV ] = geopotential( x_fV,y_fV,z_fV );disp('   face V   : ok');
[ gp_fVI ] = geopotential( x_fVI,y_fVI,z_fVI );disp('   face VI  : ok');

figure(1)
plot_cs5(n,nn,gp_fI,gp_fII,gp_fIII,gp_fIV,gp_fV,gp_fVI)

%% *** equation de moment *************************************************

%% calcul de la convection V . nabla V 
disp('3- calcul de la convection...')
for p=1:3
    [gr_I,gr_II,gr_III,gr_IV,gr_V,gr_VI]=gr(V_fI(:,:,p),V_fII(:,:,p),V_fIII(:,:,p),V_fIV(:,:,p),V_fV(:,:,p),V_fVI(:,:,p),n,nn);
    convec_fI(:,:,p) = V_fI(:,:,1).*gr_I(:,:,1)+V_fI(:,:,2).*gr_I(:,:,2)+V_fI(:,:,3).*gr_I(:,:,3);
    convec_fII(:,:,p) = V_fII(:,:,1).*gr_II(:,:,1)+V_fII(:,:,2).*gr_II(:,:,2)+V_fII(:,:,3).*gr_II(:,:,3);
    convec_fIII(:,:,p) = V_fIII(:,:,1).*gr_III(:,:,1)+V_fIII(:,:,2).*gr_III(:,:,2)+V_fIII(:,:,3).*gr_III(:,:,3);
    convec_fIV(:,:,p) = V_fIV(:,:,1).*gr_IV(:,:,1)+V_fIV(:,:,2).*gr_IV(:,:,2)+V_fIV(:,:,3).*gr_IV(:,:,3);
    convec_fV(:,:,p) = V_fV(:,:,1).*gr_V(:,:,1)+V_fV(:,:,2).*gr_V(:,:,2)+V_fV(:,:,3).*gr_V(:,:,3);
    convec_fVI(:,:,p) = V_fVI(:,:,1).*gr_VI(:,:,1)+V_fVI(:,:,2).*gr_VI(:,:,2)+V_fVI(:,:,3).*gr_VI(:,:,3);
end

%% calcul de la force de Coriolis f.kxV (attention au signe!).
disp('4- calcul de la force de Coriolis...')
cor_fI = coriolis(V_fI,x_fI,y_fI,z_fI);
cor_fII = coriolis(V_fII,x_fII,y_fII,z_fII);
cor_fIII = coriolis(V_fIII,x_fIII,y_fIII,z_fIII);
cor_fIV = coriolis(V_fIV,x_fIV,y_fIV,z_fIV);
cor_fV = coriolis(V_fV,x_fV,y_fV,z_fV);
cor_fVI = coriolis(V_fVI,x_fVI,y_fVI,z_fVI);

%% influence du geopotentiel : nabla gp
disp('5- influence du geopotentiel...')
[gr_fI,gr_fII,gr_fIII,gr_fIV,gr_fV,gr_fVI]=gr(gp_fI,gp_fII,gp_fIII,gp_fIV,gp_fV,gp_fVI,n,nn);

%% second membre total :
disp('6- bilan...')
F_fI=-convec_fI-cor_fI-gr_fI;
F_fII=-convec_fII-cor_fII-gr_fII;
F_fIII=-convec_fIII-cor_fIII-gr_fIII;
F_fIV=-convec_fIV-cor_fIV-gr_fIV;
F_fV=-convec_fV-cor_fV-gr_fV;
F_fVI=-convec_fVI-cor_fVI-gr_fVI;

%% residu
disp('7- calcul du residu...')
e_I=max(max(max(abs(F_fI))));
e_II=max(max(max(abs(F_fII))));
e_III=max(max(max(abs(F_fIII))));
e_IV=max(max(max(abs(F_fIV))));
e_V=max(max(max(abs(F_fV))));
e_VI=max(max(max(abs(F_fVI))));

e1=max([e_I,e_II,e_III,e_IV,e_V,e_VI]);

%% *** equation de conservation *******************************************
disp('8- equation de conservation...')
%% membre 1
mem1_fI=V_fI(:,:,1).*gr_fI(:,:,1)+V_fI(:,:,2).*gr_fI(:,:,2)+V_fI(:,:,3).*gr_fI(:,:,3);
mem1_fII=V_fII(:,:,1).*gr_fII(:,:,1)+V_fII(:,:,2).*gr_fII(:,:,2)+V_fII(:,:,3).*gr_fII(:,:,3);
mem1_fIII=V_fIII(:,:,1).*gr_fIII(:,:,1)+V_fIII(:,:,2).*gr_fIII(:,:,2)+V_fIII(:,:,3).*gr_fIII(:,:,3);
mem1_fIV=V_fIV(:,:,1).*gr_fIV(:,:,1)+V_fIV(:,:,2).*gr_fIV(:,:,2)+V_fIV(:,:,3).*gr_fIV(:,:,3);
mem1_fV=V_fV(:,:,1).*gr_fV(:,:,1)+V_fV(:,:,2).*gr_fV(:,:,2)+V_fV(:,:,3).*gr_fV(:,:,3);
mem1_fVI=V_fVI(:,:,1).*gr_fVI(:,:,1)+V_fVI(:,:,2).*gr_fVI(:,:,2)+V_fVI(:,:,3).*gr_fVI(:,:,3);

%% membre 2
[div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=div(V_fI,V_fII,V_fIII,V_fIV,V_fV,V_fVI,n,nn);
mem2_fI=gp_fI.*div_fI;
mem2_fII=gp_fII.*div_fII;
mem2_fIII=gp_fIII.*div_fIII;
mem2_fIV=gp_fIV.*div_fIV;
mem2_fV=gp_fV.*div_fV;
mem2_fVI=gp_fVI.*div_fVI;

G_fI=-mem1_fI-mem2_fI;
G_fII=-mem1_fII-mem2_fII;
G_fIII=-mem1_fIII-mem2_fIII;
G_fIV=-mem1_fIV-mem2_fIV;
G_fV=-mem1_fV-mem2_fV;
G_fVI=-mem1_fVI-mem2_fVI;

figure(2)
plot_cs5(n,nn,G_fI,G_fII,G_fIII,G_fIV,G_fV,G_fVI)

%% residu
disp('9- calcul du residu.')
e_I=max(max(max(abs(G_fI))));
e_II=max(max(max(abs(G_fII))));
e_III=max(max(max(abs(G_fIII))));
e_IV=max(max(max(abs(G_fIV))));
e_V=max(max(max(abs(G_fV))));
e_VI=max(max(max(abs(G_fVI))));

e2=max([e_I,e_II,e_III,e_IV,e_V,e_VI]);



