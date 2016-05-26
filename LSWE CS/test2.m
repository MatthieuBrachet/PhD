% *************************************************************************
% test Linear Shallow Water
% 
% author :
%     - Matthieu Brachet
% *************************************************************************
clc; clear all; close all
%% *** global data ********************************************************
global nn n dxi radius u0;
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr test;

%% *** options ************************************************************
opt_ftr = 2;
test = 0; % (test=0 : personal test - 0 : divergence nulle; test=1 : galewsky regular; test=2 : galewsky par morceaux)
video = 'no';
ppp=4;

%% *** physical data ******************************************************

%% *** space data *********************************************************
n=40;
nn=n+2;
mod72;

%% *** time data **********************************************************
cfl=0.1;
u0=80;
ddt=cfl*radius*dxi/u0;
Tmax=12;
tmax=Tmax*86400;
itemax=0;
%% BLOC 1
t=0; iter=0;
[ funfI] = fun (x_fI, y_fI, z_fI );
[ funfII ] = fun (x_fII, y_fII, z_fII );
[ funfIII ] = fun (x_fIII, y_fIII, z_fIII );
[ funfIV ] = fun (x_fIV, y_fIV, z_fIV );
[ funfV ] = fun (x_fV, y_fV, z_fV );
[ funfVI ] = fun (x_fVI, y_fVI, z_fVI );

funfI=funfI(:,:,1:3);
funfII=funfII(:,:,1:3);
funfIII=funfIII(:,:,1:3);
funfIV=funfIV(:,:,1:3);
funfV=funfV(:,:,1:3);
funfVI=funfVI(:,:,1:3);

[div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
        div72(funfI(:,:,1:3),funfII(:,:,1:3),funfIII(:,:,1:3)...
        ,funfIV(:,:,1:3),funfV(:,:,1:3),funfVI(:,:,1:3),n,nn);
figure(1)
plot_cs11(n,nn,div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI)

%% BLOC 2

[mfunfI,dfunI]=fun6(x_fI,y_fI,z_fI); % FONCTION VECTORIELLE SUR LES 6 FACES
[mfunfII,dfunII]=fun6(x_fII,y_fII,z_fII);
[mfunfIII,dfunIII]=fun6(x_fIII,y_fIII,z_fIII);
[mfunfIV,dfunIV]=fun6(x_fIV,y_fIV,z_fIV);
[mfunfV,dfunV]=fun6(x_fV,y_fV,z_fV);
[mfunfVI,dfunVI]=fun6(x_fVI,y_fVI,z_fVI);

[mdiv_fI,mdiv_fII,mdiv_fIII,mdiv_fIV,mdiv_fV,mdiv_fVI]=...
    div72(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI,n,nn);

figure(2)
plot_cs11(n,nn,mdiv_fI,mdiv_fII,mdiv_fIII,mdiv_fIV,mdiv_fV,mdiv_fVI)
