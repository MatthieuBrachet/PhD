clc; clear all; close all;

global dxi deta dga
global radius
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;

N=32;
% initialisation
make_cs_grid(N);

nhs_max=50;
sym=0;
err=1;
k=0;


%% BLOC SANS SYMETRIES
[A,err_i] = compute_A(nhs_max);
res = solve_weights( A,err_i,k,err,sym );

%% BLOC AVEC SYMETRIES
[A,err_i] = compute_A_sym(nhs_max);
eps_w=solve_weights(A,err_i,k,err,1);

max(max(res-eps_w))./max(max(eps_w));

%% GRAPHES
figure(1)
subplot(121)
contourf(real(res))
colorbar

subplot(122)
contourf(imag(res))
colorbar
title('without symetries')

figure(2)
subplot(121)
contourf(real(eps_w))
colorbar

subplot(122)
contourf(imag(eps_w))
colorbar
title('with symetries')

fig_placier

