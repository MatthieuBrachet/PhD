clc; clear all; close all;

global dxi deta dga
global radius
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;

N=16;
% initialisation
make_cs_grid(N);

nhs_max=15; % combre d'harmoniques à corriger
%[A,err_i] = compute_A_sym(nhs_max);
[A,err_i] = compute_A(nhs_max);

sym=0; % utiliser la symétrie
err=1; % poids epsilon (1) ou poids epsilon+dga (0)
k=0; % ???
res = solve_weights( A,err_i,k,err,sym );
%res=real(res);

nhs=2*nhs_max; % nombre d'harmoniques sphériques à tester
weights=dxi*deta*(dga+res); % poids pour l'intégration
detail=0; % visualisation de l'erreur
sph4=0; % harmoniques a tester
err_i = view_int_sph( nhs,weights,detail,sph4 )

err = int_funs_fornberg( weights )

figure(1)
semilogy(err_i)
title('test SH')

figure(2)
semilogy(err)
title('test Fornberg')

figure(3)
surf(real(dga))
title('dga')

figure(4)
contourf(real(res))
title('add weath')

fig_placier
