clear all;clc; close all
global n nn;
global radius;
global dxi deta dga;
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;

N=16;
make_cs_grid(N);

% COMPUTATION OF THE WEIGHTS

%%% with weights with some properties of symmetry
% weights=create_weights(nn,4);

%%% with the basic formula
weights=dxi*deta*dga;

%%% with the method 1
% [eps_w,opt_val]=eps_weights(nhs_max);
% weights=dxi*deta*(dga+eps_w);

%%% with the method 2
% nhs_max=0;
% k=120;
% [A,err_i]=compute_A_sym(nhs_max);
% eps_w=solve_weights(A,err_i,k,1,1);
% %eps_w=zeros(nn,nn);
% weights=dxi*deta*(dga+eps_w);

% SHOW INTEGRAL VALUES FOR SPH OF ORDER nhs
nhs = 15;
err_i=view_int_sph(nhs,weights,1,0);
y=int_funs_fornberg(weights);


