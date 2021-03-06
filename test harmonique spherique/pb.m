clc; clear all; close all

global n nn;
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI;
global opt_ftr;

%% *** Benchmarks data ****************************************************
NN=[10:10:100];
E=[];
for i=1:length(NN);
    clc; n=NN(i)
    nn=n+2;
    l=10;
    m=2;   % be carreful : [-l <= m <= l]
    %% *** Options ************************************************************
    opt_ftr=10;
    
    %% ************************************************************************
    mod_1b
    
    [h_fI]=fun(x_fI,y_fI,z_fI,m,l);
    [h_fII]=fun(x_fII,y_fII,z_fII,m,l);
    [h_fIII]=fun(x_fIII,y_fIII,z_fIII,m,l);
    [h_fIV]=fun(x_fIV,y_fIV,z_fIV,m,l);
    [h_fV]=fun(x_fV,y_fV,z_fV,m,l);
    [h_fVI]=fun(x_fVI,y_fVI,z_fVI,m,l);

     [grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
        gr(h_fI,h_fII,h_fIII,h_fIV,h_fV,h_fVI,n,nn);

     [grade_I] = gr_harm(m,l,x_fI,y_fI,z_fI);
     [grade_II] = gr_harm(m,l,x_fII,y_fII,z_fII);
     [grade_III] = gr_harm(m,l,x_fIII,y_fIII,z_fIII);
     [grade_IV] = gr_harm(m,l,x_fIV,y_fIV,z_fIV);
     [grade_V] = gr_harm(m,l,x_fV,y_fV,z_fV);
     [grade_VI] = gr_harm(m,l,x_fVI,y_fVI,z_fVI);
 
    %% * CALCUL DE L ERREUR ***************************************************

    err_fI=max(max(max(abs(grade_I-grad_I))));
    err_fII=max(max(max(abs(grade_II-grad_II))));
    err_fIII=max(max(max(abs(grade_III-grad_III))));
    err_fIV=max(max(max(abs(grade_IV-grad_IV))));
    err_fV=max(max(max(abs(grade_V-grad_V))));
    err_fVI=max(max(max(abs(grade_VI-grad_VI))));
    
    e=max([err_fI; err_fII; err_fIII; err_fIV; err_fV; err_fVI]);
    E=[E e];
end



figure(1)
semilogy(NN,E)
title('erreur')