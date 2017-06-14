clc; clear all; close all;
format long

global n nn
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global dxi corner
global opt_ftr scheme filtre

filtre='classic';
opt_ftr='redonnet10';
scheme='compact4';

test=0;
n=127;
mod101

[fun_I    ,int] = fun_quad(x_fI  ,y_fI  ,z_fI   ,test);
[fun_II   ,int] = fun_quad(x_fII ,y_fII ,z_fII  ,test);
[fun_III  ,int] = fun_quad(x_fIII,y_fIII,z_fIII ,test);
[fun_IV   ,int] = fun_quad(x_fIV ,y_fIV ,z_fIV  ,test);
[fun_V    ,int] = fun_quad(x_fV  ,y_fV  ,z_fV   ,test);
[fun_VI   ,int] = fun_quad(x_fVI ,y_fVI ,z_fVI  ,test);

Alpha=0:0.001:.5;

for i=1:length(Alpha)
    clc; corner=Alpha(i)
    
    [nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=nrm_alpha(fun_I,fun_II,fun_III,fun_IV,fun_V,fun_VI,n,nn,corner);
    err(i)=abs(nrmI-int/6);
end

figure(1)
plot(Alpha,err)